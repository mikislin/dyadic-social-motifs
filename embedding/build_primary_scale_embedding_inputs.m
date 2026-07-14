function Input = build_primary_scale_embedding_inputs(repoRoot, primaryScales, anchorManifest, eventSummary, transformAudit, featureDict, params)
%BUILD_PRIMARY_SCALE_EMBEDDING_INPUTS Materialize run-07 summary matrices.
%
% This helper consumes the run_06 primary scale-specific anchor manifest.
% Matrix construction uses feature MAT files, frame masks, run_06 transform
% stats, and predefined configuration. Provenance labels are not used by this
% loader; the embedding fitter strips condition/cohort-style labels before
% saving row metadata and retains arena only for post-fit audit overlays.

repoRoot = string(repoRoot);
primaryScales = sortrows(primaryScales, 'chunk_sec');

stats = local_stats_from_transform_audit(transformAudit);
scoreFeatures = string(featureDict.Name(featureDict.ClusteringCandidate == 1));
assert(~isempty(scoreFeatures), 'build_primary_scale_embedding_inputs:NoFeatures', ...
    'Feature dictionary does not mark any clustering-candidate features.');

nScale = height(primaryScales);
Scale = repmat(struct( ...
    'scale_index', [], ...
    'primary_scale_rank', [], ...
    'chunk_sec', [], ...
    'hierarchical_role', "", ...
    'rowMeta', table(), ...
    'X', [], ...
    'dimMeta', table(), ...
    'summaryAudit', table()), nScale, 1);

rowManifest = table();
featureRows = table();
matrixAudit = table();
nextRow = 1;

for s = 1:nScale
    scaleIndex = double(primaryScales.scale_index(s));
    rows = anchorManifest.scale_index == scaleIndex;
    A = sortrows(anchorManifest(rows, :), {'raw_index', 'anchor_time_s'});
    assert(~isempty(A), 'build_primary_scale_embedding_inputs:NoAnchorsForScale', ...
        'No run_07 anchors were available for scale_index=%g.', scaleIndex);

    [Sc, ChunkSet] = local_materialize_scale(repoRoot, A, stats, params);
    [Xsummary, dimMeta, summaryAudit] = summarize_multiresolution_chunks(Sc, ChunkSet, ...
        'featureNames', scoreFeatures, ...
        'nTemporalBins', params.summary_temporal_bins, ...
        'nDctCoeffs', params.summary_dct_coeffs, ...
        'includeDct', params.summary_include_dct, ...
        'useScaledFeatures', params.summary_use_scaled_features);
    dimMeta = local_add_feature_metadata(dimMeta, featureDict);

    [Xevent, eventDimMeta] = local_event_feature_matrix(A, eventSummary, params);
    if ~isempty(eventDimMeta)
        eventDimMeta.summary_dim_index = eventDimMeta.summary_dim_index + size(Xsummary, 2);
    end

    X = [double(Xsummary), double(Xevent)];
    dimMeta = [dimMeta; eventDimMeta]; %#ok<AGROW>
    dimMeta.embedding_feature_index = (1:height(dimMeta))';
    dimMeta = movevars(dimMeta, 'embedding_feature_index', 'Before', 1);

    rowIds = nextRow:(nextRow + height(A) - 1);
    nextRow = nextRow + height(A);
    A.embedding_row_id = rowIds(:);
    A.run07_matrix_role = repmat("scale_specific_primary_anchor_embedding_row", height(A), 1);
    A.labels_used_for_embedding_matrix = repmat("none", height(A), 1);
    A.arena_used_for_embedding_matrix = false(height(A), 1);
    A.condition_used_for_embedding_matrix = false(height(A), 1);

    Scale(s).scale_index = scaleIndex;
    Scale(s).primary_scale_rank = local_table_double(primaryScales, 'primary_scale_rank', s, s);
    Scale(s).chunk_sec = double(primaryScales.chunk_sec(s));
    Scale(s).hierarchical_role = local_table_string(primaryScales, 'hierarchical_role', s, "");
    Scale(s).rowMeta = A;
    Scale(s).X = X;
    Scale(s).dimMeta = dimMeta;
    Scale(s).summaryAudit = summaryAudit;

    dimOut = dimMeta;
    dimOut.scale_index = repmat(scaleIndex, height(dimOut), 1);
    dimOut.chunk_sec = repmat(double(primaryScales.chunk_sec(s)), height(dimOut), 1);
    featureRows = [featureRows; dimOut]; %#ok<AGROW>
    rowManifest = [rowManifest; A]; %#ok<AGROW>
    matrixAudit = [matrixAudit; local_matrix_audit(Scale(s), X, Xsummary, Xevent)]; %#ok<AGROW>
end

Input = struct();
Input.scale = Scale;
Input.rowManifest = sortrows(rowManifest, 'embedding_row_id');
Input.featureDictionary = sortrows(featureRows, {'scale_index', 'embedding_feature_index'});
Input.matrixAudit = matrixAudit;
Input.labels_used_for_matrix = "none";
Input.arena_used_for_matrix = false;
Input.condition_used_for_matrix = false;
end

function [Sc, ChunkSet] = local_materialize_scale(repoRoot, A, stats, params)
n = height(A);
L = round(A.chunk_frames(1));
assert(all(round(A.chunk_frames) == L), ...
    'build_primary_scale_embedding_inputs:MixedChunkFrames', ...
    'A scale-specific anchor table must have one chunk_frames value.');

firstDyad = local_load_dyad(repoRoot, string(A.feature_output_file(1)));
firstSeq = prepare_dyad_timeseries(firstDyad, 'stats', stats);
local_assert_fps(firstSeq, params);
D = size(firstSeq.Xscaled, 2);

Sc = struct();
Sc.chunkSec = double(A.chunk_sec(1));
Sc.nFrames = L;
Sc.X = nan(n, L, D, 'single');
Sc.Xraw = nan(n, L, D, 'single');
Sc.valid = false(n, L);
Sc.meta = A;
Sc.sessionIndex = A.session_index;
Sc.obsNames = firstSeq.obsNames;
Sc.channelMeta = firstSeq.channelMeta;
Sc.featureNames = firstSeq.featureNames;
Sc.featureMeta = firstSeq.featureMeta;

paths = unique(string(A.feature_output_file), 'stable')';
for p = paths
    idxRows = find(string(A.feature_output_file) == p);
    dyad = local_load_dyad(repoRoot, p);
    Seq = prepare_dyad_timeseries(dyad, 'stats', stats);
    local_assert_fps(Seq, params);
    assert(size(Seq.Xscaled, 2) == D, ...
        'build_primary_scale_embedding_inputs:ChannelCountMismatch', ...
        'Feature channel count changed across sessions.');
    for rr = idxRows(:)'
        st = max(1, round(A.start_frame(rr)));
        en = min(numel(Seq.time), round(A.stop_frame(rr)));
        if st > en
            continue
        end
        chunkIdx = st:en;
        if numel(chunkIdx) ~= L
            continue
        end
        Sc.X(rr, :, :) = single(Seq.Xscaled(chunkIdx, :));
        Sc.Xraw(rr, :, :) = single(Seq.X(chunkIdx, :));
        Sc.valid(rr, :) = logical(Seq.validMask(chunkIdx))';
    end
end

ChunkSet = struct();
ChunkSet.obsNames = firstSeq.obsNames;
ChunkSet.channelMeta = firstSeq.channelMeta;
ChunkSet.featureNames = firstSeq.featureNames;
ChunkSet.featureMeta = firstSeq.featureMeta;
ChunkSet.anchorMode = char(params.anchor_mode);
ChunkSet.scale = Sc;
end

function dyad = local_load_dyad(repoRoot, featureOutputFile)
featurePath = resolve_repo_path(repoRoot, featureOutputFile);
S = load(featurePath, 'dyad', 'status');
assert(isfield(S, 'dyad'), 'build_primary_scale_embedding_inputs:MissingDyad', ...
    'Missing dyad in feature file: %s', featurePath);
if isfield(S, 'status')
    assert(string(S.status) == "success", ...
        'build_primary_scale_embedding_inputs:BadFeatureStatus', ...
        'Feature file does not report success: %s', featurePath);
end
dyad = S.dyad;
end

function stats = local_stats_from_transform_audit(T)
required = ["robust_center", "robust_scale", "impute_value"];
missing = setdiff(required, string(T.Properties.VariableNames));
assert(isempty(missing), 'build_primary_scale_embedding_inputs:BadTransformAudit', ...
    'Transform audit missing required columns: %s', strjoin(missing, ', '));
if ismember('channel_index', T.Properties.VariableNames)
    T = sortrows(T, 'channel_index');
end
stats = struct();
stats.center = double(T.robust_center(:))';
stats.scale = double(T.robust_scale(:))';
stats.impute = double(T.impute_value(:))';
if ismember('n_fit_rows_total', T.Properties.VariableNames)
    stats.n_fit_rows = double(T.n_fit_rows_total(1));
else
    stats.n_fit_rows = NaN;
end
if ismember('stats_fit_scope', T.Properties.VariableNames)
    stats.fit_scope = string(T.stats_fit_scope(1));
else
    stats.fit_scope = "run06_condition_blind_transform_audit";
end
end

function [Xevent, dimMeta] = local_event_feature_matrix(A, eventSummary, params)
Xevent = zeros(height(A), 0);
dimMeta = table();
if ~logical(params.include_event_summary_features) || isempty(eventSummary)
    return
end

eventNames = params.event_summary_feature_names;
if isstring(eventNames) && numel(eventNames) == 1 && eventNames == "all"
    eventNames = local_default_event_feature_names();
end
eventNames = eventNames(ismember(eventNames, string(eventSummary.Properties.VariableNames)));
if isempty(eventNames)
    return
end

E = local_match_event_rows(A, eventSummary);
Xevent = nan(height(A), numel(eventNames));
for j = 1:numel(eventNames)
    Xevent(:, j) = double(E.(eventNames(j)));
end

n = numel(eventNames);
dimMeta = table();
dimMeta.summary_dim_index = (1:n)';
dimMeta.scale_index = repmat(double(A.scale_index(1)), n, 1);
dimMeta.chunk_sec = repmat(double(A.chunk_sec(1)), n, 1);
dimMeta.obs_name = eventNames(:);
dimMeta.base_feature = eventNames(:);
dimMeta.channel_type = repmat("event_summary", n, 1);
dimMeta.summary_kind = repmat("primary_event_summary", n, 1);
dimMeta.temporal_bin = nan(n, 1);
dimMeta.dct_coefficient = nan(n, 1);
dimMeta.source_tensor = repmat("primary_chunk_event_summary_audit", n, 1);
dimMeta.mean_frame_valid_fraction_for_channel = repmat(mean(E.event_valid_fraction, 'omitnan'), n, 1);
dimMeta.uses_frame_mask = true(n, 1);
dimMeta.labels_used_for_summary = repmat("none", n, 1);
dimMeta.feature_family = repmat("event_summary", n, 1);
dimMeta.unit = local_event_units(eventNames);
end

function E = local_match_event_rows(A, eventSummary)
if ismember('primary_anchor_global_id', A.Properties.VariableNames) && ...
        ismember('primary_anchor_global_id', eventSummary.Properties.VariableNames)
    [tf, loc] = ismember(A.primary_anchor_global_id, eventSummary.primary_anchor_global_id);
else
    keyA = local_anchor_key(A.scale_index, A.raw_index, A.anchor_frame);
    keyE = local_anchor_key(eventSummary.scale_index, eventSummary.raw_index, eventSummary.anchor_frame);
    [tf, loc] = ismember(keyA, keyE);
end
assert(all(tf), 'build_primary_scale_embedding_inputs:EventRowsMissing', ...
    'Primary event summary rows could not be matched to selected anchors.');
E = eventSummary(loc, :);
end

function key = local_anchor_key(scaleIndex, rawIndex, anchorFrame)
key = string(round(double(scaleIndex(:)))) + "_" + ...
    string(round(double(rawIndex(:)))) + "_" + ...
    string(round(double(anchorFrame(:))));
end

function names = local_default_event_feature_names()
names = ["event_valid_fraction", ...
    "contact_dwell_fraction", "contact_transition_count", "first_contact_latency_s", ...
    "mutual_approach_dwell_fraction", "withdrawal_dwell_fraction", ...
    "approach_withdraw_transition_count", "asym_positive_dwell_fraction", ...
    "asym_negative_dwell_fraction", "role_asymmetry_bias_mean", ...
    "centroid_distance_mean_mm", "centroid_distance_min_mm", ...
    "centroid_distance_delta_mm", "radial_speed_mean_mm_s", ...
    "radial_speed_sign_change_count", "approach1_positive_fraction", ...
    "approach2_positive_fraction", "role_approach_imbalance", ...
    "mutual_facing_mean", "heading_abs_change_deg", "heading_large_turn_count"];
end

function units = local_event_units(names)
units = repmat("unitless", numel(names), 1);
for i = 1:numel(names)
    if endsWith(names(i), "_mm")
        units(i) = "mm";
    elseif endsWith(names(i), "_mm_s")
        units(i) = "mm_per_s";
    elseif contains(names(i), "latency") || contains(names(i), "duration")
        units(i) = "s";
    elseif contains(names(i), "deg")
        units(i) = "degrees";
    elseif contains(names(i), "count")
        units(i) = "count";
    elseif contains(names(i), "fraction")
        units(i) = "fraction";
    end
end
end

function dimMeta = local_add_feature_metadata(dimMeta, featureDict)
dimMeta.feature_family = repmat("unknown", height(dimMeta), 1);
dimMeta.unit = repmat("", height(dimMeta), 1);
for i = 1:height(dimMeta)
    row = find(string(featureDict.Name) == string(dimMeta.base_feature(i)), 1);
    if ~isempty(row)
        dimMeta.feature_family(i) = string(featureDict.Family(row));
        dimMeta.unit(i) = string(featureDict.Unit(row));
    end
end
end

function audit = local_matrix_audit(S, X, Xsummary, Xevent)
audit = table();
audit.scale_index = S.scale_index;
audit.chunk_sec = S.chunk_sec;
audit.hierarchical_role = string(S.hierarchical_role);
audit.n_rows = size(X, 1);
audit.n_multiresolution_dims = size(Xsummary, 2);
audit.n_event_summary_dims = size(Xevent, 2);
audit.n_total_dims = size(X, 2);
audit.finite_value_fraction = mean(isfinite(X), 'all');
audit.labels_used_for_matrix = "none";
audit.arena_used_for_matrix = false;
audit.condition_used_for_matrix = false;
audit.matrix_rule = "condition_blind_multiresolution_summaries_plus_optional_event_summary_features";
end

function local_assert_fps(Seq, params)
assert(abs(double(Seq.fps) - double(params.fps)) <= 1e-10, ...
    'build_primary_scale_embedding_inputs:FpsMismatch', ...
    'Loaded feature session fps %.15g does not match run_07 config fps %.15g.', ...
    double(Seq.fps), double(params.fps));
end

function value = local_table_string(T, varName, rowIdx, defaultValue)
if ismember(varName, T.Properties.VariableNames)
    value = string(T.(varName)(rowIdx));
else
    value = string(defaultValue);
end
end

function value = local_table_double(T, varName, rowIdx, defaultValue)
if ismember(varName, T.Properties.VariableNames)
    value = double(T.(varName)(rowIdx));
else
    value = defaultValue;
end
end
