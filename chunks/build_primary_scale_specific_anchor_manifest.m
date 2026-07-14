function [anchorManifest, bankSummary] = build_primary_scale_specific_anchor_manifest(repoRoot, sessionTable, primaryScales, stats, varargin)
%BUILD_PRIMARY_SCALE_SPECIFIC_ANCHOR_MANIFEST Define final primary anchors.
%
% This helper builds a scale-specific anchor manifest for the stability- and
% dimension-supported run_06 scales. It does not use condition, cohort, arena,
% genotype, drug, or outcome labels for eligibility or materialization.

p = inputParser;
p.addParameter('strideSec', 0.25, @(x)isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('anchorMode', "past", @(x)any(string(x) == ["center","past"]));
p.addParameter('minValidFrac', 0.85, @(x)isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
p.addParameter('minFeatureFiniteFrac', 0.95, @(x)isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
p.addParameter('requireAnchorFrameValid', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('maxAnchorsPerScale', 1000, @(x)isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('minAnchorsPerSession', 2, @(x)isnumeric(x) && isscalar(x) && x >= 0);
p.parse(varargin{:});
P = p.Results;
P.anchorMode = string(P.anchorMode);
P.requireAnchorFrameValid = logical(P.requireAnchorFrameValid);

anchorManifest = local_empty_anchor_manifest();
bankSummary = local_empty_bank_summary();
if isempty(primaryScales)
    return
end

requiredScaleCols = ["scale_index", "chunk_sec"];
missingScaleCols = setdiff(requiredScaleCols, string(primaryScales.Properties.VariableNames));
assert(isempty(missingScaleCols), 'build_primary_scale_specific_anchor_manifest:MissingScaleColumn', ...
    'primaryScales is missing required columns: %s', strjoin(missingScaleCols, ', '));

scaleSec = double(primaryScales.chunk_sec(:));
nScale = numel(scaleSec);
candidatesByScale = cell(nScale, 1);
for s = 1:nScale
    candidatesByScale{s} = local_empty_candidate_table();
end

for i = 1:height(sessionTable)
    dyad = local_load_dyad(sessionTable, i, string(repoRoot));
    Seq = prepare_dyad_timeseries(dyad, 'stats', stats);
    for s = 1:nScale
        A = find_condition_blind_chunk_anchors(Seq, scaleSec(s), ...
            'strideSec', P.strideSec, ...
            'anchorMode', P.anchorMode, ...
            'minValidFrac', P.minValidFrac, ...
            'minFeatureFiniteFrac', P.minFeatureFiniteFrac, ...
            'requireAnchorFrameValid', P.requireAnchorFrameValid);
        if isempty(A)
            continue
        end
        L = max(1, round(scaleSec(s) * Seq.fps));
        [leftFrames, rightFrames] = local_chunk_geometry(L, P.anchorMode);
        C = table();
        C.scale_index = repmat(double(primaryScales.scale_index(s)), height(A), 1);
        C.primary_scale_rank = repmat(s, height(A), 1);
        C.chunk_sec = repmat(scaleSec(s), height(A), 1);
        C.chunk_frames = repmat(L, height(A), 1);
        C.feature_row_index = repmat(i, height(A), 1);
        C.session_index = repmat(double(sessionTable.session_index(i)), height(A), 1);
        C.raw_index = repmat(local_table_double(sessionTable, 'raw_index', i, i), height(A), 1);
        C.anchor_frame = A.anchor_frame;
        C.anchor_time_s = A.anchor_time_s;
        C.start_frame = A.anchor_frame - leftFrames;
        C.stop_frame = A.anchor_frame + rightFrames;
        C.start_time_s = Seq.time(C.start_frame);
        C.stop_time_s = Seq.time(C.stop_frame);
        C.valid_frac = A.min_scale_valid_frac;
        C.feature_finite_frac = A.min_scale_feature_finite_frac;
        C.anchor_frame_valid = A.anchor_frame_valid;
        C.chunk_is_valid = C.valid_frac >= P.minValidFrac & C.feature_finite_frac >= P.minFeatureFiniteFrac;
        candidatesByScale{s} = [candidatesByScale{s}; C]; %#ok<AGROW>
    end
end

for s = 1:nScale
    C = candidatesByScale{s};
    selected = select_condition_blind_anchor_subset(C, P.maxAnchorsPerScale, ...
        'minAnchorsPerSession', P.minAnchorsPerSession);
    selected = local_add_primary_scale_metadata(selected, primaryScales, s);
    selected = local_add_session_provenance(selected, sessionTable);
    selected.primary_scale_specific_anchor_id = (1:height(selected))';
    selected.anchor_selection_scope = repmat("scale_specific_primary_bank", height(selected), 1);
    selected.labels_used_for_primary_anchor_selection = repmat("none", height(selected), 1);
    selected.arena_used_for_primary_anchor_selection = false(height(selected), 1);
    selected.condition_used_for_primary_anchor_selection = false(height(selected), 1);
    selected.primary_anchor_rule = repmat( ...
        "scale_specific_frame_mask_feature_time_even_no_condition_labels", height(selected), 1);
    anchorManifest = [anchorManifest; selected]; %#ok<AGROW>

    one = table();
    one.scale_index = double(primaryScales.scale_index(s));
    one.chunk_sec = scaleSec(s);
    one.hierarchical_role = local_scale_string(primaryScales, 'hierarchical_role', s, "");
    one.rank_within_role = local_scale_double(primaryScales, 'rank_within_role', s, NaN);
    one.stable_and_dimension_supported = local_scale_logical(primaryScales, 'stable_and_dimension_supported', s, true);
    one.n_scale_specific_candidate_anchors = height(C);
    one.n_primary_scale_specific_anchors = height(selected);
    one.n_sessions_with_candidates = numel(unique(C.session_index));
    one.n_sessions_with_primary_anchors = numel(unique(selected.session_index));
    one.mean_selected_valid_frac = mean(selected.valid_frac, 'omitnan');
    one.min_selected_valid_frac = min(selected.valid_frac, [], 'omitnan');
    one.mean_selected_feature_finite_frac = mean(selected.feature_finite_frac, 'omitnan');
    one.min_selected_feature_finite_frac = min(selected.feature_finite_frac, [], 'omitnan');
    one.max_anchors_per_scale = P.maxAnchorsPerScale;
    one.primary_anchor_bank_rule = "scale_specific_primary_bank_after_stability_dimension_audit";
    one.labels_used_for_primary_anchor_bank = "none";
    one.arena_used_for_primary_anchor_bank = false;
    one.condition_used_for_primary_anchor_bank = false;
    bankSummary = [bankSummary; one]; %#ok<AGROW>
end

if ~isempty(anchorManifest)
    anchorManifest = sortrows(anchorManifest, {'scale_index', 'raw_index', 'anchor_time_s'});
    anchorManifest.primary_anchor_global_id = (1:height(anchorManifest))';
end
end

function C = local_empty_candidate_table()
C = table('Size', [0 16], ...
    'VariableTypes', {'double','double','double','double','double','double','double', ...
    'double','double','double','double','double','double','double','logical','logical'}, ...
    'VariableNames', {'scale_index','primary_scale_rank','chunk_sec','chunk_frames', ...
    'feature_row_index','session_index','raw_index','anchor_frame','anchor_time_s', ...
    'start_frame','stop_frame','start_time_s','stop_time_s','valid_frac', ...
    'anchor_frame_valid','chunk_is_valid'});
C.feature_finite_frac = zeros(0, 1);
C = movevars(C, 'feature_finite_frac', 'Before', 'anchor_frame_valid');
end

function T = local_empty_anchor_manifest()
T = table();
end

function T = local_empty_bank_summary()
T = table();
end

function T = local_add_primary_scale_metadata(T, primaryScales, s)
if isempty(T)
    return
end
T.hierarchical_role = repmat(local_scale_string(primaryScales, 'hierarchical_role', s, ""), height(T), 1);
T.rank_within_role = repmat(local_scale_double(primaryScales, 'rank_within_role', s, NaN), height(T), 1);
T.bootstrap_selection_frequency = repmat(local_scale_double(primaryScales, 'bootstrap_selection_frequency', s, NaN), height(T), 1);
T.passes_stability_threshold = repmat(local_scale_logical(primaryScales, 'passes_stability_threshold', s, true), height(T), 1);
T.passes_dimension_guard = repmat(local_scale_logical(primaryScales, 'passes_dimension_guard', s, true), height(T), 1);
T.stable_and_dimension_supported = repmat(local_scale_logical(primaryScales, 'stable_and_dimension_supported', s, true), height(T), 1);
end

function T = local_add_session_provenance(T, sessionTable)
if isempty(T)
    return
end
n = height(T);
T.session_id = strings(n, 1);
T.session_file = strings(n, 1);
T.arena = strings(n, 1);
T.arena_label = strings(n, 1);
T.cohort_id = strings(n, 1);
T.cohort_label = strings(n, 1);
T.condition_id = strings(n, 1);
T.condition_label = strings(n, 1);
T.qc_set = strings(n, 1);
T.feature_output_file = strings(n, 1);
for i = 1:n
    r = T.feature_row_index(i);
    T.session_id(i) = local_table_string(sessionTable, 'session_id', r, "");
    T.session_file(i) = local_table_string(sessionTable, 'session_file', r, "");
    T.arena(i) = local_table_string(sessionTable, 'arena', r, "");
    T.arena_label(i) = local_table_string(sessionTable, 'arena_label', r, "");
    T.cohort_id(i) = local_table_string(sessionTable, 'cohort_id', r, "");
    T.cohort_label(i) = local_table_string(sessionTable, 'cohort_label', r, "");
    T.condition_id(i) = local_table_string(sessionTable, 'condition_id', r, "");
    T.condition_label(i) = local_table_string(sessionTable, 'condition_label', r, "");
    T.qc_set(i) = local_table_string(sessionTable, 'qc_set', r, "");
    T.feature_output_file(i) = local_table_string(sessionTable, 'feature_output_file', r, "");
end
T.provenance_role = repmat("provenance_only_not_used_for_primary_anchor_bank", n, 1);
end

function dyad = local_load_dyad(sessionTable, rowIdx, repoRoot)
featurePath = resolve_repo_path(repoRoot, string(sessionTable.feature_output_file(rowIdx)));
S = load(featurePath, 'dyad', 'status');
assert(isfield(S, 'dyad'), 'build_primary_scale_specific_anchor_manifest:MissingDyad', ...
    'Missing dyad in feature file: %s', featurePath);
if isfield(S, 'status')
    assert(string(S.status) == "success", ...
        'build_primary_scale_specific_anchor_manifest:BadFeatureStatus', ...
        'Feature file does not report success: %s', featurePath);
end
dyad = S.dyad;
end

function [leftFrames, rightFrames] = local_chunk_geometry(L, anchorMode)
switch string(anchorMode)
    case "center"
        leftFrames = floor((L - 1) / 2);
        rightFrames = ceil((L - 1) / 2);
    case "past"
        leftFrames = L - 1;
        rightFrames = 0;
end
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

function value = local_scale_string(T, varName, rowIdx, defaultValue)
if ismember(varName, T.Properties.VariableNames)
    value = string(T.(varName)(rowIdx));
else
    value = string(defaultValue);
end
end

function value = local_scale_double(T, varName, rowIdx, defaultValue)
if ismember(varName, T.Properties.VariableNames)
    value = double(T.(varName)(rowIdx));
else
    value = defaultValue;
end
end

function value = local_scale_logical(T, varName, rowIdx, defaultValue)
if ismember(varName, T.Properties.VariableNames)
    value = logical(T.(varName)(rowIdx));
else
    value = logical(defaultValue);
end
end
