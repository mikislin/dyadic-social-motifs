function ScaleShard = build_scale_summary_shard_from_anchor_manifest(repoRoot, sessionTable, anchorTable, scaleRow, stats, varargin)
%BUILD_SCALE_SUMMARY_SHARD_FROM_ANCHOR_MANIFEST Stream one run_06 scale shard.
%
% This helper preserves the run_06 condition-blind representation logic while
% avoiding an all-scale dense ChunkSet. It loads one feature session at a time,
% materializes only that session's chunks for one temporal scale, summarizes
% them, and releases the frame-level tensor before moving to the next session.

p = inputParser;
p.addParameter('fps', 80, @(x)isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('anchorMode', "past", @(x)any(string(x) == ["center","past"]));
p.addParameter('minValidFrac', 0.85, @(x)isnumeric(x) && isscalar(x));
p.addParameter('minFeatureFiniteFrac', 0.95, @(x)isnumeric(x) && isscalar(x));
p.addParameter('featureNames', strings(0,1), @(x)isstring(x) || iscell(x));
p.addParameter('summaryTemporalBins', 6, @(x)isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('summaryDctCoeffs', 4, @(x)isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('summaryIncludeDct', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('summaryUseScaledFeatures', true, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});
P = p.Results;
P.anchorMode = lower(string(P.anchorMode));

assert(istable(sessionTable) && ~isempty(sessionTable), ...
    'build_scale_summary_shard:EmptySessionTable', ...
    'sessionTable must contain at least one session.');
assert(istable(anchorTable) && ~isempty(anchorTable), ...
    'build_scale_summary_shard:EmptyAnchorTable', ...
    'anchorTable must contain at least one anchor.');
assert(istable(scaleRow) && height(scaleRow) == 1, ...
    'build_scale_summary_shard:BadScaleRow', ...
    'scaleRow must be a one-row scale table.');

requiredAnchor = ["anchor_id", "feature_row_index", "session_index", "anchor_frame", "anchor_time_s"];
missingAnchor = setdiff(requiredAnchor, string(anchorTable.Properties.VariableNames));
assert(isempty(missingAnchor), 'build_scale_summary_shard:MissingAnchorColumn', ...
    'anchorTable missing required columns: %s', strjoin(missingAnchor, ', '));

scaleIndex = local_table_double(scaleRow, 'scale_index', 1, 1);
chunkSec = local_table_double(scaleRow, 'chunk_sec', 1, NaN);
assert(isfinite(chunkSec) && chunkSec > 0, ...
    'build_scale_summary_shard:BadChunkSec', ...
    'scaleRow.chunk_sec must be positive.');
chunkFrames = max(1, round(chunkSec * P.fps));
[leftFrames, rightFrames] = local_chunk_geometry(chunkFrames, P.anchorMode);

meta = local_scale_meta(anchorTable, sessionTable, scaleRow, chunkFrames, P.anchorMode);
nAnchors = height(anchorTable);

Xsummary = [];
dimMeta = table();
dimensionAudit = table();
centerScaled = [];
sessionRows = table();

obsNames = {};
channelMeta = table();
featureNames = {};
featureMeta = table();
nObs = NaN;

sessionIds = unique(anchorTable.feature_row_index(:), 'stable')';
for sess = sessionIds
    rows = find(anchorTable.feature_row_index == sess);
    if isempty(rows)
        continue
    end

    dyad = local_load_dyad(sessionTable, sess, repoRoot);
    Seq = prepare_dyad_timeseries(dyad, 'stats', stats);
    local_assert_seq_fps(Seq, P.fps);

    if isempty(obsNames)
        obsNames = Seq.obsNames;
        channelMeta = Seq.channelMeta;
        featureNames = Seq.featureNames;
        featureMeta = Seq.featureMeta;
        nObs = size(Seq.Xscaled, 2);
        centerScaled = nan(nAnchors, nObs, 'single');
    else
        assert(size(Seq.Xscaled, 2) == nObs, ...
            'build_scale_summary_shard:ChannelCountMismatch', ...
            'Session %d has a different transformed channel count.', sess);
        assert(all(string(Seq.channelMeta.ObsName) == string(channelMeta.ObsName)), ...
            'build_scale_summary_shard:ChannelOrderMismatch', ...
            'Session %d transformed channel order differs from the first session.', sess);
    end

    nSessAnchors = numel(rows);
    Xscaled = nan(nSessAnchors, chunkFrames, nObs, 'single');
    Xraw = nan(nSessAnchors, chunkFrames, nObs, 'single');
    V = false(nSessAnchors, chunkFrames);

    for ii = 1:nSessAnchors
        rr = rows(ii);
        anchorFrame = round(anchorTable.anchor_frame(rr));
        idx = (anchorFrame - leftFrames):(anchorFrame + rightFrames);
        if idx(1) < 1 || idx(end) > numel(Seq.time)
            continue
        end

        Xraw(ii, :, :) = single(Seq.X(idx, :));
        Xscaled(ii, :, :) = single(Seq.Xscaled(idx, :));
        V(ii, :) = logical(Seq.validMask(idx));

        validFrac = mean(V(ii, :));
        finiteFrac = mean(isfinite(Seq.X(idx, :)), 'all');
        meta.start_frame(rr) = idx(1);
        meta.stop_frame(rr) = idx(end);
        meta.start_time_s(rr) = Seq.time(idx(1));
        meta.stop_time_s(rr) = Seq.time(idx(end));
        meta.valid_frac(rr) = validFrac;
        meta.feature_finite_frac(rr) = finiteFrac;
        meta.chunk_is_valid(rr) = validFrac >= P.minValidFrac && finiteFrac >= P.minFeatureFiniteFrac;

        if P.anchorMode == "past"
            centerIdx = chunkFrames;
        else
            centerIdx = floor((chunkFrames - 1) / 2) + 1;
        end
        centerScaled(rr, :) = single(Seq.Xscaled(idx(centerIdx), :));
    end

    ScSess = struct();
    ScSess.chunkSec = chunkSec;
    ScSess.nFrames = chunkFrames;
    ScSess.X = Xscaled;
    ScSess.Xraw = Xraw;
    ScSess.valid = V;
    ScSess.meta = meta(rows, :);

    ChunkStub = struct();
    ChunkStub.channelMeta = channelMeta;
    ChunkStub.featureNames = featureNames;

    [XsessSummary, oneDimMeta, oneAudit] = summarize_multiresolution_chunks(ScSess, ChunkStub, ...
        'featureNames', string(P.featureNames), ...
        'nTemporalBins', P.summaryTemporalBins, ...
        'nDctCoeffs', P.summaryDctCoeffs, ...
        'includeDct', P.summaryIncludeDct, ...
        'useScaledFeatures', P.summaryUseScaledFeatures);

    if isempty(Xsummary)
        Xsummary = nan(nAnchors, size(XsessSummary, 2), 'single');
        dimMeta = oneDimMeta;
        dimensionAudit = oneAudit;
    else
        assert(size(XsessSummary, 2) == size(Xsummary, 2), ...
            'build_scale_summary_shard:SummaryDimMismatch', ...
            'Session %d summary dimension count differs from previous sessions.', sess);
    end
    Xsummary(rows, :) = single(XsessSummary);

    oneSession = table();
    oneSession.scale_index = scaleIndex;
    oneSession.chunk_sec = chunkSec;
    oneSession.feature_row_index = sess;
    oneSession.session_index = local_table_double(sessionTable, 'session_index', sess, sess);
    oneSession.raw_index = local_table_double(sessionTable, 'raw_index', sess, NaN);
    oneSession.session_id = local_table_string(sessionTable, 'session_id', sess, "");
    oneSession.n_anchor_rows = nSessAnchors;
    oneSession.n_valid_chunks = nnz(meta.chunk_is_valid(rows));
    oneSession.mean_valid_frac = mean(meta.valid_frac(rows), 'omitnan');
    oneSession.labels_used_for_summary = "none";
    sessionRows = [sessionRows; oneSession]; %#ok<AGROW>
end

assert(~isempty(Xsummary), 'build_scale_summary_shard:NoSummariesBuilt', ...
    'No scale summaries were built for scale %.4f.', chunkSec);

meta = local_finalize_scale_meta(meta, P);
dimensionAudit.n_chunks = nAnchors;
dimensionAudit.n_pca_fit_chunks = NaN;
dimensionAudit.n_pca_fit_chunks_used = NaN;
dimensionAudit.summary_source = "streamed_session_scale_shard";
dimensionAudit.labels_used_for_summary = "none";

scale = struct();
scale.chunkSec = chunkSec;
scale.nFrames = chunkFrames;
scale.X = [];
scale.Xraw = [];
scale.valid = [];
scale.meta = meta;
scale.sessionIndex = anchorTable.session_index;
scale.obsNames = obsNames;
scale.channelMeta = channelMeta;
scale.featureNames = featureNames;
scale.featureMeta = featureMeta;
scale.Xsummary = Xsummary;
scale.summaryDimMeta = dimMeta;
scale.dimensionAudit = dimensionAudit;
scale.centerScaled = centerScaled;
scale.summaryShardMode = "streamed_session_scale_summary";

ScaleShard = struct();
ScaleShard.shard_version = "run06_scale_summary_shard_v1";
ScaleShard.scale_index = scaleIndex;
ScaleShard.chunk_sec = chunkSec;
ScaleShard.chunk_frames = chunkFrames;
ScaleShard.n_anchor_rows = nAnchors;
ScaleShard.summary_temporal_bins = P.summaryTemporalBins;
ScaleShard.summary_dct_coeffs = P.summaryDctCoeffs;
ScaleShard.summary_include_dct = logical(P.summaryIncludeDct);
ScaleShard.summary_use_scaled_features = logical(P.summaryUseScaledFeatures);
ScaleShard.scale = scale;
ScaleShard.session_summary = sessionRows;
ScaleShard.created_at = string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
ScaleShard.labels_used_for_shard = "none";
ScaleShard.condition_used_for_shard = false;
ScaleShard.arena_used_for_shard = false;
end

function meta = local_scale_meta(anchorTable, sessionTable, scaleRow, L, anchorMode)
n = height(anchorTable);
meta = anchorTable;
meta.scale_index = repmat(local_table_double(scaleRow, 'scale_index', 1, 1), n, 1);
meta.chunk_sec = repmat(local_table_double(scaleRow, 'chunk_sec', 1, NaN), n, 1);
meta.chunk_frames = repmat(L, n, 1);
meta.anchor_mode = repmat(string(anchorMode), n, 1);
meta.temporal_band = repmat(local_table_string(scaleRow, 'temporal_band', 1, local_scale_band_label(meta.chunk_sec(1))), n, 1);
meta.hierarchical_role = repmat(local_table_string(scaleRow, 'hierarchical_role', 1, meta.temporal_band(1)), n, 1);
meta.chunk_id = ((1:n)' + (meta.scale_index - 1) * n);
meta.start_frame = nan(n, 1);
meta.stop_frame = nan(n, 1);
meta.start_time_s = nan(n, 1);
meta.stop_time_s = nan(n, 1);
meta.valid_frac = nan(n, 1);
meta.feature_finite_frac = nan(n, 1);
meta.chunk_is_valid = false(n, 1);
meta.condition_role = repmat("provenance_only_not_used_for_chunking_or_scale_selection", n, 1);
meta.feature_output_file = strings(n, 1);
meta.qc_set = strings(n, 1);
meta.valid_frame_fraction = nan(n, 1);
meta.badframe_mask_fraction = nan(n, 1);

for i = 1:n
    r = meta.feature_row_index(i);
    meta.feature_output_file(i) = local_table_string(sessionTable, 'feature_output_file', r, "");
    meta.qc_set(i) = local_table_string(sessionTable, 'qc_set', r, "");
    meta.valid_frame_fraction(i) = local_table_double(sessionTable, 'valid_frame_fraction', r, NaN);
    meta.badframe_mask_fraction(i) = local_table_double(sessionTable, 'badframe_mask_fraction', r, NaN);
end
end

function meta = local_finalize_scale_meta(meta, P)
meta.inclusion_flag = meta.chunk_is_valid;
meta.inclusion_rule = repmat( ...
    "predefined_condition_blind_anchor_and_chunk_validity_thresholds", ...
    height(meta), 1);
meta.validity_threshold = repmat(P.minValidFrac, height(meta), 1);
meta.feature_finite_threshold = repmat(P.minFeatureFiniteFrac, height(meta), 1);
end

function dyad = local_load_dyad(sessionTable, rowIdx, repoRoot)
pathText = string(sessionTable.feature_output_file(rowIdx));
featurePath = resolve_repo_path(repoRoot, pathText);
S = load(featurePath, 'dyad', 'status');
assert(isfield(S, 'dyad'), 'build_scale_summary_shard:MissingDyad', ...
    'Missing dyad in %s.', featurePath);
if isfield(S, 'status')
    assert(string(S.status) == "success", ...
        'build_scale_summary_shard:BadFeatureStatus', ...
        'Feature MAT %s does not report success.', featurePath);
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

function local_assert_seq_fps(Seq, expectedFps)
assert(abs(double(Seq.fps) - double(expectedFps)) <= 1e-10, ...
    'build_scale_summary_shard:FpsMismatch', ...
    'Loaded feature session fps %.15g does not match run_06 config fps %.15g.', ...
    double(Seq.fps), double(expectedFps));
end

function label = local_scale_band_label(sec)
if sec < 0.8
    label = "micro";
elseif sec < 2.5
    label = "motif";
else
    label = "context";
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
