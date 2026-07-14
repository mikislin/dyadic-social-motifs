function [ChunkSet, anchorManifest] = build_chunkset_from_anchor_manifest(sessionTable, anchorTable, scaleTable, stats, opts)
%BUILD_CHUNKSET_FROM_ANCHOR_MANIFEST Materialize a multiscale ChunkSet.
%
% sessionTable rows define feature MAT files and provenance. anchorTable rows
% define a condition-blind anchor subset with feature_row_index values into
% sessionTable. Provenance labels are copied into output metadata only.

arguments
    sessionTable table
    anchorTable table
    scaleTable table
    stats struct
    opts.repoRoot (1,1) string = string(pwd)
    opts.anchorMode (1,1) string {mustBeMember(opts.anchorMode, ["center","past"])} = "center"
    opts.minValidFrac (1,1) double {mustBeGreaterThanOrEqual(opts.minValidFrac,0),mustBeLessThanOrEqual(opts.minValidFrac,1)} = 0.85
    opts.minFeatureFiniteFrac (1,1) double {mustBeGreaterThanOrEqual(opts.minFeatureFiniteFrac,0),mustBeLessThanOrEqual(opts.minFeatureFiniteFrac,1)} = 0.95
end

if isempty(anchorTable)
    error('build_chunkset_from_anchor_manifest:EmptyAnchorTable', ...
        'anchorTable must contain at least one materialized anchor.');
end

requiredAnchor = ["anchor_id", "feature_row_index", "session_index", "anchor_frame", "anchor_time_s"];
missingAnchor = setdiff(requiredAnchor, string(anchorTable.Properties.VariableNames));
assert(isempty(missingAnchor), 'build_chunkset_from_anchor_manifest:MissingAnchorColumn', ...
    'anchorTable missing required columns: %s', strjoin(missingAnchor, ', '));
assert(ismember('chunk_sec', scaleTable.Properties.VariableNames), ...
    'build_chunkset_from_anchor_manifest:MissingScaleColumn', ...
    'scaleTable must contain chunk_sec.');

scaleSec = double(scaleTable.chunk_sec(:));
nScale = numel(scaleSec);
nAnchors = height(anchorTable);

firstRow = anchorTable.feature_row_index(1);
firstDyad = local_load_dyad(sessionTable, firstRow, opts.repoRoot);
firstSeq = prepare_dyad_timeseries(firstDyad, 'stats', stats);
D = size(firstSeq.Xscaled, 2);
fps0 = firstSeq.fps;

Scale = repmat(struct( ...
    'chunkSec', [], ...
    'nFrames', [], ...
    'X', [], ...
    'Xraw', [], ...
    'valid', [], ...
    'meta', table(), ...
    'sessionIndex', [], ...
    'obsNames', {{}}, ...
    'channelMeta', table(), ...
    'featureNames', {{}}, ...
    'featureMeta', table()), nScale, 1);

anchorManifest = table();
for s = 1:nScale
    L = max(1, round(scaleSec(s) * fps0));
    meta = local_scale_meta(anchorTable, sessionTable, s, scaleTable, L, opts.anchorMode);
    Scale(s).chunkSec = scaleSec(s);
    Scale(s).nFrames = L;
    Scale(s).X = nan(nAnchors, L, D, 'single');
    Scale(s).Xraw = nan(nAnchors, L, D, 'single');
    Scale(s).valid = false(nAnchors, L);
    Scale(s).meta = meta;
    Scale(s).sessionIndex = anchorTable.session_index;
    Scale(s).obsNames = firstSeq.obsNames;
    Scale(s).channelMeta = firstSeq.channelMeta;
    Scale(s).featureNames = firstSeq.featureNames;
    Scale(s).featureMeta = firstSeq.featureMeta;
end

sessions = cell(height(sessionTable), 1);
for i = 1:height(sessionTable)
    rows = find(anchorTable.feature_row_index == i);
    if isempty(rows)
        sessions{i} = local_empty_session_struct(sessionTable, i);
        continue
    end

    dyad = local_load_dyad(sessionTable, i, opts.repoRoot);
    Seq = prepare_dyad_timeseries(dyad, 'stats', stats);
    sessions{i} = local_session_struct(sessionTable, i, Seq);

    for s = 1:nScale
        L = Scale(s).nFrames;
        [left, right] = local_chunk_geometry(L, opts.anchorMode);
        for rr = rows(:)'
            a = anchorTable.anchor_frame(rr);
            idx = (a - left):(a + right);
            if idx(1) < 1 || idx(end) > numel(Seq.time)
                continue
            end
            Xraw = Seq.X(idx, :);
            Xscaled = Seq.Xscaled(idx, :);
            v = logical(Seq.validMask(idx));
            Scale(s).X(rr, :, :) = single(Xscaled);
            Scale(s).Xraw(rr, :, :) = single(Xraw);
            Scale(s).valid(rr, :) = v(:)';
            validFrac = mean(v);
            finiteFrac = mean(isfinite(Xraw), 'all');

            Scale(s).meta.start_frame(rr) = idx(1);
            Scale(s).meta.stop_frame(rr) = idx(end);
            Scale(s).meta.start_time_s(rr) = Seq.time(idx(1));
            Scale(s).meta.stop_time_s(rr) = Seq.time(idx(end));
            Scale(s).meta.valid_frac(rr) = validFrac;
            Scale(s).meta.feature_finite_frac(rr) = finiteFrac;
            Scale(s).meta.chunk_is_valid(rr) = ...
                validFrac >= opts.minValidFrac && finiteFrac >= opts.minFeatureFiniteFrac;
        end
    end
end

chunkTable = table();
for s = 1:nScale
    Scale(s).meta = local_finalize_scale_meta(Scale(s).meta, opts);
    anchorManifest = [anchorManifest; Scale(s).meta]; %#ok<AGROW>
    chunkTable = [chunkTable; Scale(s).meta]; %#ok<AGROW>
end
chunkTable = sortrows(chunkTable, 'chunk_id');
anchorManifest = sortrows(anchorManifest, {'scale_index', 'anchor_id'});

ChunkSet = struct();
ChunkSet.stats = stats;
ChunkSet.obsNames = firstSeq.obsNames;
ChunkSet.channelMeta = firstSeq.channelMeta;
ChunkSet.featureNames = firstSeq.featureNames;
ChunkSet.featureMeta = firstSeq.featureMeta;
ChunkSet.chunkSec = scaleSec(:)';
ChunkSet.strideSec = median(diff(unique(anchorTable.anchor_time_s)), 'omitnan');
if ~isfinite(ChunkSet.strideSec)
    ChunkSet.strideSec = NaN;
end
ChunkSet.anchorMode = char(opts.anchorMode);
ChunkSet.scale = Scale;
ChunkSet.chunkTable = chunkTable;
ChunkSet.sessions = sessions;
ChunkSet.nSessions = height(sessionTable);
ChunkSet.nObs = D;
ChunkSet.scaleBankMeta = scaleTable;
ChunkSet.anchorTable = anchorTable;
ChunkSet.materializedAnchorCount = nAnchors;
ChunkSet.manifestChunkRows = height(anchorManifest);
end

function meta = local_scale_meta(anchorTable, sessionTable, s, scaleTable, L, anchorMode)
n = height(anchorTable);
meta = anchorTable;
meta.scale_index = repmat(s, n, 1);
meta.chunk_sec = repmat(double(scaleTable.chunk_sec(s)), n, 1);
meta.chunk_frames = repmat(L, n, 1);
meta.anchor_mode = repmat(string(anchorMode), n, 1);
if ismember('temporal_band', scaleTable.Properties.VariableNames)
    meta.temporal_band = repmat(string(scaleTable.temporal_band(s)), n, 1);
else
    meta.temporal_band = repmat(local_scale_band_label(scaleTable.chunk_sec(s)), n, 1);
end
if ismember('hierarchical_role', scaleTable.Properties.VariableNames)
    meta.hierarchical_role = repmat(string(scaleTable.hierarchical_role(s)), n, 1);
else
    meta.hierarchical_role = meta.temporal_band;
end
meta.chunk_id = ((1:n)' + (s - 1) * n);
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

function meta = local_finalize_scale_meta(meta, opts)
meta.inclusion_flag = meta.chunk_is_valid;
meta.inclusion_rule = repmat( ...
    "predefined_condition_blind_anchor_and_chunk_validity_thresholds", ...
    height(meta), 1);
meta.validity_threshold = repmat(opts.minValidFrac, height(meta), 1);
meta.feature_finite_threshold = repmat(opts.minFeatureFiniteFrac, height(meta), 1);
end

function dyad = local_load_dyad(sessionTable, rowIdx, repoRoot)
pathText = string(sessionTable.feature_output_file(rowIdx));
featurePath = resolve_repo_path(repoRoot, pathText);
S = load(featurePath, 'dyad', 'status');
assert(isfield(S, 'dyad'), 'build_chunkset_from_anchor_manifest:MissingDyad', ...
    'Missing dyad in %s.', featurePath);
if isfield(S, 'status')
    assert(string(S.status) == "success", ...
        'build_chunkset_from_anchor_manifest:BadFeatureStatus', ...
        'Feature MAT %s does not report success.', featurePath);
end
dyad = S.dyad;
end

function S = local_session_struct(sessionTable, rowIdx, Seq)
S = local_empty_session_struct(sessionTable, rowIdx);
S.fps = Seq.fps;
S.nFrames = numel(Seq.time);
S.validFrameFractionFromSeq = mean(Seq.validMask);
end

function S = local_empty_session_struct(sessionTable, rowIdx)
S = struct();
S.sessionIndex = rowIdx;
S.raw_index = local_table_double(sessionTable, 'raw_index', rowIdx, rowIdx);
S.session_id = local_table_string(sessionTable, 'session_id', rowIdx, "");
S.session_file = local_table_string(sessionTable, 'session_file', rowIdx, "");
S.arena = local_table_string(sessionTable, 'arena', rowIdx, "");
S.arena_label = local_table_string(sessionTable, 'arena_label', rowIdx, "");
S.cohort_id = local_table_string(sessionTable, 'cohort_id', rowIdx, "");
S.condition_id = local_table_string(sessionTable, 'condition_id', rowIdx, "");
S.qc_set = local_table_string(sessionTable, 'qc_set', rowIdx, "");
S.fps = NaN;
S.nFrames = NaN;
S.validFrameFractionFromSeq = NaN;
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
