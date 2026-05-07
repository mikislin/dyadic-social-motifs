function Sweep = run_motif_k_sweep(ChunkSet, EmbedModel, selectedScales, varargin)
%RUN_MOTIF_K_SWEEP Cluster/segment/evaluate a range of motif counts.
%
% Sweep = run_motif_k_sweep(ChunkSet, EmbedModel, selectedScales, ...)
%
% Runs cluster_multiscale_chunks, build_motif_segments, and
% compute_motif_transition_stats for each K. The output is intentionally
% segment-aware: ARI alone is not enough.
%
% Required inputs follow the MultiScaleChunks conventions:
%   ChunkSet, EmbedModel, selectedScales, Cluster, Segments.
%
% Important outputs:
%   Sweep.table      one row per K
%   Sweep.Clusters   optional Cluster objects if KeepObjects = true
%   Sweep.Segments   optional Segments objects if KeepObjects = true
%   Sweep.Transitions optional transition stats if KeepObjects = true

p = inputParser;
p.addParameter('KList', 4:10, @(x)isnumeric(x) && isvector(x) && all(isfinite(x)) && all(x >= 2));
p.addParameter('NumPCs', 12, @(x)isscalar(x) && x >= 2);
p.addParameter('ClusterMethod', 'gmm', @(x)ischar(x) || isstring(x));
p.addParameter('RegularizationValue', 1e-5, @(x)isscalar(x) && x >= 0);
p.addParameter('Replicates', 5, @(x)isscalar(x) && x >= 1);
p.addParameter('RandomSeed', 1, @isscalar);
p.addParameter('StabilityBootstrapN', 10, @(x)isscalar(x) && x >= 0);
p.addParameter('AnchorSetMode', 'reference', @(x)ischar(x) || isstring(x));
p.addParameter('AnchorAlignment', 'nearest', @(x)ischar(x) || isstring(x));
p.addParameter('RequireAllSelectedScales', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('FrameRate', [], @(x)isempty(x) || (isscalar(x) && x > 0));
p.addParameter('MinPosterior', 0.75, @(x)isscalar(x) && x >= 0 && x <= 1);
p.addParameter('SmoothWindowAnchors', 3, @(x)isscalar(x) && x >= 1);
p.addParameter('MinRunAnchors', 4, @(x)isscalar(x) && x >= 1);
p.addParameter('MergeGapAnchors', 2, @(x)isscalar(x) && x >= 0);
p.addParameter('MinSegmentDurationSec', 0.5, @(x)isscalar(x) && x >= 0);
p.addParameter('KeepObjects', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('Verbose', true, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});
P = p.Results;

KList = P.KList(:)';
if any(abs(KList - round(KList)) > eps(max(1, max(abs(KList)))))
    error('run_motif_k_sweep:KListMustBeInteger', ...
        'KList must contain integer motif counts, for example 4:10 or [4 5 6 7 8 9 10].');
end
KList = unique(round(KList), 'stable');

rows = repmat(local_empty_row(), numel(KList), 1);
Clusters = cell(numel(KList),1);
SegmentsCell = cell(numel(KList),1);
TransitionsCell = cell(numel(KList),1);

for ii = 1:numel(KList)
    K = KList(ii);
    if logical(P.Verbose)
        fprintf('\n=== K sweep %d/%d | K = %d ===\n', ii, numel(KList), K);
    end

    Cluster = cluster_multiscale_chunks(ChunkSet, EmbedModel, selectedScales, ...
        'ClusterMethod', P.ClusterMethod, ...
        'NumClusters', K, ...
        'NumPCs', P.NumPCs, ...
        'RegularizationValue', P.RegularizationValue, ...
        'Replicates', P.Replicates, ...
        'RandomSeed', P.RandomSeed + ii - 1, ...
        'StabilityBootstrapN', P.StabilityBootstrapN, ...
        'AnchorSetMode', P.AnchorSetMode, ...
        'AnchorAlignment', P.AnchorAlignment, ...
        'RequireAllSelectedScales', logical(P.RequireAllSelectedScales), ...
        'Verbose', logical(P.Verbose));

    Segments = build_motif_segments(Cluster, ...
        'MinPosterior', P.MinPosterior, ...
        'SmoothWindowAnchors', P.SmoothWindowAnchors, ...
        'MinRunAnchors', P.MinRunAnchors, ...
        'MergeGapAnchors', P.MergeGapAnchors, ...
        'MinSegmentDurationSec', P.MinSegmentDurationSec, ...
        'FrameRate', P.FrameRate, ...
        'Verbose', false);

    Transitions = compute_motif_transition_stats(Segments, ...
        'NumMotifs', K, ...
        'Verbose', false);

    summ = local_cluster_summary(Cluster);
    segSumm = local_segment_summary(Segments);

    rows(ii).K = K;
    rows(ii).n_anchors = summ.nAnchors;
    rows(ii).n_segments = segSumm.n_segments;
    rows(ii).median_segment_duration_s = segSumm.median_segment_duration_s;
    rows(ii).median_segment_anchors = segSumm.median_segment_anchors;
    rows(ii).short_duration_fraction = segSumm.short_duration_fraction;
    rows(ii).meanARI = summ.meanARI;
    rows(ii).median_max_posterior = summ.medianMaxPosterior;
    rows(ii).min_occupancy_fraction = summ.minOccupancyFraction;
    rows(ii).max_occupancy_fraction = summ.maxOccupancyFraction;
    rows(ii).occupancy_entropy_bits = local_entropy(summ.clusterFrac);
    rows(ii).transition_entropy_bits = Transitions.globalEntropy;
    rows(ii).n_nonzero_transitions = nnz(Transitions.counts);

    if isfield(Transitions, 'outgoingEntropy')
        rows(ii).mean_outgoing_entropy_bits = mean(Transitions.outgoingEntropy, 'omitnan');
        rows(ii).max_outgoing_entropy_bits = max(Transitions.outgoingEntropy, [], 'omitnan');
    end

    if logical(P.KeepObjects)
        Clusters{ii} = Cluster;
        SegmentsCell{ii} = Segments;
        TransitionsCell{ii} = Transitions;
    end
end

Sweep = struct();
Sweep.table = struct2table(rows);
Sweep.KList = KList;
Sweep.params = P;
Sweep.Clusters = Clusters;
Sweep.Segments = SegmentsCell;
Sweep.Transitions = TransitionsCell;

if logical(P.Verbose)
    fprintf('\n=== K sweep summary ===\n');
    disp(Sweep.table);
end
end

function row = local_empty_row()
row = struct();
row.K = NaN;
row.n_anchors = NaN;
row.n_segments = NaN;
row.median_segment_duration_s = NaN;
row.median_segment_anchors = NaN;
row.short_duration_fraction = NaN;
row.meanARI = NaN;
row.median_max_posterior = NaN;
row.min_occupancy_fraction = NaN;
row.max_occupancy_fraction = NaN;
row.occupancy_entropy_bits = NaN;
row.transition_entropy_bits = NaN;
row.n_nonzero_transitions = NaN;
row.mean_outgoing_entropy_bits = NaN;
row.max_outgoing_entropy_bits = NaN;
end

function S = local_cluster_summary(Cluster)
S = struct();

S.nAnchors = local_getfield_flexible(Cluster, { ...
    {'summary','nAnchors'}, {'summary','n_anchors'}, {'nAnchors'}, {'NumAnchors'}}, NaN);
if isnan(S.nAnchors)
    if isfield(Cluster, 'labels')
        S.nAnchors = numel(Cluster.labels);
    elseif isfield(Cluster, 'anchorTable') && istable(Cluster.anchorTable)
        S.nAnchors = height(Cluster.anchorTable);
    elseif isfield(Cluster, 'Data') && isfield(Cluster.Data, 'anchorTable') && istable(Cluster.Data.anchorTable)
        S.nAnchors = height(Cluster.Data.anchorTable);
    end
end

S.meanARI = local_getfield_flexible(Cluster, { ...
    {'summary','meanARI'}, {'stability','meanARI'}, {'meanARI'}}, NaN);

S.medianMaxPosterior = local_getfield_flexible(Cluster, { ...
    {'summary','medianMaxPosterior'}, {'medianMaxPosterior'}}, NaN);
if isnan(S.medianMaxPosterior) && isfield(Cluster, 'maxPosterior')
    S.medianMaxPosterior = median(Cluster.maxPosterior, 'omitnan');
end

S.clusterFrac = local_getfield_flexible(Cluster, { ...
    {'summary','clusterFrac'}, {'clusterFrac'}, {'occupancyFrac'}}, []);
if isempty(S.clusterFrac)
    counts = local_getfield_flexible(Cluster, { ...
        {'summary','clusterCounts'}, {'clusterCounts'}, {'occupancyCounts'}}, []);
    if isempty(counts) && isfield(Cluster, 'labels')
        labels = Cluster.labels(:);
        labels = labels(isfinite(labels) & labels > 0);
        K = local_getfield_flexible(Cluster, {{'NumClusters'}, {'K'}}, max(labels));
        counts = accumarray(labels, 1, [K 1], @sum, 0);
    end
    if ~isempty(counts) && sum(counts) > 0
        S.clusterFrac = counts(:) ./ sum(counts);
    end
end

if isempty(S.clusterFrac)
    S.clusterFrac = NaN;
    S.minOccupancyFraction = NaN;
    S.maxOccupancyFraction = NaN;
else
    S.clusterFrac = S.clusterFrac(:);
    S.minOccupancyFraction = min(S.clusterFrac, [], 'omitnan');
    S.maxOccupancyFraction = max(S.clusterFrac, [], 'omitnan');
end
end

function S = local_segment_summary(Segments)
S = struct();
S.n_segments = local_getfield_flexible(Segments, {{'summary','n_segments'}, {'summary','nSegments'}, {'n_segments'}}, NaN);
S.median_segment_duration_s = local_getfield_flexible(Segments, {{'summary','median_segment_duration_s'}, {'summary','medianSegmentDurationSec'}}, NaN);
S.median_segment_anchors = local_getfield_flexible(Segments, {{'summary','median_segment_anchors'}, {'summary','medianSegmentAnchors'}}, NaN);
S.short_duration_fraction = local_getfield_flexible(Segments, {{'summary','short_duration_fraction'}, {'summary','shortDurationFraction'}}, NaN);

if (isnan(S.n_segments) || isnan(S.median_segment_duration_s)) && isfield(Segments, 'table') && istable(Segments.table)
    T = Segments.table;
    S.n_segments = height(T);
    if ismember('duration_s', T.Properties.VariableNames)
        S.median_segment_duration_s = median(T.duration_s, 'omitnan');
    end
    if ismember('n_anchors', T.Properties.VariableNames)
        S.median_segment_anchors = median(T.n_anchors, 'omitnan');
    end
end
end

function v = local_getfield_flexible(S, paths, fallback)
v = fallback;
for p = 1:numel(paths)
    try
        tmp = S;
        names = paths{p};
        for j = 1:numel(names)
            tmp = tmp.(names{j});
        end
        if ~isempty(tmp)
            v = tmp;
            return
        end
    catch
    end
end
end

function H = local_entropy(p)
p = p(:);
p = p(isfinite(p) & p > 0);
if isempty(p)
    H = NaN;
else
    p = p ./ sum(p);
    H = -sum(p .* log2(p));
end
end
