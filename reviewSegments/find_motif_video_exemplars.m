function Exemplars = find_motif_video_exemplars(Cluster, Segments, varargin)
%FIND_MOTIF_VIDEO_EXEMPLARS Rank segment/bout exemplars near motif medoids.
%
% Exemplars = find_motif_video_exemplars(Cluster, Segments, ...)
%
% This function does not require video files. It returns a table of the top
% N medoid-neighbor bouts per motif. If VideoTable or VideoPath is supplied,
% it also adds video path and suggested clip window columns.
%
% Name-value options
%   TopN                  10
%   WindowSec             1.0       padding around segment for video clip
%   MinPosterior          0         minimum segment median posterior
%   VideoTable            table()   session_index, video_path columns
%   VideoPath             ''        same video for all sessions
%   FrameRate             []        used if frame/time conversion needed
%   Verbose               true

p = inputParser;
p.addParameter('TopN', 10, @(x)isscalar(x) && x >= 1);
p.addParameter('WindowSec', 1.0, @(x)isscalar(x) && x >= 0);
p.addParameter('MinPosterior', 0, @(x)isscalar(x) && x >= 0 && x <= 1);
p.addParameter('VideoTable', table(), @(x)istable(x));
p.addParameter('VideoPath', '', @(x)ischar(x) || isstring(x));
p.addParameter('FrameRate', [], @(x)isempty(x) || (isscalar(x) && x > 0));
p.addParameter('Verbose', true, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});
P = p.Results;

if ~isfield(Cluster, 'Xpca') || ~isfield(Cluster, 'Data') || ~isfield(Cluster.Data, 'anchorTable')
    error('find_motif_video_exemplars:BadCluster', 'Cluster must contain Xpca and Data.anchorTable.');
end
T = Segments.table;
labelName = local_label_name(T);
K = Cluster.NumClusters;
A = Cluster.Data.anchorTable;
Xpca = Cluster.Xpca;
medoidIdx = Cluster.medoidIdx(:);

rows = struct([]);
r = 0;
for k = 1:K
    segIdx = find(double(T.(labelName)) == k);
    candidates = [];
    for ii = 1:numel(segIdx)
        ti = segIdx(ii);
        if ismember('median_posterior', T.Properties.VariableNames) && T.median_posterior(ti) < P.MinPosterior
            continue
        end
        anchorIdx = local_segment_anchor_indices(T(ti,:), A);
        if isempty(anchorIdx)
            continue
        end
        xseg = mean(Xpca(anchorIdx,:), 1, 'omitnan');
        if isfinite(medoidIdx(k)) && medoidIdx(k) >= 1 && medoidIdx(k) <= size(Xpca,1)
            d = sqrt(sum((xseg - Xpca(medoidIdx(k),:)).^2));
        else
            d = NaN;
        end
        candidates(end+1).segment_row = ti; %#ok<AGROW>
        candidates(end).distance_to_medoid = d;
        candidates(end).n_anchor_overlap = numel(anchorIdx);
    end
    if isempty(candidates)
        continue
    end
    C = struct2table(candidates);
    C = sortrows(C, 'distance_to_medoid', 'ascend');
    nTake = min(P.TopN, height(C));
    for jj = 1:nTake
        ti = C.segment_row(jj);
        r = r + 1;
        rows(r).motif_id = k; %#ok<AGROW>
        rows(r).rank = jj;
        rows(r).segment_row = ti;
        rows(r).session_index = T.session_index(ti);
        rows(r).distance_to_medoid = C.distance_to_medoid(jj);
        rows(r).n_anchor_overlap = C.n_anchor_overlap(jj);
        rows(r).start_frame = local_get_value(T, ti, 'start_frame', NaN);
        rows(r).stop_frame = local_get_value(T, ti, 'stop_frame', NaN);
        rows(r).start_time_s = local_segment_time(T, ti, 'start_time_s', rows(r).start_frame, P.FrameRate);
        rows(r).stop_time_s = local_segment_time(T, ti, 'stop_time_s', rows(r).stop_frame, P.FrameRate);
        rows(r).clip_start_time_s = max(0, rows(r).start_time_s - P.WindowSec);
        rows(r).clip_stop_time_s = rows(r).stop_time_s + P.WindowSec;
        rows(r).median_posterior = local_get_value(T, ti, 'median_posterior', NaN);
        rows(r).duration_s = local_get_value(T, ti, 'duration_s', rows(r).stop_time_s - rows(r).start_time_s);
    end
end

if isempty(rows)
    exemplarTable = table();
else
    exemplarTable = struct2table(rows);
end
exemplarTable = local_add_video_paths(exemplarTable, P.VideoTable, P.VideoPath);

Exemplars = struct();
Exemplars.table = exemplarTable;
Exemplars.params = P;

if logical(P.Verbose)
    fprintf('find_motif_video_exemplars | motifs = %d | exemplars = %d | topN = %d\n', ...
        K, height(exemplarTable), P.TopN);
end
end

function name = local_label_name(T)
vars = string(T.Properties.VariableNames);
if any(vars == "motif_id")
    name = 'motif_id';
elseif any(vars == "cluster")
    name = 'cluster';
elseif any(vars == "cluster_id")
    name = 'cluster_id';
else
    error('find_motif_video_exemplars:NoLabel', 'Segments.table needs motif_id, cluster, or cluster_id.');
end
end

function anchorIdx = local_segment_anchor_indices(segRow, A)
s = segRow.session_index(1);
mask = A.session_index == s;
if ismember('start_frame', segRow.Properties.VariableNames) && ismember('stop_frame', segRow.Properties.VariableNames)
    mask = mask & A.anchor_frame >= segRow.start_frame(1) & A.anchor_frame <= segRow.stop_frame(1);
elseif ismember('start_time_s', segRow.Properties.VariableNames) && ismember('anchor_time_s', A.Properties.VariableNames)
    mask = mask & A.anchor_time_s >= segRow.start_time_s(1) & A.anchor_time_s <= segRow.stop_time_s(1);
end
anchorIdx = find(mask);
end

function v = local_get_value(T, i, name, fallback)
if ismember(name, T.Properties.VariableNames)
    v = T.(name)(i);
else
    v = fallback;
end
end

function t = local_segment_time(T, i, name, frameVal, fps)
if ismember(name, T.Properties.VariableNames)
    t = T.(name)(i);
elseif ~isempty(fps) && isfinite(frameVal)
    t = double(frameVal) ./ double(fps);
else
    t = NaN;
end
end

function E = local_add_video_paths(E, videoTable, videoPath)
if isempty(E)
    return
end
paths = strings(height(E),1);
if strlength(string(videoPath)) > 0
    paths(:) = string(videoPath);
elseif ~isempty(videoTable) && istable(videoTable)
    if ~all(ismember({'session_index','video_path'}, videoTable.Properties.VariableNames))
        warning('find_motif_video_exemplars:BadVideoTable', 'VideoTable must contain session_index and video_path.');
    else
        for i = 1:height(E)
            idx = find(videoTable.session_index == E.session_index(i), 1, 'first');
            if ~isempty(idx)
                paths(i) = string(videoTable.video_path(idx));
            end
        end
    end
end
E.video_path = paths;
E.clip_label = strings(height(E),1);
for i = 1:height(E)
    E.clip_label(i) = sprintf('motif%02d_rank%02d_session%d_%0.2f-%0.2fs', ...
        E.motif_id(i), E.rank(i), E.session_index(i), E.clip_start_time_s(i), E.clip_stop_time_s(i));
end
end
