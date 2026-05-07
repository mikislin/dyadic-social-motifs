function Exemplars = find_motif_video_exemplars(Cluster, Segments, varargin)
%FIND_MOTIF_VIDEO_EXEMPLARS Rank motif bouts and optionally render skeleton exemplars.
%
% Exemplars = find_motif_video_exemplars(Cluster, Segments, ...)
%
% This function is designed for the MultiScaleChunks pipeline. It works even
% when raw video is not available. If preprocessed pose points are supplied,
% it can render skeleton-only exemplar clips for the top-ranked bouts.
%
% Required inputs
%   Cluster  - output of cluster_multiscale_chunks
%   Segments - output of build_motif_segments
%
% Name-value options
%   'TopN'              number of exemplars per motif (default: 10)
%   'WindowSec'         clip half-window around the exemplar anchor, unless
%                       the bout duration is longer (default: 1.0)
%   'FrameRate'         video/pose frame rate in Hz (default: inferred or 80)
%   'PoseData'          optional pose source for skeleton rendering. Accepted:
%                         DyadSet-like struct, dyadCell, numeric array,
%                         cell array, or dbase-like struct array
%   'MakeSkeletonClips' true/false, render skeleton clips if PoseData exists
%                       (default: false)
%   'OutputDir'         output folder for skeleton clips (default: '')
%   'SkeletonEdges'     E-by-2 edge list for bodypart skeleton (default: mouse)
%   'BodypartNames'     bodypart labels for inferred/default edges (default: {})
%   'MarkerSize'        marker size for skeleton video (default: 40)
%   'LineWidth'         skeleton line width (default: 2)
%   'PadSec'            extra padding around segment for clip bounds (default: 0)
%   'MaxFramesPerClip'  cap rendered frames per clip (default: 240)
%   'VideoFormat'       'mp4' or 'avi' (default: 'mp4')
%   'Verbose'           print summary (default: true)
%
% Output
%   Exemplars.table          ranked exemplar bouts
%   Exemplars.byMotif        cell array, one table per motif
%   Exemplars.settings       settings used
%   Exemplars.skeletonFiles  table of rendered files, if requested
%   Exemplars.notes          notes/warnings
%
% The ranking prioritizes segments close to the motif medoid in Cluster.Xpca,
% high posterior confidence, and sufficient bout duration.

p = inputParser;
p.addRequired('Cluster', @isstruct);
p.addRequired('Segments', @isstruct);
p.addParameter('TopN', 10, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('WindowSec', 1.0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('FrameRate', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
p.addParameter('PoseData', [], @(x) true);
p.addParameter('MakeSkeletonClips', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('OutputDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('SkeletonEdges', [], @(x) isempty(x) || (isnumeric(x) && size(x,2) == 2));
p.addParameter('BodypartNames', {}, @(x) iscellstr(x) || isstring(x) || isempty(x));
p.addParameter('MarkerSize', 40, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('LineWidth', 2, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('PadSec', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('MaxFramesPerClip', 240, @(x) isnumeric(x) && isscalar(x) && x >= 10);
p.addParameter('VideoFormat', 'mp4', @(x) any(strcmpi(string(x), ["mp4","avi"])));
p.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
p.parse(Cluster, Segments, varargin{:});
P = p.Results;
P.TopN = round(P.TopN);
P.MakeSkeletonClips = logical(P.MakeSkeletonClips);
P.OutputDir = char(P.OutputDir);
P.VideoFormat = lower(char(P.VideoFormat));
P.BodypartNames = cellstr(P.BodypartNames);
P.Verbose = logical(P.Verbose);

notes = strings(0,1);

fps = local_get_fps(P.FrameRate, Cluster, Segments);
P.FrameRate = fps;

segTable = local_get_segment_table(Segments);
anchorTable = local_get_anchor_table(Cluster, Segments);
labels = local_get_labels(Cluster, Segments, anchorTable);
post = local_get_max_posterior(Cluster, Segments, anchorTable);
Xpca = local_get_xpca(Cluster);
medoidAnchorIdx = local_get_medoid_anchor_indices(Cluster, labels, Xpca);

if isempty(segTable)
    error('find_motif_video_exemplars:MissingSegments', ...
        'Segments.table is required and was not found.');
end

requiredVars = {'motif_id','session_index'};
for i = 1:numel(requiredVars)
    if ~ismember(requiredVars{i}, segTable.Properties.VariableNames)
        error('find_motif_video_exemplars:MissingSegmentField', ...
            'Segments.table must contain %s.', requiredVars{i});
    end
end

segTable = local_add_segment_time_fields(segTable, fps);
numMotifs = local_get_num_motifs(Cluster, Segments, segTable);
rows = table();
byMotif = cell(numMotifs,1);

for k = 1:numMotifs
    motifRows = find(segTable.motif_id == k);
    if isempty(motifRows)
        byMotif{k} = table();
        continue
    end

    motifScores = nan(numel(motifRows),1);
    motifMedoidDist = nan(numel(motifRows),1);
    motifPosterior = nan(numel(motifRows),1);
    motifCenterAnchor = nan(numel(motifRows),1);
    motifCenterFrame = nan(numel(motifRows),1);
    motifCenterTime = nan(numel(motifRows),1);

    for ii = 1:numel(motifRows)
        r = motifRows(ii);
        [anchorIdx, centerAnchorIdx] = local_segment_anchor_indices(segTable(r,:), anchorTable);
        motifCenterAnchor(ii) = centerAnchorIdx;

        if ~isnan(centerAnchorIdx) && centerAnchorIdx >= 1 && centerAnchorIdx <= height(anchorTable)
            if ismember('anchor_frame', anchorTable.Properties.VariableNames)
                motifCenterFrame(ii) = anchorTable.anchor_frame(centerAnchorIdx);
            end
            if ismember('anchor_time_s', anchorTable.Properties.VariableNames)
                motifCenterTime(ii) = anchorTable.anchor_time_s(centerAnchorIdx);
            elseif isfinite(motifCenterFrame(ii))
                motifCenterTime(ii) = motifCenterFrame(ii) ./ fps;
            end
        end

        if ~isempty(anchorIdx) && ~isempty(post)
            idx = anchorIdx(anchorIdx >= 1 & anchorIdx <= numel(post));
            if ~isempty(idx)
                motifPosterior(ii) = median(post(idx), 'omitnan');
            end
        end

        if ~isempty(Xpca) && ~isnan(centerAnchorIdx) && ~isnan(medoidAnchorIdx(k))
            if centerAnchorIdx <= size(Xpca,1) && medoidAnchorIdx(k) <= size(Xpca,1)
                d = Xpca(centerAnchorIdx,:) - Xpca(medoidAnchorIdx(k),:);
                motifMedoidDist(ii) = sqrt(sum(d.^2, 'omitnan'));
            end
        end

        dur = segTable.duration_s(r);
        durScore = min(dur ./ max(P.WindowSec, 0.1), 2);
        distScore = 1 ./ (1 + motifMedoidDist(ii));
        if isnan(distScore)
            distScore = 0.5;
        end
        postScore = motifPosterior(ii);
        if isnan(postScore)
            postScore = 0.5;
        end
        motifScores(ii) = 0.55 * distScore + 0.30 * postScore + 0.15 * durScore;
    end

    [~, order] = sortrows([-motifScores, motifMedoidDist, -motifPosterior]);
    order = order(1:min(P.TopN, numel(order)));
    keepRows = motifRows(order);

    T = segTable(keepRows,:);
    T.rank = (1:height(T))';
    T.center_anchor_index = motifCenterAnchor(order);
    T.center_frame = motifCenterFrame(order);
    T.center_time_s = motifCenterTime(order);
    T.medoid_distance = motifMedoidDist(order);
    T.segment_median_posterior = motifPosterior(order);
    T.exemplar_score = motifScores(order);
    T.clip_start_frame = max(1, floor(T.center_frame - (P.WindowSec + P.PadSec) * fps));
    T.clip_end_frame = ceil(T.center_frame + (P.WindowSec + P.PadSec) * fps);
    T.clip_start_time_s = T.clip_start_frame ./ fps;
    T.clip_end_time_s = T.clip_end_frame ./ fps;
    T.clip_file = strings(height(T),1);
    byMotif{k} = T;
    rows = [rows; T]; %#ok<AGROW>
end

Exemplars = struct();
Exemplars.table = rows;
Exemplars.byMotif = byMotif;
Exemplars.settings = P;
Exemplars.notes = notes;
Exemplars.skeletonFiles = table();

if P.MakeSkeletonClips
    if isempty(P.OutputDir)
        error('find_motif_video_exemplars:OutputDirRequired', ...
            'OutputDir is required when MakeSkeletonClips is true.');
    end
    if isempty(P.PoseData)
        notes(end+1,1) = "MakeSkeletonClips requested but PoseData was empty; no clips rendered.";
    elseif isempty(rows)
        notes(end+1,1) = "No exemplar rows found; no clips rendered.";
    else
        if ~exist(P.OutputDir, 'dir')
            mkdir(P.OutputDir);
        end
        [fileTable, renderNotes] = render_motif_skeleton_exemplars(P.PoseData, rows, ...
            'OutputDir', P.OutputDir, ...
            'FrameRate', fps, ...
            'SkeletonEdges', P.SkeletonEdges, ...
            'BodypartNames', P.BodypartNames, ...
            'MarkerSize', P.MarkerSize, ...
            'LineWidth', P.LineWidth, ...
            'MaxFramesPerClip', P.MaxFramesPerClip, ...
            'VideoFormat', P.VideoFormat, ...
            'Verbose', P.Verbose);
        Exemplars.skeletonFiles = fileTable;
        Exemplars.notes = [notes; renderNotes(:)];
        if ~isempty(fileTable) && ismember('clip_file', fileTable.Properties.VariableNames)
            [tf, loc] = ismember(Exemplars.table.exemplar_id, fileTable.exemplar_id);
            if any(tf)
                Exemplars.table.clip_file(tf) = string(fileTable.clip_file(loc(tf)));
            end
        end
    end
end

if ~ismember('exemplar_id', Exemplars.table.Properties.VariableNames) && ~isempty(Exemplars.table)
    Exemplars.table.exemplar_id = compose('motif%02d_rank%02d', Exemplars.table.motif_id, Exemplars.table.rank);
end

if P.Verbose
    fprintf('find_motif_video_exemplars | motifs = %d | exemplars = %d | topN = %d', ...
        numMotifs, height(Exemplars.table), P.TopN);
    if P.MakeSkeletonClips
        fprintf(' | skeleton clips = %d', height(Exemplars.skeletonFiles));
    end
    fprintf('\n');
end
end

function fps = local_get_fps(frameRate, Cluster, Segments)
fps = frameRate;
if ~isempty(fps) && isfinite(fps) && fps > 0
    fps = double(fps);
    return
end
fps = [];
objs = {Cluster, Segments};
for i = 1:numel(objs)
    obj = objs{i};
    if isstruct(obj)
        if isfield(obj, 'FrameRate') && isnumeric(obj.FrameRate) && isscalar(obj.FrameRate)
            fps = obj.FrameRate; return
        end
        if isfield(obj, 'settings') && isstruct(obj.settings) && isfield(obj.settings, 'FrameRate')
            fps = obj.settings.FrameRate; return
        end
    end
end
fps = 80;
end

function T = local_get_segment_table(Segments)
T = table();
if isfield(Segments, 'table') && istable(Segments.table)
    T = Segments.table;
elseif istable(Segments)
    T = Segments;
end
end

function A = local_get_anchor_table(Cluster, Segments)
A = table();
if isfield(Cluster, 'Data') && isstruct(Cluster.Data) && isfield(Cluster.Data, 'anchorTable') && istable(Cluster.Data.anchorTable)
    A = Cluster.Data.anchorTable;
elseif isfield(Cluster, 'anchorTable') && istable(Cluster.anchorTable)
    A = Cluster.anchorTable;
elseif isfield(Segments, 'anchorLabelTable') && istable(Segments.anchorLabelTable)
    A = Segments.anchorLabelTable;
end
end

function labels = local_get_labels(Cluster, Segments, anchorTable)
labels = [];
fields = {'labels','label','cluster_id','motif_id'};
for f = 1:numel(fields)
    if isfield(Cluster, fields{f})
        labels = Cluster.(fields{f})(:); return
    end
end
if istable(anchorTable)
    for f = 1:numel(fields)
        if ismember(fields{f}, anchorTable.Properties.VariableNames)
            labels = anchorTable.(fields{f})(:); return
        end
    end
end
if isfield(Segments, 'anchorLabelTable') && istable(Segments.anchorLabelTable)
    T = Segments.anchorLabelTable;
    for f = 1:numel(fields)
        if ismember(fields{f}, T.Properties.VariableNames)
            labels = T.(fields{f})(:); return
        end
    end
end
end

function post = local_get_max_posterior(Cluster, Segments, anchorTable)
post = [];
if isfield(Cluster, 'maxPosterior')
    post = Cluster.maxPosterior(:); return
end
if isfield(Cluster, 'posterior')
    P = Cluster.posterior;
    if ismatrix(P)
        post = max(P, [], 2); return
    end
end
if isfield(Cluster, 'posteriors')
    P = Cluster.posteriors;
    if ismatrix(P)
        post = max(P, [], 2); return
    end
end
if istable(anchorTable)
    candidates = {'max_posterior','posterior','segment_median_posterior'};
    for i = 1:numel(candidates)
        if ismember(candidates{i}, anchorTable.Properties.VariableNames)
            post = anchorTable.(candidates{i})(:); return
        end
    end
end
if isfield(Segments, 'anchorLabelTable') && istable(Segments.anchorLabelTable)
    T = Segments.anchorLabelTable;
    candidates = {'max_posterior','posterior'};
    for i = 1:numel(candidates)
        if ismember(candidates{i}, T.Properties.VariableNames)
            post = T.(candidates{i})(:); return
        end
    end
end
end

function X = local_get_xpca(Cluster)
X = [];
if isfield(Cluster, 'Xpca') && isnumeric(Cluster.Xpca)
    X = Cluster.Xpca;
elseif isfield(Cluster, 'X_pca') && isnumeric(Cluster.X_pca)
    X = Cluster.X_pca;
elseif isfield(Cluster, 'Data') && isstruct(Cluster.Data) && isfield(Cluster.Data, 'Xpca')
    X = Cluster.Data.Xpca;
end
end

function medoidIdx = local_get_medoid_anchor_indices(Cluster, labels, Xpca)
numMotifs = max(labels, [], 'omitnan');
medoidIdx = nan(numMotifs,1);
if isfield(Cluster, 'medoidIdx') && numel(Cluster.medoidIdx) >= numMotifs
    medoidIdx = double(Cluster.medoidIdx(:));
    medoidIdx = medoidIdx(1:numMotifs);
    return
end
if isempty(labels) || isempty(Xpca)
    return
end
for k = 1:numMotifs
    idx = find(labels == k);
    if isempty(idx)
        continue
    end
    X = Xpca(idx,:);
    c = median(X, 1, 'omitnan');
    d = sqrt(sum((X - c).^2, 2, 'omitnan'));
    [~, jj] = min(d);
    medoidIdx(k) = idx(jj);
end
end

function numMotifs = local_get_num_motifs(Cluster, Segments, segTable)
if isfield(Cluster, 'NumClusters') && ~isempty(Cluster.NumClusters)
    numMotifs = double(Cluster.NumClusters);
elseif isfield(Cluster, 'K') && ~isempty(Cluster.K)
    numMotifs = double(Cluster.K);
elseif isfield(Segments, 'motifSummary') && istable(Segments.motifSummary) && ismember('motif_id', Segments.motifSummary.Properties.VariableNames)
    numMotifs = max(Segments.motifSummary.motif_id);
else
    numMotifs = max(segTable.motif_id);
end
end

function T = local_add_segment_time_fields(T, fps)
if ~ismember('duration_s', T.Properties.VariableNames)
    if ismember('start_time_s', T.Properties.VariableNames) && ismember('end_time_s', T.Properties.VariableNames)
        T.duration_s = T.end_time_s - T.start_time_s;
    elseif ismember('start_frame', T.Properties.VariableNames) && ismember('end_frame', T.Properties.VariableNames)
        T.duration_s = (T.end_frame - T.start_frame + 1) ./ fps;
    elseif ismember('n_anchors', T.Properties.VariableNames)
        T.duration_s = T.n_anchors .* 0.1;
    else
        T.duration_s = nan(height(T),1);
    end
end
if ~ismember('start_frame', T.Properties.VariableNames)
    if ismember('start_anchor_frame', T.Properties.VariableNames)
        T.start_frame = T.start_anchor_frame;
    elseif ismember('start_time_s', T.Properties.VariableNames)
        T.start_frame = round(T.start_time_s .* fps);
    else
        T.start_frame = nan(height(T),1);
    end
end
if ~ismember('end_frame', T.Properties.VariableNames)
    if ismember('end_anchor_frame', T.Properties.VariableNames)
        T.end_frame = T.end_anchor_frame;
    elseif ismember('end_time_s', T.Properties.VariableNames)
        T.end_frame = round(T.end_time_s .* fps);
    else
        T.end_frame = nan(height(T),1);
    end
end
if ~ismember('start_time_s', T.Properties.VariableNames)
    T.start_time_s = T.start_frame ./ fps;
end
if ~ismember('end_time_s', T.Properties.VariableNames)
    T.end_time_s = T.end_frame ./ fps;
end
if ~ismember('segment_id', T.Properties.VariableNames)
    T.segment_id = (1:height(T))';
end
end

function [idx, centerIdx] = local_segment_anchor_indices(segRow, anchorTable)
idx = [];
centerIdx = nan;
if isempty(anchorTable) || height(anchorTable) == 0
    return
end
sess = segRow.session_index;
mask = anchorTable.session_index == sess;
if ismember('start_anchor_index', segRow.Properties.VariableNames) && ismember('end_anchor_index', segRow.Properties.VariableNames)
    a = segRow.start_anchor_index;
    b = segRow.end_anchor_index;
    if isfinite(a) && isfinite(b)
        idx = (max(1,a):min(height(anchorTable),b))';
    end
elseif ismember('anchor_index_start', segRow.Properties.VariableNames) && ismember('anchor_index_end', segRow.Properties.VariableNames)
    a = segRow.anchor_index_start;
    b = segRow.anchor_index_end;
    if isfinite(a) && isfinite(b)
        idx = (max(1,a):min(height(anchorTable),b))';
    end
elseif ismember('start_frame', segRow.Properties.VariableNames) && ismember('end_frame', segRow.Properties.VariableNames) && ismember('anchor_frame', anchorTable.Properties.VariableNames)
    mask = mask & anchorTable.anchor_frame >= segRow.start_frame & anchorTable.anchor_frame <= segRow.end_frame;
    idx = find(mask);
elseif ismember('start_time_s', segRow.Properties.VariableNames) && ismember('end_time_s', segRow.Properties.VariableNames) && ismember('anchor_time_s', anchorTable.Properties.VariableNames)
    mask = mask & anchorTable.anchor_time_s >= segRow.start_time_s & anchorTable.anchor_time_s <= segRow.end_time_s;
    idx = find(mask);
end
if isempty(idx)
    return
end
if ismember('start_frame', segRow.Properties.VariableNames) && ismember('end_frame', segRow.Properties.VariableNames) && ismember('anchor_frame', anchorTable.Properties.VariableNames)
    c = 0.5 * (segRow.start_frame + segRow.end_frame);
    [~, jj] = min(abs(anchorTable.anchor_frame(idx) - c));
    centerIdx = idx(jj);
elseif ismember('start_time_s', segRow.Properties.VariableNames) && ismember('end_time_s', segRow.Properties.VariableNames) && ismember('anchor_time_s', anchorTable.Properties.VariableNames)
    c = 0.5 * (segRow.start_time_s + segRow.end_time_s);
    [~, jj] = min(abs(anchorTable.anchor_time_s(idx) - c));
    centerIdx = idx(jj);
else
    centerIdx = idx(round((numel(idx)+1)/2));
end
end
