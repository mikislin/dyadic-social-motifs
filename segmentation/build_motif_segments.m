function Segments = build_motif_segments(Cluster, varargin)
%BUILD_MOTIF_SEGMENTS Convert anchor-level motif labels into bout-level segments.
%
% Segments = build_motif_segments(Cluster, Name, Value)
%
% Required input
%   Cluster : output of cluster_multiscale_chunks. The function uses only Cluster
%             and expects Cluster.Data.anchorTable, Cluster.labels, and preferably
%             Cluster.posteriors / Cluster.maxPosterior.
%
% Name-value options
%   MinPosterior              : anchors below this confidence are set to 0 before smoothing (default 0)
%   SmoothWindowAnchors       : odd moving-mode window in anchor units (default 3)
%   MinRunAnchors             : short same-label runs shorter than this are absorbed (default 3)
%   MergeGapAnchors           : merge same-label bouts separated by <= this many anchors (default 2)
%   MinSegmentDurationSec     : final bouts shorter than this are absorbed when possible (default 0.5)
%   FrameRate                 : fallback FPS if anchor_time_s is unavailable (default [])
%   FrameCountBySession       : optional vector/map/table of frame counts for frameLabelCell (default [])
%   Verbose                   : print summary (default true)
%
% Outputs
%   Segments.table            : bout table
%   Segments.frameLabelCell   : per-session frame labels when frame counts can be inferred
%   Segments.framePosteriorCell : per-session frame posterior confidence when possible
%   Segments.anchorLabelTable : anchor table with raw/smoothed labels and confidence
%   Segments.summary          : global summary table
%   Segments.motifSummary     : motif-level summary table
%   Segments.sessionSummary   : session-level summary table

p = inputParser;
p.addRequired('Cluster', @(x)isstruct(x));
p.addParameter('MinPosterior', 0, @(x)isscalar(x) && x >= 0 && x <= 1);
p.addParameter('SmoothWindowAnchors', 3, @(x)isscalar(x) && x >= 1);
p.addParameter('MinRunAnchors', 3, @(x)isscalar(x) && x >= 1);
p.addParameter('MergeGapAnchors', 2, @(x)isscalar(x) && x >= 0);
p.addParameter('MinSegmentDurationSec', 0.5, @(x)isscalar(x) && x >= 0);
p.addParameter('FrameRate', [], @(x)isempty(x) || (isscalar(x) && x > 0));
p.addParameter('FrameCountBySession', [], @(x)isempty(x) || isnumeric(x) || istable(x) || isa(x,'containers.Map'));
p.addParameter('Verbose', true, @(x)islogical(x) || isnumeric(x));
p.parse(Cluster, varargin{:});
P = p.Results;
P.SmoothWindowAnchors = max(1, round(P.SmoothWindowAnchors));
if mod(P.SmoothWindowAnchors,2) == 0
    P.SmoothWindowAnchors = P.SmoothWindowAnchors + 1;
end
P.MinRunAnchors = max(1, round(P.MinRunAnchors));
P.MergeGapAnchors = max(0, round(P.MergeGapAnchors));

[A, rawLabel, post, maxPost] = local_unpack_cluster(Cluster);
N = height(A);
assert(numel(rawLabel) == N, 'Cluster labels do not match anchor table height.');

[~, ord] = sortrows([A.session_index(:), A.anchor_frame(:)], [1 2]);
invOrd = nan(N,1); invOrd(ord) = (1:N)';
Aord = A(ord,:);
rawOrd = rawLabel(ord);
maxPostOrd = maxPost(ord);
postOrd = post(ord,:);

labelConf = rawOrd;
labelConf(maxPostOrd < P.MinPosterior) = 0;

smoothOrd = zeros(size(labelConf));
sessions = unique(Aord.session_index, 'stable');
for i = 1:numel(sessions)
    idx = find(Aord.session_index == sessions(i));
    z = labelConf(idx);
    z = local_moving_mode(z, P.SmoothWindowAnchors);
    z = local_absorb_short_runs(z, P.MinRunAnchors);
    z = local_merge_same_label_gaps(z, P.MergeGapAnchors);
    smoothOrd(idx) = z;
end

anchorLabelTable = Aord;
anchorLabelTable.raw_label = rawOrd;
anchorLabelTable.smoothed_label = smoothOrd;
anchorLabelTable.max_posterior = maxPostOrd;
if ~isempty(postOrd)
    K = size(postOrd,2);
    for k = 1:K
        anchorLabelTable.(sprintf('posterior_%02d', k)) = postOrd(:,k);
    end
else
    K = max(rawOrd);
end

segTable = local_segments_from_anchors(anchorLabelTable, P);
segTable = local_absorb_short_duration_segments(segTable, P.MinSegmentDurationSec);
[frameLabelCell, framePosteriorCell] = local_project_to_frames(anchorLabelTable, segTable, K, P);

Segments = struct();
Segments.table = segTable;
Segments.frameLabelCell = frameLabelCell;
Segments.framePosteriorCell = framePosteriorCell;
Segments.anchorLabelTable = anchorLabelTable;
Segments.summary = local_global_summary(segTable, anchorLabelTable);
Segments.motifSummary = local_motif_summary(segTable, K);
Segments.sessionSummary = local_session_summary(segTable, anchorLabelTable);
Segments.params = rmfield(P, 'Cluster');
Segments.rawOrder = ord;
Segments.inverseRawOrder = invOrd;

if P.Verbose
    fprintf('build_motif_segments | anchors = %d | segments = %d | motifs = %d\n', ...
        N, height(Segments.table), K);
    if ~isempty(Segments.summary)
        disp(Segments.summary);
    end
end
end

function [A, labels, post, maxPost] = local_unpack_cluster(Cluster)
assert(isfield(Cluster,'Data') && isfield(Cluster.Data,'anchorTable') && istable(Cluster.Data.anchorTable), ...
    'Cluster.Data.anchorTable is required. Rerun cluster_multiscale_chunks after the Fig3 patch.');
A = Cluster.Data.anchorTable;
assert(ismember('session_index', A.Properties.VariableNames), 'anchorTable.session_index is required.');
assert(ismember('anchor_frame', A.Properties.VariableNames), 'anchorTable.anchor_frame is required.');
if ~ismember('anchor_time_s', A.Properties.VariableNames)
    A.anchor_time_s = nan(height(A),1);
end

assert(isfield(Cluster,'labels'), 'Cluster.labels is required.');
labels = double(Cluster.labels(:));

if isfield(Cluster,'posteriors') && ~isempty(Cluster.posteriors)
    post = double(Cluster.posteriors);
else
    K = max(labels);
    post = zeros(numel(labels), K);
    post(sub2ind(size(post), (1:numel(labels))', labels)) = 1;
end

if isfield(Cluster,'maxPosterior') && ~isempty(Cluster.maxPosterior)
    maxPost = double(Cluster.maxPosterior(:));
else
    maxPost = max(post, [], 2);
end
end

function zOut = local_moving_mode(z, win)
z = z(:); zOut = z;
if win <= 1 || numel(z) <= 2, return, end
h = floor(win/2);
for i = 1:numel(z)
    a = max(1, i-h); b = min(numel(z), i+h);
    zz = z(a:b); zz = zz(zz > 0);
    if isempty(zz), continue, end
    vals = unique(zz);
    cnt = arrayfun(@(v)sum(zz==v), vals);
    [~,ix] = max(cnt);
    zOut(i) = vals(ix);
end
end

function z = local_absorb_short_runs(z, minLen)
if minLen <= 1 || isempty(z), return, end
changed = true;
while changed
    changed = false;
    R = local_runs(z);
    for r = 1:height(R)
        if R.label(r) == 0 || R.length(r) >= minLen, continue, end
        left = r - 1; right = r + 1;
        leftOk = left >= 1 && R.label(left) > 0;
        rightOk = right <= height(R) && R.label(right) > 0;
        newLab = 0;
        if leftOk && rightOk && R.label(left) == R.label(right)
            newLab = R.label(left);
        elseif leftOk && rightOk
            if R.length(left) >= R.length(right), newLab = R.label(left); else, newLab = R.label(right); end
        elseif leftOk
            newLab = R.label(left);
        elseif rightOk
            newLab = R.label(right);
        end
        if newLab > 0
            z(R.start(r):R.stop(r)) = newLab;
            changed = true;
            break
        end
    end
end
end

function z = local_merge_same_label_gaps(z, maxGap)
if maxGap <= 0 || isempty(z), return, end
changed = true;
while changed
    changed = false;
    R = local_runs(z);
    for r = 2:height(R)-1
        if R.label(r) ~= 0 && R.length(r) > maxGap, continue, end
        if R.label(r-1) > 0 && R.label(r-1) == R.label(r+1) && R.length(r) <= maxGap
            z(R.start(r):R.stop(r)) = R.label(r-1);
            changed = true;
            break
        end
    end
end
end

function R = local_runs(z)
z = z(:);
if isempty(z)
    R = table([],[],[],[], 'VariableNames', {'label','start','stop','length'});
    return
end
b = [true; diff(z) ~= 0];
st = find(b);
sp = [st(2:end)-1; numel(z)];
lab = z(st);
len = sp - st + 1;
R = table(lab, st, sp, len, 'VariableNames', {'label','start','stop','length'});
end

function T = local_segments_from_anchors(A, P)
rows = struct([]); r = 0;
sessions = unique(A.session_index, 'stable');
for si = 1:numel(sessions)
    idx = find(A.session_index == sessions(si));
    if isempty(idx), continue, end
    lab = A.smoothed_label(idx);
    fr = A.anchor_frame(idx);
    tm = A.anchor_time_s(idx);
    mp = A.max_posterior(idx);
    R = local_runs(lab);
    for j = 1:height(R)
        if R.label(j) <= 0, continue, end
        ii = idx(R.start(j):R.stop(j));
        r = r + 1;
        rows(r).session_index = sessions(si); %#ok<AGROW>
        rows(r).motif_id = R.label(j);
        rows(r).start_anchor_index = ii(1);
        rows(r).stop_anchor_index = ii(end);
        rows(r).start_frame = fr(R.start(j));
        rows(r).stop_frame = fr(R.stop(j));
        rows(r).n_anchors = numel(ii);
        rows(r).median_posterior = median(mp(R.start(j):R.stop(j)), 'omitnan');
        if all(isfinite(tm([R.start(j), R.stop(j)])))
            rows(r).start_time_s = tm(R.start(j));
            rows(r).stop_time_s = tm(R.stop(j));
            rows(r).duration_s = max(0, rows(r).stop_time_s - rows(r).start_time_s);
        elseif ~isempty(P.FrameRate)
            rows(r).start_time_s = fr(R.start(j)) ./ P.FrameRate;
            rows(r).stop_time_s = fr(R.stop(j)) ./ P.FrameRate;
            rows(r).duration_s = max(0, rows(r).stop_time_s - rows(r).start_time_s);
        else
            rows(r).start_time_s = NaN;
            rows(r).stop_time_s = NaN;
            rows(r).duration_s = NaN;
        end
        rows(r).length_frames = rows(r).stop_frame - rows(r).start_frame + 1;
    end
end
if isempty(rows)
    T = table();
else
    T = struct2table(rows);
end
end

function T = local_absorb_short_duration_segments(T, minDur)
if isempty(T) || minDur <= 0 || ~ismember('duration_s', T.Properties.VariableNames), return, end
% Conservative final filter: mark very short segments but do not silently delete them.
T.is_short_duration = isfinite(T.duration_s) & T.duration_s < minDur;
end

function [labelCell, posteriorCell] = local_project_to_frames(A, S, K, P)
sessions = unique(A.session_index, 'stable');
labelCell = cell(max(sessions),1);
posteriorCell = cell(max(sessions),1);
for si = 1:numel(sessions)
    sess = sessions(si);
    nFrames = local_frame_count(sess, A, P);
    if isempty(nFrames) || nFrames <= 0, continue, end
    lab = zeros(nFrames,1);
    conf = nan(nFrames,1);
    idxS = find(S.session_index == sess);
    for j = idxS(:)'
        a = max(1, round(S.start_frame(j)));
        b = min(nFrames, round(S.stop_frame(j)));
        if b < a, continue, end
        lab(a:b) = S.motif_id(j);
        conf(a:b) = S.median_posterior(j);
    end
    labelCell{sess} = lab;
    posteriorCell{sess} = conf;
end
end

function nFrames = local_frame_count(sess, A, P)
nFrames = [];
F = P.FrameCountBySession;
if isempty(F)
    idx = A.session_index == sess;
    if any(idx), nFrames = max(A.anchor_frame(idx)); end
elseif isnumeric(F)
    if numel(F) >= sess, nFrames = F(sess); end
elseif istable(F)
    if all(ismember({'session_index','n_frames'}, F.Properties.VariableNames))
        ii = find(F.session_index == sess, 1);
        if ~isempty(ii), nFrames = F.n_frames(ii); end
    end
elseif isa(F,'containers.Map')
    key = num2str(sess);
    if isKey(F, key), nFrames = F(key); end
end
end

function T = local_global_summary(S, A)
if isempty(S)
    T = table(0, height(A), NaN, NaN, NaN, 'VariableNames', {'n_segments','n_anchors','median_segment_duration_s','median_segment_anchors','short_duration_fraction'});
    return
end
shortFrac = NaN;
if ismember('is_short_duration', S.Properties.VariableNames)
    shortFrac = mean(S.is_short_duration, 'omitnan');
end
T = table(height(S), height(A), median(S.duration_s,'omitnan'), median(S.n_anchors,'omitnan'), shortFrac, ...
    'VariableNames', {'n_segments','n_anchors','median_segment_duration_s','median_segment_anchors','short_duration_fraction'});
end

function T = local_motif_summary(S, K)
motif_id = (1:K)';
n_segments = zeros(K,1);
n_anchors = zeros(K,1);
median_duration_s = nan(K,1);
median_posterior = nan(K,1);
for k = 1:K
    idx = S.motif_id == k;
    n_segments(k) = sum(idx);
    if any(idx)
        n_anchors(k) = sum(S.n_anchors(idx));
        median_duration_s(k) = median(S.duration_s(idx), 'omitnan');
        median_posterior(k) = median(S.median_posterior(idx), 'omitnan');
    end
end
anchor_fraction = n_anchors ./ max(sum(n_anchors), eps);
T = table(motif_id, n_segments, n_anchors, anchor_fraction, median_duration_s, median_posterior);
end

function T = local_session_summary(S, A)
sessions = unique(A.session_index, 'stable');
session_index = sessions(:);
n_segments = zeros(numel(sessions),1);
n_anchors = zeros(numel(sessions),1);
median_duration_s = nan(numel(sessions),1);
for i = 1:numel(sessions)
    sess = sessions(i);
    n_anchors(i) = sum(A.session_index == sess);
    idx = S.session_index == sess;
    n_segments(i) = sum(idx);
    if any(idx), median_duration_s(i) = median(S.duration_s(idx), 'omitnan'); end
end
T = table(session_index, n_segments, n_anchors, median_duration_s);
end
