function ScaleScore = score_multiscale_chunk_bank(ChunkSet, varargin)
%SCORE_MULTISCALE_CHUNK_BANK Score usefulness of a dense multi-scale chunk bank.
%
% This function quantifies how informative each temporal scale is for
% downstream discovery of dyadic social interaction structure. It produces
% a per-scale score table, aligned per-scale embeddings, and a recommended
% hierarchical grouping of scales into micro / motif / context bands.
%
% Required input
%   ChunkSet : output of build_multiscale_chunk_dataset
%
% Name-value pairs
%   'nPCs'               : per-scale PCA dimensions for scoring (default 6)
%   'nClusters'          : provisional kmeans clusters for coherence scoring (default 8)
%   'maxChunksPerScale'  : cap per-scale sample size for fitting (default 4000)
%   'featureNames'       : features to emphasize in trajectory coherence (default auto)
%   'rngSeed'            : RNG seed (default 1)
%   'verbose'            : print summaries (default true)
%
% Output struct fields
%   .scaleTable          : one row per scale with raw metrics and composite scores
%   .selectedTable       : selected hierarchical subset and assignments
%   .embeddingByScale    : per-scale aligned embeddings used for scoring
%   .anchorTable         : aligned anchor metadata common across all scales
%   .params              : run parameters
%
% Notes
%   - Scoring is computed on anchor-aligned chunks. This preserves scale
%     comparability and avoids conflating scale identity with scale utility.
%   - Metrics are intentionally multi-criterion: representation quality,
%     coherence, persistence, cross-scale predictive value, and redundancy.
%
% See also: select_operational_timescales, validate_scale_usefulness_scores,
%           plot_scale_usefulness_diagnostics

p = inputParser;
p.addParameter('nPCs', 6, @(x)isscalar(x) && x >= 2);
p.addParameter('nClusters', 8, @(x)isscalar(x) && x >= 2);
p.addParameter('maxChunksPerScale', 4000, @(x)isscalar(x) && x >= 200);
p.addParameter('featureNames', strings(0,1), @(x)isstring(x) || iscell(x));
p.addParameter('rngSeed', 1, @isscalar);
p.addParameter('verbose', true, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});
P = p.Results;

rng(P.rngSeed);

assert(isfield(ChunkSet, 'scale') && ~isempty(ChunkSet.scale), ...
    'ChunkSet.scale is required and must be non-empty.');

featureNames = string(P.featureNames);
if isempty(featureNames)
    featureNames = ["centroid_dist","body2body_dist","mutual_facing", ...
        "heading_diff_deg","radial_speed_12","tangential_speed_12", ...
        "approach_speed_1","approach_speed_2","in_contact","close_pair"];
end

% ---------------------------------------------------------------------
% Build anchor alignment table shared across scales.
% ---------------------------------------------------------------------
anchorTable = i_build_anchor_table(ChunkSet);
assert(~isempty(anchorTable), 'No common anchors found across scales.');

nScale = numel(ChunkSet.scale);
embedByScale = cell(nScale,1);

for s = 1:nScale
    Sc = ChunkSet.scale(s);
    embedByScale{s} = i_fit_scale_embedding(Sc, ChunkSet, anchorTable, P, featureNames);
end

% ---------------------------------------------------------------------
% Per-scale metrics.
% ---------------------------------------------------------------------
rows = table();
for s = 1:nScale
    E = embedByScale{s};
    rep = i_representation_metrics(E);
    coh = i_coherence_metrics(E, P.nClusters);
    tmp = i_temporal_metrics(E);

    [predShort, predLong, redShort, redLong] = i_crossscale_metrics(embedByScale, s);

    band = i_initial_band_from_scale(E.chunk_sec);
    rows = [rows; table( ...
        s, E.chunk_sec, E.chunk_frames, E.nAnchors, ...
        band, rep.pc1_explained, rep.cum5_explained, rep.effective_dim, ...
        coh.silhouette_like, coh.within_cluster_dispersion, coh.between_cluster_separation, ...
        tmp.lag1_embedding_corr, tmp.label_run_frames, tmp.transition_entropy, ...
        predShort, predLong, redShort, redLong, ...
        'VariableNames', { ...
        'scale_index','chunk_sec','chunk_frames','n_anchors', ...
        'initial_band','pc1_explained','cum5_explained','effective_dim', ...
        'silhouette_like','within_cluster_dispersion','between_cluster_separation', ...
        'lag1_embedding_corr','label_run_frames','transition_entropy', ...
        'predict_short_r2','predict_long_r2','redundancy_short','redundancy_long'})]; %#ok<AGROW>
end

% ---------------------------------------------------------------------
% Normalize metrics and compute composite role scores.
% ---------------------------------------------------------------------
rows = i_add_normalized_scores(rows);

ScaleScore = struct();
ScaleScore.scaleTable = rows;
ScaleScore.embeddingByScale = embedByScale;
ScaleScore.anchorTable = anchorTable;
ScaleScore.params = P;
ScaleScore.selectedTable = select_operational_timescales(ScaleScore);

if P.verbose
    disp('=== Scale usefulness: per-scale summary ===');
    disp(rows(:, {'scale_index','chunk_sec','initial_band','pc1_explained','cum5_explained', ...
        'silhouette_like','lag1_embedding_corr','label_run_frames', ...
        'predict_short_r2','predict_long_r2','composite_micro','composite_motif','composite_context'}));
    disp('=== Recommended operational subset ===');
    disp(ScaleScore.selectedTable)
end
end

function anchorTable = i_build_anchor_table(ChunkSet)
% One row per session-anchor approximately common across all scales.
%
% Exact anchor-frame intersection can be empty when different scales induce
% different left/right offsets modulo the fixed stride. To make alignment
% robust, build a master anchor-time table and match anchors within half a
% stride. This preserves scale comparability while avoiding brittle exact
% frame-index intersection.

% Estimate common stride in seconds.
strideList = [];
for s = 1:numel(ChunkSet.scale)
    meta = ChunkSet.scale(s).meta;
    for sess = unique(meta.session_index(:))'
        t = sort(meta.anchor_time_s(meta.session_index == sess));
        dt = diff(t);
        dt = dt(isfinite(dt) & dt > 0);
        if ~isempty(dt)
            strideList(end+1,1) = median(dt); %#ok<AGROW>
        end
    end
end
assert(~isempty(strideList), 'Could not estimate anchor stride from ChunkSet.scale meta.');
strideSec = median(strideList, 'omitnan');
tolSec = 0.51 * strideSec;

rows = [];
for sess = 1:ChunkSet.nSessions
    % Session-wise anchor times for each scale.
    scaleTimes = cell(numel(ChunkSet.scale),1);
    tMin = -inf;
    tMax = inf;
    for s = 1:numel(ChunkSet.scale)
        meta = ChunkSet.scale(s).meta;
        ts = sort(unique(meta.anchor_time_s(meta.session_index == sess)));
        scaleTimes{s} = ts(:);
        if isempty(ts)
            tMin = inf;
            tMax = -inf;
            break
        end
        tMin = max(tMin, ts(1));
        tMax = min(tMax, ts(end));
    end
    if ~(isfinite(tMin) && isfinite(tMax) && tMax >= tMin)
        continue
    end

    % Use the longest-scale anchors as a stable reference grid.
    refMeta = ChunkSet.scale(end).meta;
    refTimes = sort(unique(refMeta.anchor_time_s(refMeta.session_index == sess)));
    refTimes = refTimes(refTimes >= tMin & refTimes <= tMax);
    if isempty(refTimes)
        continue
    end

    keep = false(numel(refTimes),1);
    refFrame = nan(numel(refTimes),1);
    for i = 1:numel(refTimes)
        ti = refTimes(i);
        ok = true;
        for s = 1:numel(scaleTimes)
            ts = scaleTimes{s};
            if isempty(ts) || min(abs(ts - ti)) > tolSec
                ok = false;
                break
            end
        end
        if ok
            keep(i) = true;
            refFrame(i) = 1 + round(ti * ChunkSet.sessions{sess}.fps);
        end
    end

    if any(keep)
        T = table();
        T.session_index = repmat(sess, nnz(keep), 1);
        T.anchor_frame = refFrame(keep);
        T.anchor_time_s = refTimes(keep);
        rows = [rows; T]; %#ok<AGROW>
    end
end

anchorTable = rows;
if ~isempty(anchorTable)
    anchorTable = sortrows(unique(anchorTable, 'rows'), {'session_index','anchor_time_s'});
    anchorTable.anchor_id = (1:height(anchorTable))';
end
end

function E = i_fit_scale_embedding(Sc, ChunkSet, anchorTable, P, featureNames)
meta = Sc.meta;
idx = i_align_anchor_rows(meta, anchorTable, ChunkSet.sessions, Sc.chunkSec);

X = Sc.Xraw(idx,:,:);
valid = Sc.valid(idx,:);
Xflat = flatten_chunk_tensor(Sc, ChunkSet, X, valid, featureNames);

% Sample for model fit.
nA = size(Xflat,1);
if nA > P.maxChunksPerScale
    fitIdx = sort(randperm(nA, P.maxChunksPerScale));
else
    fitIdx = 1:nA;
end

[Xproc, prep] = i_preprocess_chunk_matrix(Xflat(fitIdx,:));
[coeff, scoreFit, latent, ~, explained, mu] = pca(Xproc, 'NumComponents', min([P.nPCs, size(Xproc,2), size(Xproc,1)-1])); %#ok<ASGLU>

XallProc = i_apply_preprocess_chunk_matrix(Xflat, prep);
scoreAll = (XallProc - mu) * coeff;
scoreAll = scoreAll(:, 1:min(P.nPCs, size(scoreAll,2)));

% Provisional clustering for coherence and run metrics.
K = min(P.nClusters, max(2, floor(sqrt(size(scoreAll,1)/40))));
labels = kmeans(scoreAll, K, 'Replicates', 3, 'MaxIter', 200, 'Display', 'off');

E = struct();
E.chunk_sec = Sc.chunkSec;
E.chunk_frames = Sc.nFrames;
E.anchorTable = anchorTable;
E.featureNames = featureNames;
E.Xflat = Xflat;
E.score = scoreAll;
E.coeff = coeff;
E.mu = mu;
E.explained = explained(:);
E.labels = labels(:);
E.preprocess = prep;
E.nAnchors = size(scoreAll,1);
end

function Xflat = flatten_chunk_tensor(Sc, ChunkSet, XrawOverride, validOverride, featureNames)
%FLATTEN_CHUNK_TENSOR Flatten selected chunk features into one row per chunk.
if nargin < 3 || isempty(XrawOverride)
    XrawOverride = Sc.Xraw;
end
if nargin < 4 || isempty(validOverride)
    validOverride = Sc.valid;
end
if nargin < 5 || isempty(featureNames)
    featureNames = string(ChunkSet.featureNames);
else
    featureNames = string(featureNames);
end

featIdx = find(ismember(string(ChunkSet.featureNames), featureNames));
assert(~isempty(featIdx), 'No requested featureNames were found in ChunkSet.featureNames.');

X = XrawOverride(:,:,featIdx);
valid = validOverride;
[n, L, d] = size(X);
Xflat = nan(n, L*d);
for i = 1:n
    Xi = squeeze(X(i,:,:));
    if size(Xi,1) ~= L
        Xi = Xi';
    end
    vi = valid(i,:).';
    Xi(~vi, :) = NaN;
    Xflat(i,:) = reshape(Xi.', 1, []);
end
end


function idx = i_align_anchor_rows(meta, anchorTable, sessions, chunkSec)
% Align anchorTable rows to the nearest chunk rows for this scale.
idx = nan(height(anchorTable),1);
for sess = unique(anchorTable.session_index(:))'
    qa = find(anchorTable.session_index == sess);
    qm = find(meta.session_index == sess);
    if isempty(qa) || isempty(qm)
        continue
    end
    tq = anchorTable.anchor_time_s(qa);
    tm = meta.anchor_time_s(qm);

    if numel(tm) >= 2
        dt = diff(sort(unique(tm)));
        dt = dt(isfinite(dt) & dt > 0);
        if isempty(dt)
            tolSec = 0.5 / sessions{sess}.fps;
        else
            tolSec = 0.51 * median(dt, 'omitnan');
        end
    else
        tolSec = 0.5 / sessions{sess}.fps;
    end

    for ii = 1:numel(qa)
        [d, j] = min(abs(tm - tq(ii)));
        if isfinite(d) && d <= tolSec
            idx(qa(ii)) = qm(j);
        end
    end
end
assert(all(isfinite(idx)), 'Anchor alignment failed for scale %.4f.', chunkSec);
idx = idx(:);
end

function [Xproc, prep] = i_preprocess_chunk_matrix(X)
% Winsorize, impute, robust-scale.
D = size(X,2);
Xw = X;
med = nan(1,D);
sc = nan(1,D);
for j = 1:D
    x = Xw(:,j);
    ok = isfinite(x);
    if nnz(ok) >= 10
        q = quantile(x(ok), [0.005 0.995]);
        x(ok) = min(max(x(ok), q(1)), q(2));
        med(j) = median(x(ok));
        sc(j) = iqr(x(ok));
        if ~(isfinite(sc(j)) && sc(j) > 0)
            sc(j) = std(x(ok), 0);
        end
        if ~(isfinite(sc(j)) && sc(j) > 0)
            sc(j) = 1;
        end
        x(~ok) = med(j);
        Xw(:,j) = (x - med(j)) ./ sc(j);
    else
        x(~ok) = 0;
        Xw(:,j) = x;
        med(j) = 0;
        sc(j) = 1;
    end
end
prep = struct('median', med, 'scale', sc);
Xproc = Xw;
end

function Xproc = i_apply_preprocess_chunk_matrix(X, prep)
Xproc = X;
for j = 1:size(X,2)
    x = Xproc(:,j);
    x(~isfinite(x)) = prep.median(j);
    Xproc(:,j) = (x - prep.median(j)) ./ max(prep.scale(j), eps);
end
end

function rep = i_representation_metrics(E)
expl = E.explained(:);
if isempty(expl)
    expl = zeros(1,1);
end
rep = struct();
rep.pc1_explained = expl(1);
rep.cum5_explained = sum(expl(1:min(5,numel(expl))));
l = expl / max(sum(expl), eps);
rep.effective_dim = 1 / max(sum(l.^2), eps);
end

function coh = i_coherence_metrics(E, K)
X = E.score(:,1:min(4,size(E.score,2)));
lab = E.labels(:);
coh = struct();
if size(X,1) < 20 || numel(unique(lab)) < 2
    coh.silhouette_like = NaN;
    coh.within_cluster_dispersion = NaN;
    coh.between_cluster_separation = NaN;
    return
end

mu = mean(X,1);
within = 0;
between = 0;
count = 0;
cent = zeros(max(lab), size(X,2));
for k = unique(lab)'
    idx = lab == k;
    Xk = X(idx,:);
    ck = mean(Xk,1);
    cent(k,:) = ck;
    within = within + mean(sum((Xk - ck).^2, 2));
    nk = sum(idx);
    between = between + nk * sum((ck - mu).^2);
    count = count + 1;
end
between = between / max(size(X,1),1);
within = within / max(count,1);
coh.within_cluster_dispersion = within;
coh.between_cluster_separation = between;
coh.silhouette_like = (between - within) / max(between + within, eps);
end

function tmp = i_temporal_metrics(E)
X = E.score(:,1:min(3,size(E.score,2)));
lab = E.labels(:);
A = E.anchorTable;
rows = [];
runs = [];
transEnt = [];
for s = unique(A.session_index(:))'
    idx = find(A.session_index == s);
    idx = idx(:);
    if numel(idx) < 5
        continue
    end
    Xs = X(idx,:);
    labs = lab(idx);
    r = corr(Xs(1:end-1,1), Xs(2:end,1), 'Rows', 'pairwise');
    if isfinite(r)
        rows(end+1,1) = r; %#ok<AGROW>
    end
    rr = i_run_lengths(labs);
    runs(end+1,1) = mean(rr); %#ok<AGROW>
    transEnt(end+1,1) = i_transition_entropy(labs); %#ok<AGROW>
end

tmp = struct();
tmp.lag1_embedding_corr = mean(rows, 'omitnan');
tmp.label_run_frames = mean(runs, 'omitnan');
tmp.transition_entropy = mean(transEnt, 'omitnan');
end

function [predShort, predLong, redShort, redLong] = i_crossscale_metrics(embedByScale, s)
predShort = NaN; predLong = NaN; redShort = NaN; redLong = NaN;
if s > 1
    [predShort, redShort] = i_pair_predict_redundancy(embedByScale{s}, embedByScale{s-1});
end
if s < numel(embedByScale)
    [predLong, redLong] = i_pair_predict_redundancy(embedByScale{s}, embedByScale{s+1});
end
end

function [r2mean, redundancy] = i_pair_predict_redundancy(Ea, Eb)
Xa = Ea.score(:,1:min(4,size(Ea.score,2)));
Xb = Eb.score(:,1:min(4,size(Eb.score,2)));
assert(size(Xa,1) == size(Xb,1), 'Aligned scale embeddings must share anchor count.');

r2 = nan(1,size(Xb,2));
for j = 1:size(Xb,2)
    y = Xb(:,j);
    X = [ones(size(Xa,1),1), Xa];
    beta = X \ y;
    yhat = X * beta;
    ssRes = sum((y - yhat).^2);
    ssTot = sum((y - mean(y)).^2);
    r2(j) = 1 - ssRes / max(ssTot, eps);
end
r2mean = mean(r2, 'omitnan');

cc = corr(Xa, Xb, 'Rows', 'pairwise');
redundancy = median(abs(cc(:)), 'omitnan');
end

function rows = i_add_normalized_scores(rows)
% Larger-is-better metrics.
rows.z_pc1 = i_z(rows.pc1_explained);
rows.z_cum5 = i_z(rows.cum5_explained);
rows.z_effdim = i_z(rows.effective_dim);
rows.z_sil = i_z(rows.silhouette_like);
rows.z_between = i_z(rows.between_cluster_separation);
rows.z_persist = i_z(rows.lag1_embedding_corr);
rows.z_run = i_z(rows.label_run_frames);
rows.z_pred_short = i_z(rows.predict_short_r2);
rows.z_pred_long = i_z(rows.predict_long_r2);
rows.z_red_short = i_z(rows.redundancy_short);
rows.z_red_long = i_z(rows.redundancy_long);
rows.z_transEntropy = i_z(rows.transition_entropy);

% Role-specific composites.
rows.composite_micro = ...
    0.35 * rows.z_pc1 + ...
    0.20 * rows.z_cum5 + ...
    0.15 * rows.z_sil + ...
    0.20 * (-rows.z_run) + ...
    0.10 * rows.z_pred_long;

rows.composite_motif = ...
    0.20 * rows.z_pc1 + ...
    0.15 * rows.z_cum5 + ...
    0.20 * rows.z_sil + ...
    0.15 * rows.z_between + ...
    0.15 * rows.z_persist + ...
    0.10 * rows.z_run + ...
    0.05 * rows.z_pred_short;

rows.composite_context = ...
    0.10 * rows.z_pc1 + ...
    0.15 * rows.z_cum5 + ...
    0.15 * rows.z_effdim + ...
    0.20 * rows.z_persist + ...
    0.20 * rows.z_run + ...
    0.10 * rows.z_pred_short + ...
    0.10 * (-rows.z_transEntropy);

% Global usefulness with redundancy penalty.
redPenalty = mean([rows.z_red_short, rows.z_red_long], 2, 'omitnan');
redPenalty(~isfinite(redPenalty)) = 0;
rows.composite_global = max([rows.composite_micro, rows.composite_motif, rows.composite_context], [], 2) - 0.25 * redPenalty;
end

function z = i_z(x)
x = double(x(:));
mu = mean(x, 'omitnan');
sd = std(x, 0, 'omitnan');
if ~(isfinite(sd) && sd > 0)
    z = zeros(size(x));
else
    z = (x - mu) ./ sd;
end
end

function band = i_initial_band_from_scale(scaleSec)
if scaleSec < 0.9
    band = "micro";
elseif scaleSec < 2.6
    band = "motif";
else
    band = "context";
end
end

function rr = i_run_lengths(lab)
lab = lab(:);
if isempty(lab)
    rr = NaN;
    return
end
chg = [true; diff(lab) ~= 0];
start = find(chg);
stop = [start(2:end)-1; numel(lab)];
rr = stop - start + 1;
end

function H = i_transition_entropy(lab)
lab = lab(:);
if numel(lab) < 2
    H = NaN;
    return
end
u = unique(lab);
K = numel(u);
M = zeros(K,K);
for t = 1:numel(lab)-1
    i = find(u == lab(t), 1);
    j = find(u == lab(t+1), 1);
    M(i,j) = M(i,j) + 1;
end
P = M ./ max(sum(M,2), eps);
Hrow = -sum(P .* log2(max(P, eps)), 2, 'omitnan');
H = mean(Hrow, 'omitnan');
end
