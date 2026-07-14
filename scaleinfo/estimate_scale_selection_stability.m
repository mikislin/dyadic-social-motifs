function Stability = estimate_scale_selection_stability(ScaleScore, varargin)
%ESTIMATE_SCALE_SELECTION_STABILITY Bootstrap operational scale selection.
%
% This is a condition-blind fixed-representation stability audit. It
% resamples sessions from the common anchor table, recomputes scale-quality
% metrics from already-fitted scale embeddings, and reruns the predefined
% temporal-role selector. Provenance labels are not used.

p = inputParser;
p.addParameter('nBootstraps', 40, @(x)isscalar(x) && x >= 1);
p.addParameter('nMicro', 3, @(x)isscalar(x) && x >= 1);
p.addParameter('nMotif', 4, @(x)isscalar(x) && x >= 1);
p.addParameter('nContext', 2, @(x)isscalar(x) && x >= 1);
p.addParameter('minLogGap', 0.12, @(x)isscalar(x) && x >= 0);
p.addParameter('stabilityThreshold', 0.70, @(x)isscalar(x) && x >= 0 && x <= 1);
p.addParameter('rngSeed', 1, @isscalar);
p.parse(varargin{:});
P = p.Results;

assert(isfield(ScaleScore, 'scaleTable') && istable(ScaleScore.scaleTable), ...
    'estimate_scale_selection_stability:MissingScaleTable', ...
    'ScaleScore.scaleTable is required.');
assert(isfield(ScaleScore, 'embeddingByScale') && iscell(ScaleScore.embeddingByScale), ...
    'estimate_scale_selection_stability:MissingEmbeddings', ...
    'ScaleScore.embeddingByScale is required.');
assert(isfield(ScaleScore, 'anchorTable') && istable(ScaleScore.anchorTable), ...
    'estimate_scale_selection_stability:MissingAnchorTable', ...
    'ScaleScore.anchorTable is required.');

rng(P.rngSeed);
T0 = ScaleScore.scaleTable;
A = ScaleScore.anchorTable;
sessions = unique(A.session_index(:), 'stable');
assert(~isempty(sessions), 'estimate_scale_selection_stability:NoSessions', ...
    'ScaleScore.anchorTable contains no sessions.');

nScale = height(T0);
nBoot = floor(P.nBootstraps);
selectedCount = zeros(nScale, 1);
selectedRank = nan(nScale, nBoot);
roleRank = nan(nScale, nBoot);
roleScore = nan(nScale, nBoot);

for b = 1:nBoot
    sampledSessions = sessions(randi(numel(sessions), numel(sessions), 1));
    rowIdx = local_bootstrap_rows(A, sampledSessions);
    Tb = local_bootstrap_scale_table(T0, ScaleScore.embeddingByScale, rowIdx);

    Sel = select_operational_timescales(struct('scaleTable', Tb), ...
        'nMicro', P.nMicro, ...
        'nMotif', P.nMotif, ...
        'nContext', P.nContext, ...
        'minLogGap', P.minLogGap);

    [tf, loc] = ismember(T0.scale_index, Sel.scale_index);
    selectedCount = selectedCount + double(tf);
    selectedRank(tf, b) = Sel.rank_within_role(loc(tf));

    [rankB, scoreB] = local_role_ranks_and_scores(Tb);
    roleRank(:, b) = rankB(:);
    roleScore(:, b) = scoreB(:);
end

origSelected = false(nScale, 1);
origRank = nan(nScale, 1);
if isfield(ScaleScore, 'selectedTable') && istable(ScaleScore.selectedTable) && ~isempty(ScaleScore.selectedTable)
    [tf, loc] = ismember(T0.scale_index, ScaleScore.selectedTable.scale_index);
    origSelected(tf) = true;
    origRank(tf) = ScaleScore.selectedTable.rank_within_role(loc(tf));
end

[origRoleRank, origRoleScore] = local_role_ranks_and_scores(T0);

Stability = T0(:, {'scale_index','chunk_sec','initial_band'});
Stability.selected_in_original_run = origSelected;
Stability.original_selection_rank_within_role = origRank;
Stability.original_role_rank = origRoleRank;
Stability.original_role_score = origRoleScore;
Stability.selection_count = selectedCount;
Stability.selection_frequency = selectedCount ./ nBoot;
Stability.median_selected_rank_when_selected = local_row_median(selectedRank);
Stability.median_role_rank = local_row_median(roleRank);
Stability.role_rank_p10 = local_row_prctile(roleRank, 10);
Stability.role_rank_p90 = local_row_prctile(roleRank, 90);
Stability.mean_bootstrap_role_score = mean(roleScore, 2, 'omitnan');
Stability.sd_bootstrap_role_score = std(roleScore, 0, 2, 'omitnan');
Stability.n_bootstraps = repmat(nBoot, nScale, 1);
Stability.stability_threshold = repmat(P.stabilityThreshold, nScale, 1);
Stability.passes_stability_threshold = Stability.selection_frequency >= P.stabilityThreshold;
Stability.bootstrap_unit = repmat("session_blocks_from_common_anchor_table", nScale, 1);
Stability.bootstrap_representation = repmat("fixed_multiresolution_pca_scores", nScale, 1);
Stability.labels_used_for_stability = repmat("none", nScale, 1);
Stability.condition_used_for_stability = false(nScale, 1);
Stability.arena_used_for_stability = false(nScale, 1);
end

function rowIdx = local_bootstrap_rows(A, sampledSessions)
rowIdx = zeros(0, 1);
for i = 1:numel(sampledSessions)
    rows = find(A.session_index == sampledSessions(i));
    rowIdx = [rowIdx; rows(:)]; %#ok<AGROW>
end
if isempty(rowIdx)
    rowIdx = (1:height(A))';
end
end

function T = local_bootstrap_scale_table(T0, embedByScale, rowIdx)
T = T0;
nScale = height(T0);
for s = 1:nScale
    E = embedByScale{s};
    coh = local_coherence_metrics(E.score(rowIdx, :), E.labels(rowIdx));
    tmp = local_temporal_metrics(E, rowIdx);
    T.silhouette_like(s) = coh.silhouette_like;
    T.within_cluster_dispersion(s) = coh.within_cluster_dispersion;
    T.between_cluster_separation(s) = coh.between_cluster_separation;
    T.lag1_embedding_corr(s) = tmp.lag1_embedding_corr;
    T.label_run_frames(s) = tmp.label_run_frames;
    T.transition_entropy(s) = tmp.transition_entropy;
end

for s = 1:nScale
    if s > 1
        [T.predict_short_r2(s), T.redundancy_short(s)] = ...
            local_pair_predict_redundancy(embedByScale{s}, embedByScale{s - 1}, rowIdx);
    end
    if s < nScale
        [T.predict_long_r2(s), T.redundancy_long(s)] = ...
            local_pair_predict_redundancy(embedByScale{s}, embedByScale{s + 1}, rowIdx);
    end
end
T = local_add_normalized_scores(T);
end

function coh = local_coherence_metrics(score, labels)
X = score(:, 1:min(4, size(score, 2)));
lab = labels(:);
coh = struct('silhouette_like', NaN, 'within_cluster_dispersion', NaN, 'between_cluster_separation', NaN);
if size(X, 1) < 20 || numel(unique(lab)) < 2
    return
end
mu = mean(X, 1, 'omitnan');
within = 0;
between = 0;
count = 0;
for k = unique(lab)'
    idx = lab == k;
    Xk = X(idx, :);
    if isempty(Xk)
        continue
    end
    ck = mean(Xk, 1, 'omitnan');
    within = within + mean(sum((Xk - ck) .^ 2, 2), 'omitnan');
    between = between + nnz(idx) * sum((ck - mu) .^ 2, 'omitnan');
    count = count + 1;
end
between = between / max(size(X, 1), 1);
within = within / max(count, 1);
coh.within_cluster_dispersion = within;
coh.between_cluster_separation = between;
coh.silhouette_like = (between - within) / max(between + within, eps);
end

function tmp = local_temporal_metrics(E, rowIdx)
X = E.score(:, 1:min(3, size(E.score, 2)));
lab = E.labels(:);
A = E.anchorTable;
sampledSessions = A.session_index(rowIdx);
sessions = unique(sampledSessions, 'stable');
lagRows = [];
runs = [];
transEnt = [];
for s = sessions(:)'
    idx = rowIdx(sampledSessions == s);
    [~, ord] = sort(A.anchor_time_s(idx));
    idx = idx(ord);
    if numel(idx) < 5
        continue
    end
    Xs = X(idx, :);
    labs = lab(idx);
    r = corr(Xs(1:end-1, 1), Xs(2:end, 1), 'Rows', 'pairwise');
    if isfinite(r)
        lagRows(end + 1, 1) = r; %#ok<AGROW>
    end
    runs(end + 1, 1) = mean(local_run_lengths(labs), 'omitnan'); %#ok<AGROW>
    transEnt(end + 1, 1) = local_transition_entropy(labs); %#ok<AGROW>
end
tmp = struct();
tmp.lag1_embedding_corr = mean(lagRows, 'omitnan');
tmp.label_run_frames = mean(runs, 'omitnan');
tmp.transition_entropy = mean(transEnt, 'omitnan');
end

function [r2mean, redundancy] = local_pair_predict_redundancy(Ea, Eb, rowIdx)
Xa = Ea.score(rowIdx, 1:min(4, size(Ea.score, 2)));
Xb = Eb.score(rowIdx, 1:min(4, size(Eb.score, 2)));
r2 = nan(1, size(Xb, 2));
for j = 1:size(Xb, 2)
    y = Xb(:, j);
    X = [ones(size(Xa, 1), 1), Xa];
    beta = local_ridge_solve(X, y);
    yhat = X * beta;
    ssRes = sum((y - yhat) .^ 2);
    ssTot = sum((y - mean(y, 'omitnan')) .^ 2);
    r2(j) = 1 - ssRes / max(ssTot, eps);
end
r2mean = mean(r2, 'omitnan');
cc = corr(Xa, Xb, 'Rows', 'pairwise');
redundancy = median(abs(cc(:)), 'omitnan');
end

function beta = local_ridge_solve(X, y)
beta = pinv(X) * y;
end

function rows = local_add_normalized_scores(rows)
rows.z_pc1 = local_z(rows.pc1_explained);
rows.z_cum5 = local_z(rows.cum5_explained);
rows.z_cum_embedding_pcs = local_z(rows.cum_embedding_pcs_explained);
rows.z_effdim = local_z(rows.effective_dim);
rows.z_sil = local_z(rows.silhouette_like);
rows.z_between = local_z(rows.between_cluster_separation);
rows.z_persist = local_z(rows.lag1_embedding_corr);
rows.z_run = local_z(rows.label_run_frames);
rows.z_pred_short = local_z(rows.predict_short_r2);
rows.z_pred_long = local_z(rows.predict_long_r2);
rows.z_red_short = local_z(rows.redundancy_short);
rows.z_red_long = local_z(rows.redundancy_long);
rows.z_transEntropy = local_z(rows.transition_entropy);

rows.composite_micro = ...
    0.35 * rows.z_pc1 + ...
    0.20 * rows.z_cum_embedding_pcs + ...
    0.15 * rows.z_sil + ...
    0.20 * (-rows.z_run) + ...
    0.10 * rows.z_pred_long;

rows.composite_motif = ...
    0.20 * rows.z_pc1 + ...
    0.15 * rows.z_cum_embedding_pcs + ...
    0.20 * rows.z_sil + ...
    0.15 * rows.z_between + ...
    0.15 * rows.z_persist + ...
    0.10 * rows.z_run + ...
    0.05 * rows.z_pred_short;

rows.composite_context = ...
    0.10 * rows.z_pc1 + ...
    0.10 * rows.z_cum_embedding_pcs + ...
    0.20 * rows.z_effdim + ...
    0.20 * rows.z_persist + ...
    0.20 * rows.z_run + ...
    0.10 * rows.z_pred_short + ...
    0.10 * (-rows.z_transEntropy);

redPenalty = mean([rows.z_red_short, rows.z_red_long], 2, 'omitnan');
redPenalty(~isfinite(redPenalty)) = 0;
rows.composite_global = max([rows.composite_micro, rows.composite_motif, rows.composite_context], [], 2) - 0.25 * redPenalty;
end

function [rankWithinRole, scoreWithinRole] = local_role_ranks_and_scores(T)
rankWithinRole = nan(height(T), 1);
scoreWithinRole = nan(height(T), 1);
bands = ["micro","motif","context"];
scoreNames = ["composite_micro","composite_motif","composite_context"];
for b = 1:numel(bands)
    idx = find(T.initial_band == bands(b));
    if isempty(idx)
        continue
    end
    score = T.(scoreNames(b))(idx);
    [~, ord] = sort(score, 'descend', 'MissingPlacement', 'last');
    rankWithinRole(idx(ord)) = (1:numel(idx))';
    scoreWithinRole(idx) = score;
end
end

function z = local_z(x)
x = double(x(:));
mu = mean(x, 'omitnan');
sd = std(x, 0, 'omitnan');
if ~(isfinite(sd) && sd > 0)
    z = zeros(size(x));
else
    z = (x - mu) ./ sd;
end
end

function rr = local_run_lengths(lab)
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

function H = local_transition_entropy(lab)
lab = lab(:);
if numel(lab) < 2
    H = NaN;
    return
end
u = unique(lab);
K = numel(u);
M = zeros(K, K);
for t = 1:(numel(lab) - 1)
    i = find(u == lab(t), 1);
    j = find(u == lab(t + 1), 1);
    M(i, j) = M(i, j) + 1;
end
P = M ./ max(sum(M, 2), eps);
Hrow = -sum(P .* log2(max(P, eps)), 2, 'omitnan');
H = mean(Hrow, 'omitnan');
end

function med = local_row_median(X)
med = median(X, 2, 'omitnan');
end

function q = local_row_prctile(X, pct)
q = nan(size(X, 1), 1);
for i = 1:size(X, 1)
    xi = X(i, :);
    xi = xi(isfinite(xi));
    if ~isempty(xi)
        q(i) = prctile(xi, pct);
    end
end
end
