function Stability = estimate_pca_loading_stability(ScaleScore, varargin)
%ESTIMATE_PCA_LOADING_STABILITY Split-half stability of run_06 PCA loadings.
%
% The audit resamples chunks within each scale using no biological labels.
% It reports sign-invariant PC1 agreement and PCA subspace similarity.

p = inputParser;
p.addParameter('nSplits', 20, @(x)isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('nPCs', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1));
p.addParameter('stabilityThreshold', 0.75, @(x)isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
p.addParameter('rngSeed', 123, @(x)isnumeric(x) && isscalar(x));
p.parse(varargin{:});
P = p.Results;

assert(isstruct(ScaleScore) && isfield(ScaleScore, 'embeddingByScale'), ...
    'estimate_pca_loading_stability:BadScaleScore', ...
    'ScaleScore.embeddingByScale is required.');

rng(P.rngSeed);
nScale = numel(ScaleScore.embeddingByScale);
Stability = table();

for s = 1:nScale
    E = ScaleScore.embeddingByScale{s};
    X = double(E.Xflat);
    n = size(X, 1);
    if isempty(P.nPCs)
        nPcs = min(size(E.coeff, 2), 12);
    else
        nPcs = min(P.nPCs, size(E.coeff, 2));
    end
    nPcs = max(1, nPcs);

    subspace = nan(P.nSplits, 1);
    pc1corr = nan(P.nSplits, 1);
    for b = 1:P.nSplits
        ord = randperm(n);
        nHalf = floor(n / 2);
        a = ord(1:nHalf);
        c = ord((nHalf + 1):(2 * nHalf));
        if numel(a) <= nPcs + 1 || numel(c) <= nPcs + 1
            continue
        end
        [Ca, okA] = local_fit_coeff(X(a,:), nPcs);
        [Cb, okB] = local_fit_coeff(X(c,:), nPcs);
        if ~(okA && okB)
            continue
        end
        k = min([size(Ca,2), size(Cb,2), nPcs]);
        sv = svd(Ca(:,1:k)' * Cb(:,1:k));
        subspace(b) = mean(min(max(sv, 0), 1), 'omitnan');
        pc1corr(b) = abs(dot(Ca(:,1), Cb(:,1)));
    end

    one = table();
    one.scale_index = local_scale_index(E, s);
    one.chunk_sec = E.chunk_sec;
    one.initial_band = local_scale_band(ScaleScore, s);
    one.n_chunks = n;
    one.n_split_repeats = P.nSplits;
    one.n_pcs_compared = nPcs;
    one.median_subspace_similarity = median(subspace, 'omitnan');
    one.subspace_similarity_p10 = local_prctile(subspace, 10);
    one.subspace_similarity_p90 = local_prctile(subspace, 90);
    one.median_pc1_abs_correlation = median(pc1corr, 'omitnan');
    one.pc1_abs_correlation_p10 = local_prctile(pc1corr, 10);
    one.pc1_abs_correlation_p90 = local_prctile(pc1corr, 90);
    one.loading_stability_threshold = P.stabilityThreshold;
    one.passes_loading_stability_threshold = one.median_subspace_similarity >= P.stabilityThreshold;
    one.bootstrap_unit = "chunks_split_half_within_scale";
    one.labels_used_for_loading_stability = "none";
    one.arena_used_for_loading_stability = false;
    one.condition_used_for_loading_stability = false;
    Stability = [Stability; one]; %#ok<AGROW>
end
end

function [coeff, ok] = local_fit_coeff(X, nPcs)
[Xproc, okCols] = local_preprocess(X);
maxPC = min([sum(okCols), size(Xproc,1)-1, nPcs]);
if maxPC < 1
    coeff = zeros(size(X,2), 0);
    ok = false;
    return
end
warnState = warning('off', 'stats:pca:ColRankDefX');
cleanup = onCleanup(@() warning(warnState)); %#ok<NASGU>
[cSmall, ~, ~, ~, ~] = pca(Xproc(:, okCols), 'NumComponents', maxPC);
coeff = zeros(size(X,2), maxPC);
coeff(okCols, :) = cSmall(:, 1:maxPC);
ok = true;
end

function [Xproc, okCols] = local_preprocess(X)
D = size(X, 2);
Xproc = zeros(size(X));
okCols = false(1, D);
for j = 1:D
    x = X(:, j);
    ok = isfinite(x);
    if nnz(ok) < 5
        continue
    end
    q = quantile(x(ok), [0.005 0.995]);
    x(ok) = min(max(x(ok), q(1)), q(2));
    med = median(x(ok), 'omitnan');
    sc = iqr(x(ok));
    if ~(isfinite(sc) && sc > 1e-6)
        sc = std(x(ok), 0, 'omitnan');
    end
    if ~(isfinite(sc) && sc > 1e-6)
        continue
    end
    x(~ok) = med;
    Xproc(:, j) = (x - med) ./ sc;
    okCols(j) = true;
end
end

function idx = local_scale_index(E, fallback)
idx = fallback;
if isfield(E, 'dimensionAudit') && istable(E.dimensionAudit) && ...
        ismember('scale_index', E.dimensionAudit.Properties.VariableNames)
    idx = double(E.dimensionAudit.scale_index(1));
end
end

function band = local_scale_band(ScaleScore, s)
band = "";
if isfield(ScaleScore, 'scaleTable') && istable(ScaleScore.scaleTable) && ...
        ismember('initial_band', ScaleScore.scaleTable.Properties.VariableNames)
    band = string(ScaleScore.scaleTable.initial_band(s));
end
end

function value = local_prctile(x, p)
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = prctile(x, p);
end
end
