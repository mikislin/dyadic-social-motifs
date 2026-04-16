function selection = select_features_for_clustering(X, featureNames, opts)
%SELECT_FEATURES_FOR_CLUSTERING Seurat-inspired variable-feature selection.
%
% Steps
%   1. Filter features with too much missingness or too little variation.
%   2. Robust-scale each remaining feature.
%   3. Score variability using MAD and IQR.
%   4. Keep top variable features.
%   5. Prune highly correlated features.

arguments
    X double
    featureNames cell
    opts.maxMissingFrac (1,1) double {mustBeGreaterThanOrEqual(opts.maxMissingFrac,0),mustBeLessThanOrEqual(opts.maxMissingFrac,1)} = 0.10
    opts.minMAD (1,1) double {mustBeGreaterThanOrEqual(opts.minMAD,0)} = 1e-6
    opts.minIQR (1,1) double {mustBeGreaterThanOrEqual(opts.minIQR,0)} = 1e-6
    opts.maxKeep (1,1) double {mustBeInteger,mustBePositive} = 120
    opts.minKeep (1,1) double {mustBeInteger,mustBePositive} = 20
    opts.maxPairCorr (1,1) double {mustBePositive,mustBeLessThan(opts.maxPairCorr,1.01)} = 0.90
end

assert(size(X,2) == numel(featureNames), 'featureNames must match X columns.');

p = size(X,2);
missingFrac = mean(isnan(X), 1);
MAD = mad(X, 1, 1);
IQRv = iqr(X);

baseKeep = missingFrac <= opts.maxMissingFrac & MAD > opts.minMAD & IQRv > opts.minIQR;
baseIdx = find(baseKeep);

scores = -inf(1,p);
Xscaled = nan(size(X));
for j = baseIdx
    mu = median(X(:,j), 'omitnan');
    sc = mad(X(:,j), 1);
    if isempty(sc) || sc == 0 || isnan(sc)
        continue;
    end
    Xscaled(:,j) = (X(:,j) - mu) ./ sc;
    scores(j) = mad(Xscaled(:,j), 1) + iqr(Xscaled(:,j));
end

[~, order] = sort(scores, 'descend', 'MissingPlacement', 'last');
order = order(isfinite(scores(order)));

nCandidate = min(max(opts.minKeep, sum(baseKeep)), opts.maxKeep);
candidates = order(1:min(nCandidate, numel(order)));

% Correlation pruning, keeping higher-scoring features first.
keep = false(1,p);
chosen = [];
for j = candidates
    if isempty(chosen)
        keep(j) = true;
        chosen(end+1) = j; %#ok<AGROW>
        continue;
    end
    r = corr(X(:,j), X(:,chosen), 'Type', 'Spearman', 'Rows', 'pairwise');
    if all(abs(r) < opts.maxPairCorr | isnan(r))
        keep(j) = true;
        chosen(end+1) = j; %#ok<AGROW>
    end
end

selection = struct();
selection.keepMask = keep;
selection.keepIdx = find(keep);
selection.keepNames = featureNames(keep);
selection.scores = scores(:);
selection.table = table(featureNames(:), missingFrac(:), MAD(:), IQRv(:), scores(:), keep(:), ...
    'VariableNames', {'Feature','MissingFrac','MAD','IQR','VariableScore','Keep'});
selection.Xselected = X(:, keep);
end
