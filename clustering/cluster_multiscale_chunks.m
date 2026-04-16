function Cluster = cluster_multiscale_chunks(ChunkSet, EmbedModel, selectedScales, varargin)
%CLUSTER_MULTISCALE_CHUNKS Cluster anchor-aligned multiscale chunk embeddings.
%
% Inputs
%   ChunkSet        : multiscale chunk dataset
%   EmbedModel      : output of fit_multiscale_chunk_embedding
%   selectedScales  : numeric vector or selected-scales table/struct
%
% Output
%   Cluster         : struct with integrated data, labels, posteriors,
%                     medoids, stability, and evaluation metadata

p = inputParser;
p.addParameter('ClusterMethod', 'gmm', @(x)ischar(x) || isstring(x));
p.addParameter('NumClusters', 6, @(x)isscalar(x) && x >= 2);
p.addParameter('NumPCs', 10, @(x)isscalar(x) && x >= 2);
p.addParameter('StabilityBootstrapN', 10, @(x)isscalar(x) && x >= 0);
p.addParameter('Verbose', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('RegularizationValue', 1e-5, @(x)isscalar(x) && x >= 0);
p.addParameter('Replicates', 5, @(x)isscalar(x) && x >= 1);
p.addParameter('RandomSeed', 1, @isscalar);

% Forwarded build_anchor_multiscale_matrix options
p.addParameter('AnchorSetMode', 'reference', @(x)ischar(x) || isstring(x));
p.addParameter('AnchorAlignment', 'nearest', @(x)ischar(x) || isstring(x));
p.addParameter('RequireAllSelectedScales', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('AddScalePresenceIndicators', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('ImputeMissingWithScaleMean', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('RetentionRule', 'auto', @(x)ischar(x) || isstring(x));
p.addParameter('MinScalesPresent', [], @(x)isempty(x) || (isscalar(x) && x >= 1));
p.addParameter('MaxMissingFrac', [], @(x)isempty(x) || (isscalar(x) && x >= 0 && x <= 1));
p.addParameter('ReferenceScaleSec', [], @(x)isempty(x) || isscalar(x));
p.addParameter('AnchorToleranceFrames', [], @(x)isempty(x) || (isscalar(x) && x >= 0));
p.parse(varargin{:});
P = p.Results;

rng(P.RandomSeed);

% Be robust to either old or new build_anchor argument naming
commonArgs = { ...
    'AnchorSetMode', P.AnchorSetMode, ...
    'AnchorAlignment', P.AnchorAlignment, ...
    'RequireAllSelectedScales', logical(P.RequireAllSelectedScales), ...
    'AddScalePresenceIndicators', logical(P.AddScalePresenceIndicators), ...
    'ImputeMissingWithScaleMean', logical(P.ImputeMissingWithScaleMean), ...
    'RetentionRule', P.RetentionRule, ...
    'MinScalesPresent', P.MinScalesPresent, ...
    'ReferenceScaleSec', P.ReferenceScaleSec, ...
    'AnchorToleranceFrames', P.AnchorToleranceFrames, ...
    'Verbose', P.Verbose};

try
    Data = build_anchor_multiscale_matrix(ChunkSet, EmbedModel, selectedScales, ...
        commonArgs{:}, ...
        'MaxMissingFrac', P.MaxMissingFrac);
catch ME
    if contains(ME.message, 'MaxMissingFrac') || contains(ME.message, 'recognized parameter')
        Data = build_anchor_multiscale_matrix(ChunkSet, EmbedModel, selectedScales, ...
            commonArgs{:}, ...
            'MaxMissingScaleFraction', P.MaxMissingFrac);
    else
        rethrow(ME);
    end
end

X = Data.X;
assert(~isempty(X) && size(X,1) >= P.NumClusters, ...
    'Too few retained anchors for requested NumClusters.');
assert(all(isfinite(X), 'all'), ...
    'cluster_multiscale_chunks:NonfiniteX', ...
    'Integrated matrix contains NaN/Inf values.');

nPCs = min([P.NumPCs, size(X,2), size(X,1)-1]);
assert(nPCs >= 1, 'cluster_multiscale_chunks:TooFewPCs', ...
    'Not enough retained anchors or dimensions for PCA.');

[coeff, score, latent, ~, explained, mu] = pca(X, 'NumComponents', nPCs);
Xpca = score(:,1:nPCs);

method = lower(string(P.ClusterMethod));
labels = [];
post = [];
maxPost = [];
model = [];

switch method
    case "gmm"
        gm = fitgmdist(Xpca, P.NumClusters, ...
            'RegularizationValue', P.RegularizationValue, ...
            'Replicates', P.Replicates, ...
            'Options', statset('MaxIter', 1000, 'Display', 'off'));
        labels = cluster(gm, Xpca);
        post = posterior(gm, Xpca);
        maxPost = max(post, [], 2);
        model = gm;

    case "kmeans"
        [labels, C] = kmeans(Xpca, P.NumClusters, ...
            'Replicates', P.Replicates, ...
            'MaxIter', 1000, ...
            'Display', 'off');
        post = local_hard_posteriors(labels, P.NumClusters);
        maxPost = ones(size(labels));
        model = struct('Centroids', C, 'Method', 'kmeans');

    otherwise
        error('cluster_multiscale_chunks:BadMethod', ...
            'Unsupported ClusterMethod: %s', P.ClusterMethod);
end

labels = labels(:);

occupancy = accumarray(labels, 1, [P.NumClusters 1], @sum, 0);
occupancyFrac = occupancy / max(sum(occupancy), eps);

medoidIdx = local_compute_medoids(Xpca, labels, P.NumClusters);
centroidsPCA = local_compute_centroids(Xpca, labels, P.NumClusters);

stability = local_bootstrap_stability(Xpca, labels, P, method);

Cluster = struct();
Cluster.Data = Data;
Cluster.X = X;
Cluster.Xpca = Xpca;
Cluster.PCA = struct('coeff', coeff, 'score', score, 'latent', latent, ...
    'explained', explained, 'mu', mu, 'nPCsUsed', nPCs);

Cluster.labels = labels;
Cluster.posteriors = post;
Cluster.maxPosterior = maxPost(:);
Cluster.model = model;

Cluster.NumClusters = P.NumClusters;
Cluster.ClusterMethod = char(method);
Cluster.occupancy = occupancy;
Cluster.occupancyFrac = occupancyFrac;
Cluster.medoidIdx = medoidIdx;
Cluster.centroidsPCA = centroidsPCA;
Cluster.medoids = local_build_medoid_table(Data, medoidIdx);
Cluster.stability = stability;
Cluster.params = P;

assert(isfield(Cluster, 'maxPosterior') && numel(Cluster.maxPosterior) == size(X,1), ...
    'cluster_multiscale_chunks:MissingMaxPosterior', ...
    'Cluster.maxPosterior missing from final output struct.');
assert(isfield(Cluster, 'posteriors') && isequal(size(Cluster.posteriors), [size(X,1), P.NumClusters]), ...
    'cluster_multiscale_chunks:MissingPosteriors', ...
    'Cluster.posteriors missing or wrong size.');

if P.Verbose
    medianMaxPost = median(Cluster.maxPosterior, 'omitnan');
    fprintf('cluster_multiscale_chunks | anchors = %d | PCs = %d | clusters = %d | median max posterior = %.3f\n', ...
        size(X,1), nPCs, P.NumClusters, medianMaxPost);

    fprintf('  cluster occupancy frac = ');
    fprintf('%.4f ', Cluster.occupancyFrac);
    fprintf('\n');

    medConf = accumarray(Cluster.labels, Cluster.maxPosterior, [P.NumClusters 1], @median, NaN);
    fprintf('  median posterior by cluster = ');
    fprintf('%.4f ', medConf);
    fprintf('\n');

    if isfield(stability, 'meanARI')
        fprintf('  stability mean ARI = %.4f\n', stability.meanARI);
    end
end

end

function post = local_hard_posteriors(labels, K)
N = numel(labels);
post = zeros(N, K);
idx = sub2ind([N K], (1:N)', labels(:));
post(idx) = 1;
end

function medoidIdx = local_compute_medoids(X, labels, K)
medoidIdx = nan(K,1);
for k = 1:K
    idx = find(labels == k);
    if isempty(idx)
        continue
    end
    Xk = X(idx,:);
    ctr = mean(Xk, 1, 'omitnan');
    d2 = sum((Xk - ctr).^2, 2);
    [~, ix] = min(d2);
    medoidIdx(k) = idx(ix);
end
end

function centroids = local_compute_centroids(X, labels, K)
centroids = nan(K, size(X,2));
for k = 1:K
    idx = labels == k;
    if any(idx)
        centroids(k,:) = mean(X(idx,:), 1, 'omitnan');
    end
end
end

function T = local_build_medoid_table(Data, medoidIdx)
K = numel(medoidIdx);
cluster = (1:K)';
anchorIndex = medoidIdx(:);
sessionIndex = nan(K,1);
anchorFrame = nan(K,1);
timeSec = nan(K,1);

assert(isfield(Data, 'anchorTable') && istable(Data.anchorTable), ...
    'cluster_multiscale_chunks:MissingAnchorTable', ...
    'Data.anchorTable is required.');

for k = 1:K
    if ~isnan(anchorIndex(k)) && anchorIndex(k) >= 1 && anchorIndex(k) <= height(Data.anchorTable)
        sessionIndex(k) = Data.anchorTable.session_index(anchorIndex(k));
        anchorFrame(k) = Data.anchorTable.anchor_frame(anchorIndex(k));
        if ismember('anchor_time_s', Data.anchorTable.Properties.VariableNames)
            timeSec(k) = Data.anchorTable.anchor_time_s(anchorIndex(k));
        end
    end
end

T = table(cluster, anchorIndex, sessionIndex, anchorFrame, timeSec, ...
    'VariableNames', {'Cluster','AnchorIndex','SessionIndex','AnchorFrame','AnchorTimeSec'});
end

function stability = local_bootstrap_stability(X, labelsRef, P, method)
stability = struct();
stability.ARI = nan(P.StabilityBootstrapN,1);

if P.StabilityBootstrapN <= 0
    stability.meanARI = NaN;
    return
end

N = size(X,1);
for b = 1:P.StabilityBootstrapN
    idx = randsample(N, N, true);
    Xb = X(idx,:);
    try
        switch method
            case "gmm"
                gm = fitgmdist(Xb, P.NumClusters, ...
                    'RegularizationValue', P.RegularizationValue, ...
                    'Replicates', 1, ...
                    'Options', statset('MaxIter', 500, 'Display', 'off'));
                lb = cluster(gm, Xb);
            otherwise
                lb = kmeans(Xb, P.NumClusters, ...
                    'Replicates', 1, ...
                    'MaxIter', 500, ...
                    'Display', 'off');
        end
        stability.ARI(b) = local_adjusted_rand_index(labelsRef(idx), lb);
    catch
        stability.ARI(b) = NaN;
    end
end
stability.meanARI = mean(stability.ARI, 'omitnan');
end

function ari = local_adjusted_rand_index(a, b)
a = a(:);
b = b(:);
[~,~,a] = unique(a);
[~,~,b] = unique(b);

ua = unique(a);
ub = unique(b);
cont = zeros(numel(ua), numel(ub));
for i = 1:numel(ua)
    for j = 1:numel(ub)
        cont(i,j) = sum(a == ua(i) & b == ub(j));
    end
end

n = sum(cont, 'all');
comb2 = @(x) x .* (x - 1) ./ 2;

index = sum(comb2(cont), 'all');
rowSum = sum(cont, 2);
colSum = sum(cont, 1);

expected = sum(comb2(rowSum)) * sum(comb2(colSum)) / max(comb2(n), eps);
maxIndex = 0.5 * (sum(comb2(rowSum)) + sum(comb2(colSum)));
ari = (index - expected) / max(maxIndex - expected, eps);
end