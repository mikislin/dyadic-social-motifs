function EmbedModel = fit_multiscale_chunk_embedding(ChunkSet, varargin)
%FIT_MULTISCALE_CHUNK_EMBEDDING Learn multi-scale chunk embeddings from chunk tensors.
%
% This function is the next layer after build_multiscale_chunk_dataset.
% It fits a compact chunk representation while preserving the current
% dyad variable conventions and scale structure.
%
% Strategy
% --------
% 1) For each scale separately, flatten chunk tensors [N x L x D] to [N x L*D].
% 2) Mask invalid timepoints using ChunkSet.scale(s).valid.
% 3) Robust-center and scale each flattened dimension across chunks.
% 4) Run PCA per scale.
% 5) Retain, for each scale, at least MinGlobalPCsPerScale PCs and enough PCs
%    to explain GlobalVarianceToKeep percent of per-scale variance.
% 6) Concatenate the retained per-scale scores across all scales.
% 7) Run a second-stage global PCA over the concatenated retained scores.
%
% Inputs
% ------
% ChunkSet : output of build_multiscale_chunk_dataset
%
% Name-value pairs
% ----------------
% 'nPCsPerScale' : max PCs retained per scale before global stage (default 12)
% 'globalNPCs' : max PCs retained after concatenation (default 16)
% 'minChunksPerScale' : minimum chunks required to fit a scale (default 50)
% 'winsorQuantiles' : [low high] quantiles before PCA (default [0.005 0.995])
% 'featureFamilyWeights': struct of optional weights by feature family
% 'booleanWeight' : additional multiplier for boolean channels (default 0.5)
% 'scaleWeightMode' : 'equal', 'sqrt_dim', or 'none' (default 'equal')
% 'globalVarianceToKeep' : cumulative per-scale variance threshold for the
%                          first-stage retention before the second PCA
%                          (default 96)
% 'minGlobalPCsPerScale' : minimum retained PCs per scale for the second PCA
%                          (default 10)
% 'verbose' : print progress (default true)
%
% Output
% ------
% EmbedModel.scale(s)
%   .chunkSec
%   .nInputDims
%   .nPCs
%   .nRetainedForGlobal
%   .coeff
%   .mu
%   .scale
%   .winsorLow
%   .winsorHigh
%   .explained
%   .score
%   .scoreForGlobal
%   .chunk_id
%   .meta
% EmbedModel.chunkTable
%   Combined metadata with unified embedding columns.
% EmbedModel.embedding
%   Final chunk embedding matrix after second-stage global PCA.
% EmbedModel.embeddingNames
%   Names of final embedding columns.
%
% See also: validate_multiscale_chunk_embedding, plot_multiscale_chunk_embedding_diagnostics

p = inputParser;
p.addParameter('nPCsPerScale', 12, @(x)isscalar(x) && x >= 1);
p.addParameter('globalNPCs', 16, @(x)isscalar(x) && x >= 1);
p.addParameter('minChunksPerScale', 50, @(x)isscalar(x) && x >= 1);
p.addParameter('winsorQuantiles', [0.005 0.995], @(x)isnumeric(x) && numel(x)==2 && x(1)>=0 && x(2)<=1 && x(1)<x(2));
p.addParameter('featureFamilyWeights', struct(), @isstruct);
p.addParameter('booleanWeight', 0.5, @(x)isscalar(x) && x > 0);
p.addParameter('scaleWeightMode', 'equal', @(x)ischar(x) || isstring(x));
p.addParameter('globalVarianceToKeep', 96, @(x)isscalar(x) && x > 0 && x <= 100);
p.addParameter('minGlobalPCsPerScale', 10, @(x)isscalar(x) && x >= 1);
p.addParameter('verbose', true, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});
P = p.Results;

scaleWeightMode = lower(string(P.scaleWeightMode));
assert(any(scaleWeightMode == ["equal","sqrt_dim","none"]), ...
    'scaleWeightMode must be ''equal'', ''sqrt_dim'', or ''none''.');

nScale = numel(ChunkSet.scale);
ScaleModel = repmat(struct( ...
    'chunkSec', [], ...
    'nInputDims', [], ...
    'nPCs', [], ...
    'nRetainedForGlobal', [], ...
    'coeff', [], ...
    'mu', [], ...
    'scale', [], ...
    'winsorLow', [], ...
    'winsorHigh', [], ...
    'explained', [], ...
    'score', [], ...
    'scoreForGlobal', [], ...
    'chunk_id', [], ...
    'meta', table(), ...
    'featureWeights', []), nScale, 1);

allTables = cell(nScale,1);
allChunkIds = cell(nScale,1);
allScores = cell(nScale,1);
allNames = cell(nScale,1);
allGlobalScores = cell(nScale,1);
allGlobalNames = cell(nScale,1);

for s = 1:nScale
    Sc = ChunkSet.scale(s);
    meta = Sc.meta;
    nChunks = size(Sc.X,1);
    L = size(Sc.X,2);
    D = size(Sc.X,3);

    if P.verbose
        fprintf('Embedding scale %d/%d (%.4gs) | nChunks = %d | L = %d | D = %d\n', ...
            s, nScale, Sc.chunkSec, nChunks, L, D);
    end

    ScaleModel(s).chunkSec = Sc.chunkSec;
    ScaleModel(s).nInputDims = L * D;

    if nChunks < P.minChunksPerScale
        warning('fit_multiscale_chunk_embedding:TooFewChunks', ...
            'Skipping scale %.4gs because nChunks=%d < minChunksPerScale=%d.', ...
            Sc.chunkSec, nChunks, P.minChunksPerScale);
        allTables{s} = table();
        allChunkIds{s} = zeros(0,1);
        allScores{s} = zeros(0,0);
        allNames{s} = strings(0,1);
        allGlobalScores{s} = zeros(0,0);
        allGlobalNames{s} = strings(0,1);
        continue
    end

    [Xflat, baseDimMeta] = flatten_chunk_tensor(Sc, ChunkSet);
    featureWeights = build_chunk_feature_weights(baseDimMeta, ChunkSet.channelMeta, P.featureFamilyWeights, P.booleanWeight);
    Xflat = Xflat .* featureWeights';
    [Xproc, prepStats] = preprocess_chunk_matrix(Xflat, 'winsorQuantiles', P.winsorQuantiles);

    maxPC = min([P.nPCsPerScale, size(Xproc,1)-1, size(Xproc,2)]);
    if maxPC < 1
        warning('fit_multiscale_chunk_embedding:NoPCs', ...
            'Skipping scale %.4gs because there are too few usable dimensions.', Sc.chunkSec);
        allTables{s} = table();
        allChunkIds{s} = zeros(0,1);
        allScores{s} = zeros(0,0);
        allNames{s} = strings(0,1);
        allGlobalScores{s} = zeros(0,0);
        allGlobalNames{s} = strings(0,1);
        continue
    end

    [coeff, score, ~, ~, explained] = pca(Xproc, 'NumComponents', maxPC);
    score = apply_scale_weight(score, scaleWeightMode);
    nKeepGlobal = choose_num_pcs_for_global(explained, P.globalVarianceToKeep, P.minGlobalPCsPerScale, maxPC);
    scoreForGlobal = score(:, 1:nKeepGlobal);

    ScaleModel(s).nPCs = maxPC;
    ScaleModel(s).nRetainedForGlobal = nKeepGlobal;
    ScaleModel(s).coeff = coeff;
    ScaleModel(s).mu = prepStats.mu;
    ScaleModel(s).scale = prepStats.scale;
    ScaleModel(s).winsorLow = prepStats.winsorLow;
    ScaleModel(s).winsorHigh = prepStats.winsorHigh;
    ScaleModel(s).explained = explained(:);
    ScaleModel(s).score = score;
    ScaleModel(s).scoreForGlobal = scoreForGlobal;
    ScaleModel(s).chunk_id = meta.chunk_id;
    ScaleModel(s).meta = meta;
    ScaleModel(s).featureWeights = featureWeights;

    varNames = compose('scale%02d_pc%02d', s, 1:maxPC)';
    varNamesGlobal = compose('scale%02d_gpc%02d', s, 1:nKeepGlobal)';

    T = meta;
    for j = 1:maxPC
        T.(varNames{j}) = score(:,j);
    end
    for j = 1:nKeepGlobal
        T.(varNamesGlobal{j}) = scoreForGlobal(:,j);
    end

    allTables{s} = T;
    allChunkIds{s} = meta.chunk_id;
    allScores{s} = score;
    allNames{s} = string(varNames(:));
    allGlobalScores{s} = scoreForGlobal;
    allGlobalNames{s} = string(varNamesGlobal(:));

    if P.verbose
        fprintf('  retained for global PCA: %d PCs (cum explained = %.2f%%)\n', ...
            nKeepGlobal, sum(explained(1:nKeepGlobal)));
    end
end

chunkTable = join_scale_embeddings(ChunkSet.chunkTable, allTables, [allNames; allGlobalNames]);
[embedding, embeddingNames, globalModel, globalStruct] = build_global_embedding(chunkTable, allGlobalNames, P.globalNPCs);

EmbedModel = struct();
EmbedModel.scale = ScaleModel;
EmbedModel.chunkTable = chunkTable;
EmbedModel.embedding = embedding;
EmbedModel.embeddingNames = embeddingNames;
EmbedModel.globalModel = globalModel;
EmbedModel.global = globalStruct; % compatibility with plotting utilities
EmbedModel.params = P;
end

function [Xflat, dimMeta] = flatten_chunk_tensor(Sc, ChunkSet)
% Flatten [N x L x D] chunk tensor to [N x (L*D)] and mask invalid entries.
X = double(Sc.X);
V = logical(Sc.valid);
[nChunks, L, D] = size(X);

for t = 1:L
    bad = ~V(:,t);
    if any(bad)
        X(bad,t,:) = NaN;
    end
end

Xflat = reshape(X, nChunks, L * D);
obsNames = string(ChunkSet.obsNames(:));
channelMeta = ChunkSet.channelMeta;
baseFeature = strings(L*D,1);
channelType = strings(L*D,1);
timeIndex = zeros(L*D,1);
obsName = strings(L*D,1);
family = strings(L*D,1);
idx = 0;

for t = 1:L
    for d = 1:D
        idx = idx + 1;
        obsName(idx) = obsNames(d);
        baseFeature(idx) = string(channelMeta.BaseFeature(d));
        channelType(idx) = string(channelMeta.ChannelType(d));
        row = strcmp(string(ChunkSet.featureNames), string(channelMeta.BaseFeature(d)));
        if any(row)
            family(idx) = string(ChunkSet.featureMeta.Family{find(row,1,'first')});
        else
            family(idx) = "unknown";
        end
        timeIndex(idx) = t;
    end
end

dimMeta = table(obsName, baseFeature, channelType, family, timeIndex, ...
    'VariableNames', {'ObsName','BaseFeature','ChannelType','Family','TimeIndex'});
end

function weights = build_chunk_feature_weights(dimMeta, ~, familyWeightsStruct, booleanWeight)
weights = ones(height(dimMeta),1);
for i = 1:height(dimMeta)
    if dimMeta.ChannelType(i) == "boolean"
        weights(i) = weights(i) * booleanWeight;
    end
    fam = matlab.lang.makeValidName(char(dimMeta.Family(i)));
    if isfield(familyWeightsStruct, fam)
        weights(i) = weights(i) * familyWeightsStruct.(fam);
    end
end
muw = mean(weights);
if isfinite(muw) && muw > 0
    weights = weights ./ muw;
end
end

function [Xproc, stats] = preprocess_chunk_matrix(X, varargin)
p = inputParser;
p.addParameter('winsorQuantiles', [0.005 0.995], @(x)isnumeric(x) && numel(x)==2);
p.parse(varargin{:});
P = p.Results;

Xproc = X;
D = size(X,2);
stats = struct();
stats.mu = zeros(1,D);
stats.scale = ones(1,D);
stats.winsorLow = nan(1,D);
stats.winsorHigh = nan(1,D);

for d = 1:D
    x = Xproc(:,d);
    ok = isfinite(x);
    if nnz(ok) < 5
        Xproc(:,d) = 0;
        continue
    end
    q = quantile(x(ok), P.winsorQuantiles);
    x(ok) = min(max(x(ok), q(1)), q(2));
    stats.winsorLow(d) = q(1);
    stats.winsorHigh(d) = q(2);
    med = median(x(ok));
    sc = iqr(x(ok));
    if ~(isfinite(sc) && sc > 0)
        sc = std(x(ok), 0);
    end
    if ~(isfinite(sc) && sc > 0)
        sc = 1;
    end
    x(~ok) = med;
    x = (x - med) ./ sc;
    Xproc(:,d) = x;
    stats.mu(d) = med;
    stats.scale(d) = sc;
end

keep = any(abs(Xproc) > 0, 1);
if any(keep)
    Xproc = Xproc(:, keep);
    stats.mu = stats.mu(keep);
    stats.scale = stats.scale(keep);
    stats.winsorLow = stats.winsorLow(keep);
    stats.winsorHigh = stats.winsorHigh(keep);
end
end

function score = apply_scale_weight(score, mode)
switch mode
    case "equal"
        if ~isempty(score)
            score = score ./ sqrt(size(score,2));
        end
    case "sqrt_dim"
        if ~isempty(score)
            score = score ./ (size(score,2) ^ 0.25);
        end
    case "none"
        % no-op
end
end

function nKeep = choose_num_pcs_for_global(explained, varToKeep, minPCs, maxPC)
if isempty(explained)
    nKeep = 0;
    return
end
cumExpl = cumsum(explained(:));
idx = find(cumExpl >= varToKeep, 1, 'first');
if isempty(idx)
    idx = numel(explained);
end
nKeep = max(minPCs, idx);
nKeep = min(nKeep, maxPC);
nKeep = min(nKeep, numel(explained));
end

function chunkTable = join_scale_embeddings(baseChunkTable, allTables, nameGroups)
chunkTable = baseChunkTable;
for s = 1:numel(allTables)
    Ts = allTables{s};
    if isempty(Ts)
        continue
    end
    keepNames = ["chunk_id"; string(nameGroups{s}(:))];
    keepNames = keepNames(ismember(keepNames, string(Ts.Properties.VariableNames)));
    TsSmall = Ts(:, cellstr(keepNames));
    chunkTable = outerjoin(chunkTable, TsSmall, 'Keys', 'chunk_id', 'MergeKeys', true, 'Type', 'left');
end
chunkTable = sortrows(chunkTable, 'chunk_id');
end

function [embedding, embeddingNames, globalModel, globalStruct] = build_global_embedding(chunkTable, allGlobalNames, globalNPCs)
scoreNames = strings(0,1);
for s = 1:numel(allGlobalNames)
    scoreNames = [scoreNames; allGlobalNames{s}(:)]; %#ok<AGROW>
end
scoreNames = unique(scoreNames, 'stable');

if isempty(scoreNames)
    embedding = zeros(height(chunkTable),0);
    embeddingNames = strings(0,1);
    globalModel = struct('coeff', [], 'explained', [], 'mu', [], 'score', [], 'scoreNames', scoreNames);
    globalStruct = struct('score', [], 'explained', [], 'meta', table());
    return
end

X = chunkTable{:, cellstr(scoreNames)};
for j = 1:size(X,2)
    x = X(:,j);
    med = median(x, 'omitnan');
    if ~isfinite(med)
        med = 0;
    end
    x(~isfinite(x)) = med;
    X(:,j) = x;
end

nPC = min([globalNPCs, size(X,1)-1, size(X,2)]);
if nPC < 1
    embedding = X;
    embeddingNames = scoreNames;
    globalModel = struct('coeff', [], 'explained', [], 'mu', zeros(1,size(X,2)), 'score', X, 'scoreNames', scoreNames);
    globalStruct = struct('score', X, 'explained', [], 'meta', chunkTable);
    return
end

[coeff, score, ~, ~, explained, mu] = pca(X, 'NumComponents', nPC);
embedding = score;
embeddingNames = compose('embed_pc%02d', 1:nPC)';

globalModel = struct();
globalModel.coeff = coeff;
globalModel.explained = explained(:);
globalModel.mu = mu;
globalModel.score = score;
globalModel.scoreNames = scoreNames;

globalStruct = struct();
globalStruct.score = score;
globalStruct.explained = explained(:);
globalStruct.meta = chunkTable;
end
