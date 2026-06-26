function [embedding, embeddingNames, globalModel, globalStruct] = ...
    build_global_embedding(chunkTable, allGlobalNames, globalNPCs)
%BUILD_GLOBAL_EMBEDDING Fit a second-stage PCA over per-scale scores.

scoreNames = strings(0,1);
for s = 1:numel(allGlobalNames)
    scoreNames = [scoreNames; allGlobalNames{s}(:)]; %#ok<AGROW>
end
scoreNames = unique(scoreNames, 'stable');

if isempty(scoreNames)
    embedding = zeros(height(chunkTable),0);
    embeddingNames = strings(0,1);
    globalModel = struct('coeff', [], 'explained', [], 'mu', [], ...
        'score', [], 'scoreNames', scoreNames);
    globalStruct = struct('score', [], 'explained', [], 'meta', table());
    return
end

present = ismember(scoreNames, string(chunkTable.Properties.VariableNames));
if any(~present)
    missingNames = scoreNames(~present);
    warning('build_global_embedding:MissingGlobalScores', ...
        'Dropping %d missing global score columns before second-stage PCA. First missing column: %s', ...
        numel(missingNames), missingNames(1));
    scoreNames = scoreNames(present);
end
if isempty(scoreNames)
    embedding = zeros(height(chunkTable),0);
    embeddingNames = strings(0,1);
    globalModel = struct('coeff', [], 'explained', [], 'mu', zeros(1,size(chunkTable,2)), ...
        'score', [], 'scoreNames', scoreNames);
    globalStruct = struct('score', [], 'explained', [], 'meta', chunkTable);
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
    globalModel = struct('coeff', [], 'explained', [], ...
        'mu', zeros(1,size(X,2)), 'score', X, 'scoreNames', scoreNames);
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
globalModel.inputNames = scoreNames;
globalModel.embeddingNames = embeddingNames;

globalStruct = struct();
globalStruct.score = score;
globalStruct.explained = explained(:);
globalStruct.meta = chunkTable;
end
