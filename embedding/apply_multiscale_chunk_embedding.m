function EmbedOut = apply_multiscale_chunk_embedding(ChunkSet, EmbedModel)
%APPLY_MULTISCALE_CHUNK_EMBEDDING Apply an existing embedding model to a new ChunkSet.
%
% Mirrors fit_multiscale_chunk_embedding but uses stored per-scale PCA models.

nScale = numel(EmbedModel.scale);
allTables = cell(nScale,1);
allNames = cell(nScale,1);

for s = 1:nScale
    ScModel = EmbedModel.scale(s);
    if ScModel.nPCs < 1 || isempty(ScModel.coeff)
        allTables{s} = table();
        allNames{s} = strings(0,1);
        continue
    end

    Sc = ChunkSet.scale(s);
    [Xflat, dimMeta] = flatten_chunk_tensor(Sc, ChunkSet); %#ok<ASGLU>
    Xflat = Xflat .* ScModel.featureWeights';
    Xproc = apply_preprocess_chunk_matrix(Xflat, ScModel);
    score = Xproc * ScModel.coeff;

    varNames = compose('scale%02d_pc%02d', s, 1:size(score,2))';
    T = Sc.meta;
    for j = 1:size(score,2)
        T.(varNames{j}) = score(:,j);
    end
    allTables{s} = T;
    allNames{s} = string(varNames(:));
end

chunkTable = join_scale_embeddings(ChunkSet.chunkTable, allTables, allNames);
[embedding, embeddingNames, ~] = build_global_embedding(chunkTable, allNames, size(EmbedModel.embedding,2));

EmbedOut = struct();
EmbedOut.chunkTable = chunkTable;
EmbedOut.embedding = embedding;
EmbedOut.embeddingNames = embeddingNames;
end

function Xproc = apply_preprocess_chunk_matrix(X, ScModel)
Xproc = X;
for d = 1:size(Xproc,2)
    x = Xproc(:,d);
    ok = isfinite(x);
    if d <= numel(ScModel.winsorLow) && isfinite(ScModel.winsorLow(d))
        ql = ScModel.winsorLow(d);
        qh = ScModel.winsorHigh(d);
        x(ok) = min(max(x(ok), ql), qh);
    end
    mu = 0; sc = 1;
    if d <= numel(ScModel.mu), mu = ScModel.mu(d); end
    if d <= numel(ScModel.scale), sc = ScModel.scale(d); end
    x(~ok) = mu;
    x = (x - mu) ./ max(sc, eps);
    Xproc(:,d) = x;
end
keep = size(Xproc,2) >= size(ScModel.coeff,1);
if keep
    Xproc = Xproc(:,1:size(ScModel.coeff,1));
end
end
