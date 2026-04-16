function Report = validate_multiscale_chunk_embedding(EmbedModel, varargin)
%VALIDATE_MULTISCALE_CHUNK_EMBEDDING Quantitative QC for chunk embeddings.
%
% Reports per-scale PCA behavior, global embedding structure, redundancy,
% and session/scale balance. This is the checkpoint before clustering or
% segmental HSMM fitting.

p = inputParser;
p.addParameter('makePlots', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('verbose', true, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});
P = p.Results;

nScale = numel(EmbedModel.scale);
scaleSummary = table();

for s = 1:nScale
    Sm = EmbedModel.scale(s);
    nChunks = size(Sm.score,1);
    nPCs = size(Sm.score,2);
    if isempty(Sm.explained)
        cum5 = NaN;
        pc1 = NaN;
    else
        pc1 = Sm.explained(1);
        cum5 = sum(Sm.explained(1:min(5, numel(Sm.explained))));
    end

    corrMed = NaN;
    if nPCs >= 2
        R = corr(Sm.score(:,1:min(5,nPCs)), 'Rows', 'pairwise');
        corrMed = median(abs(R(triu(true(size(R)),1))), 'omitnan');
    end

    scaleSummary = [scaleSummary; table( ...
        s, Sm.chunkSec, nChunks, nPCs, pc1, cum5, corrMed, ...
        'VariableNames', {'scale_index','chunk_sec','n_chunks','n_pcs','pc1_explained','cum5_explained','median_abs_pc_corr'})]; %#ok<AGROW>
end

globalSummary = table();
if ~isempty(EmbedModel.embedding)
    X = EmbedModel.embedding;
    R = corr(X, 'Rows', 'pairwise');
    globalSummary = table( ...
        size(X,1), size(X,2), ...
        EmbedModel.globalModel.explained(1), ...
        sum(EmbedModel.globalModel.explained(1:min(5,numel(EmbedModel.globalModel.explained)))), ...
        median(abs(R(triu(true(size(R)),1))), 'omitnan'), ...
        'VariableNames', {'n_chunks','n_dims','pc1_explained','cum5_explained','median_abs_corr'});
end

sessionScaleSummary = summarize_embedding_balance(EmbedModel.chunkTable);

Report = struct();
Report.scaleSummary = scaleSummary;
Report.globalSummary = globalSummary;
Report.sessionScaleSummary = sessionScaleSummary;
Report.isReadyForClustering = ~isempty(EmbedModel.embedding) && size(EmbedModel.embedding,2) >= 4;

if P.verbose
    disp('=== Chunk embedding validation: per-scale summary ===');
    disp(scaleSummary)
    disp('=== Chunk embedding validation: global summary ===');
    disp(globalSummary)
    disp('=== Chunk embedding validation: session-scale balance ===');
    disp(sessionScaleSummary)
    fprintf('Ready for clustering / segment scoring: %d\n', Report.isReadyForClustering);
end

if P.makePlots
    plot_multiscale_chunk_embedding_diagnostics(EmbedModel, Report);
end
end

function T = summarize_embedding_balance(chunkTable)
if isempty(chunkTable)
    T = table();
    return
end
sessions = unique(chunkTable.session_index, 'stable');
scales = unique(chunkTable.scale_index, 'stable');
rows = table();
for i = 1:numel(sessions)
    for s = 1:numel(scales)
        idx = chunkTable.session_index == sessions(i) & chunkTable.scale_index == scales(s);
        rows = [rows; table(sessions(i), scales(s), sum(idx), ...
            'VariableNames', {'session_index','scale_index','n_chunks'})]; %#ok<AGROW>
    end
end
T = rows;
end
