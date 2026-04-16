function Fig = plot_multiscale_chunk_embedding_diagnostics(EmbedModel, varargin)
%PLOT_MULTISCALE_CHUNK_EMBEDDING_DIAGNOSTICS
% Journal-style diagnostics for the chunk-embedding layer.
%
% This function is written to be robust to small variations in the
% EmbedModel struct fields across package iterations.
%
% Expected content (any reasonable subset is fine):
%   EmbedModel.scale(s).chunkSec
%   EmbedModel.scale(s).explained
%   EmbedModel.scale(s).score
%   EmbedModel.scale(s).meta
%   EmbedModel.global.score
%   EmbedModel.global.explained
%   EmbedModel.global.meta
%
% Example:
%   Fig = plot_multiscale_chunk_embedding_diagnostics(EmbedModel, ...
%       'maxScatterPoints', 10000, ...
%       'maxScalesForDensity', 4);

p = inputParser;
p.addParameter('maxScatterPoints', 12000, @(x)isscalar(x) && x >= 100);
p.addParameter('maxScalesForDensity', 4, @(x)isscalar(x) && x >= 1);
p.addParameter('rngSeed', 1, @isscalar);
p.parse(varargin{:});
P = p.Results;

rng(P.rngSeed);

Fig = struct();

[scaleSec, perScaleExplained, globalExplained, globalScore, globalMeta, sessionScaleCount] = ...
    local_unpack_embed_model(EmbedModel);

%% Figure 1: scree / variance overview
Fig.scree = figure('Color','w', 'Name','Chunk embedding scree', ...
    'Position',[60 60 1600 900]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% 1A: per-scale scree
nexttile;
hold on
for s = 1:numel(scaleSec)
    ex = perScaleExplained{s};
    if isempty(ex); continue; end
    nShow = min(12, numel(ex));
    plot(1:nShow, ex(1:nShow), 'o-', 'LineWidth', 1.3);
end
xlabel('PC');
ylabel('Explained variance (%)');
title('Per-scale scree');
legend(compose('%.1fs', scaleSec), 'Location','northeastoutside');
box off
grid on

% 1B: cumulative first 5 PCs by scale
nexttile;
cum5 = nan(numel(scaleSec),1);
for s = 1:numel(scaleSec)
    ex = perScaleExplained{s};
    if isempty(ex); continue; end
    cum5(s) = sum(ex(1:min(5,numel(ex))));
end
bar(categorical(cellstr(compose('%.1fs', scaleSec))), cum5);
ylabel('Cumulative explained (%)');
title('Cumulative variance in first 5 PCs');
box off

% 1C: PC1 concentration by scale
nexttile;
pc1 = nan(numel(scaleSec),1);
for s = 1:numel(scaleSec)
    ex = perScaleExplained{s};
    if isempty(ex); continue; end
    pc1(s) = ex(1);
end
bar(categorical(cellstr(compose('%.1fs', scaleSec))), pc1);
ylabel('PC1 explained (%)');
title('PC1 concentration by scale');
box off

% 1D: global scree
nexttile;
if ~isempty(globalExplained)
    nShow = min(12, numel(globalExplained));
    yyaxis left
    plot(1:nShow, globalExplained(1:nShow), 'o-', 'LineWidth', 1.5);
    ylabel('Explained (%)')

    yyaxis right
    plot(1:nShow, cumsum(globalExplained(1:nShow)), 's--', 'LineWidth', 1.5);
    ylabel('Cumulative explained (%)')
    xlabel('Global embedding PC')
    title('Global embedding scree')
    box off
    grid on
else
    axis off
    text(0.1, 0.5, 'No global explained variance found.', 'FontSize', 12);
end

%% Figure 2: global embedding colored by scale / session
Fig.embedding = figure('Color','w', 'Name','Global embedding overview', ...
    'Position',[80 80 1800 800]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

if ~isempty(globalScore) && size(globalScore,2) >= 2 && ~isempty(globalMeta)
    idx = local_subsample_rows(size(globalScore,1), P.maxScatterPoints);

    % 2A scale
    nexttile;
    scatter(globalScore(idx,1), globalScore(idx,2), 10, double(globalMeta.scale_index(idx)), 'filled');
    xlabel('Embed PC1');
    ylabel('Embed PC2');
    title('Global embedding colored by scale');
    box off
    grid on
    colorbar

    % 2B session
    nexttile;
    if ismember('session_index', globalMeta.Properties.VariableNames)
        c = double(globalMeta.session_index(idx));
    else
        c = ones(numel(idx),1);
    end
    scatter(globalScore(idx,1), globalScore(idx,2), 10, c, 'filled');
    xlabel('Embed PC1');
    ylabel('Embed PC2');
    title('Global embedding colored by session');
    box off
    grid on
    colorbar
else
    nexttile; axis off; text(0.1,0.5,'Global embedding unavailable','FontSize',12);
    nexttile; axis off
end

%% Figure 3: per-scale density fingerprints in global PC1/PC2
Fig.density = figure('Color','w', 'Name','Scale density fingerprints', ...
    'Position',[100 100 1800 1000]);

if ~isempty(globalScore) && size(globalScore,2) >= 2 && ~isempty(globalMeta)
    showScales = 1:min(P.maxScalesForDensity, numel(scaleSec));
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    for ii = 1:numel(showScales)
        s = showScales(ii);
        nexttile;
        idx = globalMeta.scale_index == s;
        if nnz(idx) < 20
            axis off
            title(sprintf('Scale %d (%.1fs)', s, scaleSec(s)));
            continue
        end
        x = globalScore(idx,1);
        y = globalScore(idx,2);
        [N,C] = hist3([x y], [60 60]); %#ok<HIST3>
        imagesc(C{1}, C{2}, N');
        axis xy
        colorbar
        xlabel('Embed PC1');
        ylabel('Embed PC2');
        title(sprintf('Scale %d (%.1fs)', s, scaleSec(s)));
    end
else
    tiledlayout(1,1); nexttile; axis off
    text(0.1,0.5,'Density view unavailable','FontSize',12);
end

%% Figure 4: chunk count by session and scale
Fig.counts = figure('Color','w', 'Name','Chunk count by session and scale', ...
    'Position',[120 120 1200 900]);

if ~isempty(sessionScaleCount)
    imagesc(sessionScaleCount);
    colorbar
    xlabel('Scale index');
    ylabel('Session index');
    title('Chunk count by session and scale');
    xticks(1:numel(scaleSec));
    xticklabels(string(1:numel(scaleSec)));
else
    axis off
    text(0.1,0.5,'Session-scale count table unavailable','FontSize',12);
end

end

% -------------------------------------------------------------------------
function [scaleSec, perScaleExplained, globalExplained, globalScore, globalMeta, sessionScaleCount] = ...
    local_unpack_embed_model(M)

scaleSec = [];
perScaleExplained = {};
globalExplained = [];
globalScore = [];
globalMeta = table();
sessionScaleCount = [];

% Per-scale content
if isfield(M, 'scale') && ~isempty(M.scale)
    S = M.scale;
    nS = numel(S);
    scaleSec = nan(nS,1);
    perScaleExplained = cell(nS,1);

    for s = 1:nS
        if isfield(S(s), 'chunkSec')
            scaleSec(s) = S(s).chunkSec;
        elseif isfield(S(s), 'chunk_sec')
            scaleSec(s) = S(s).chunk_sec;
        else
            scaleSec(s) = s;
        end

        if isfield(S(s), 'explained')
            perScaleExplained{s} = S(s).explained(:)';
        elseif isfield(S(s), 'pcaExplained')
            perScaleExplained{s} = S(s).pcaExplained(:)';
        else
            perScaleExplained{s} = [];
        end
    end
elseif isfield(M, 'scaleSummary') && istable(M.scaleSummary)
    T = M.scaleSummary;
    if ismember('chunk_sec', T.Properties.VariableNames)
        scaleSec = T.chunk_sec;
    elseif ismember('chunkSec', T.Properties.VariableNames)
        scaleSec = T.chunkSec;
    else
        scaleSec = (1:height(T))';
    end
    perScaleExplained = repmat({[]}, numel(scaleSec), 1);
end

% Global content
if isfield(M, 'global') && isstruct(M.global)
    G = M.global;
    if isfield(G, 'explained'); globalExplained = G.explained(:)'; end
    if isfield(G, 'score');     globalScore = G.score; end
    if isfield(G, 'meta');      globalMeta = G.meta; end
elseif isfield(M, 'globalExplained')
    globalExplained = M.globalExplained(:)';
    if isfield(M, 'globalScore'); globalScore = M.globalScore; end
    if isfield(M, 'globalMeta');  globalMeta = M.globalMeta; end
end

% If no globalMeta, try to assemble from per-scale meta
if isempty(globalMeta) && isfield(M, 'scale') && ~isempty(M.scale)
    rows = table();
    for s = 1:numel(M.scale)
        if isfield(M.scale(s), 'score') && ~isempty(M.scale(s).score)
            n = size(M.scale(s).score,1);
            if isfield(M.scale(s), 'meta') && istable(M.scale(s).meta)
                T = M.scale(s).meta;
            else
                T = table();
            end
            if ~ismember('scale_index', T.Properties.VariableNames)
                T.scale_index = repmat(s, height(T), 1);
            end
            if height(T) ~= n
                T = table((1:n)', repmat(s,n,1), 'VariableNames', {'row_index','scale_index'});
            end
            rows = [rows; T]; %#ok<AGROW>
        end
    end
    globalMeta = rows;
end

% Session x scale count matrix
if ~isempty(globalMeta) && ...
        ismember('session_index', globalMeta.Properties.VariableNames) && ...
        ismember('scale_index', globalMeta.Properties.VariableNames)
    sessions = unique(globalMeta.session_index, 'stable');
    scales = unique(globalMeta.scale_index, 'stable');
    sessionScaleCount = nan(numel(sessions), numel(scales));
    for i = 1:numel(sessions)
        for j = 1:numel(scales)
            sessionScaleCount(i,j) = sum(globalMeta.session_index == sessions(i) & ...
                                         globalMeta.scale_index == scales(j));
        end
    end
end

end

% -------------------------------------------------------------------------
function idx = local_subsample_rows(n, maxN)
if n <= maxN
    idx = 1:n;
else
    idx = sort(randperm(n, maxN));
end
end