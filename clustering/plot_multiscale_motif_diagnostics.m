function Fig = plot_multiscale_motif_diagnostics(Cluster, varargin)
%PLOT_MULTISCALE_MOTIF_DIAGNOSTICS Publication-style diagnostics for motif discovery.
%
% Required input
%   Cluster : output of cluster_multiscale_chunks
%
% Name-value pairs
%   'ChunkSet'             : original multiscale chunk dataset (optional, enables raw traces)
%   'ExampleSession'       : session index to show in ethogram/raw trace panel (default 1)
%   'RawFeatureNames'      : cellstr/string of raw dyad features to plot
%   'EthogramSmoothFrames' : median-filter window on labels for ethogram (default 9)
%   'MedoidWindowSec'      : window around medoid anchor for trace plots (default 4)
%   'MaxClustersToShow'    : max number of clusters in medoid grid (default all)
%
% Output
%   Fig : struct of figure handles

p = inputParser;
p.addParameter('ChunkSet', [], @(x)isstruct(x) || isempty(x));
p.addParameter('ExampleSession', 1, @(x)isscalar(x) && x >= 1);
p.addParameter('RawFeatureNames', ...
    {'centroid_dist','mutual_facing','radial_speed_12','tangential_speed_12','in_contact'}, ...
    @(x)iscell(x) || isstring(x));
p.addParameter('EthogramSmoothFrames', 9, @(x)isscalar(x) && x >= 1);
p.addParameter('MedoidWindowSec', 4, @(x)isscalar(x) && x > 0);
p.addParameter('MaxClustersToShow', inf, @(x)isscalar(x) && x >= 1);
p.parse(varargin{:});
P = p.Results;

assert(isstruct(Cluster), 'Cluster must be a struct.');
assert(isfield(Cluster, 'labels'), 'Cluster.labels is required.');
assert(isfield(Cluster, 'NumClusters'), 'Cluster.NumClusters is required.');
assert(isfield(Cluster, 'Data'), 'Cluster.Data is required.');
assert(isfield(Cluster.Data, 'anchorTable'), 'Cluster.Data.anchorTable is required.');

labels = Cluster.labels(:);
K = Cluster.NumClusters;
A = Cluster.Data.anchorTable;
nAnchors = numel(labels);

assert(height(A) == nAnchors, ...
    'anchorTable height must match number of retained anchors.');

if ismember('session_index', A.Properties.VariableNames)
    sessionIdx = A.session_index(:);
elseif ismember('sessionIdx', A.Properties.VariableNames)
    sessionIdx = A.sessionIdx(:);
else
    error('anchorTable must contain session_index or sessionIdx.');
end

if ismember('anchor_frame', A.Properties.VariableNames)
    anchorFrame = A.anchor_frame(:);
elseif ismember('anchorFrame', A.Properties.VariableNames)
    anchorFrame = A.anchorFrame(:);
else
    error('anchorTable must contain anchor_frame or anchorFrame.');
end

if ismember('anchor_time_s', A.Properties.VariableNames)
    anchorTime = A.anchor_time_s(:);
elseif ismember('anchorTimeSec', A.Properties.VariableNames)
    anchorTime = A.anchorTimeSec(:);
else
    anchorTime = anchorFrame(:);
end

rawFeatureNames = string(P.RawFeatureNames(:));
ethColors = local_ethogram_colors(K);

Fig = struct();

%% Figure 1: overview
Fig.overview = figure('Color','w', 'Name','Motif discovery overview', ...
    'Position',[40 40 1600 900]);
tiledlayout(2,3, 'TileSpacing','compact', 'Padding','compact');

nexttile;
axis off
txt = {
    sprintf('Retained anchors: %d', nAnchors)
    sprintf('Clusters: %d', K)
    sprintf('Median max posterior: %.3f', local_get_median_max_post(Cluster))
    sprintf('Stability mean ARI: %.3f', local_get_stability_mean_ari(Cluster))
    sprintf('PCs used: %d', local_get_npcs_used(Cluster))
    };
text(0.02, 0.95, txt, 'VerticalAlignment','top', 'FontSize', 12);

nexttile;
bar(1:K, Cluster.occupancyFrac(:), 'FaceColor',[0.25 0.25 0.25]);
xlabel('Cluster'); ylabel('Occupancy fraction');
title('Cluster occupancy');
box off

nexttile;
if isfield(Cluster, 'maxPosterior') && ~isempty(Cluster.maxPosterior)
    medConf = accumarray(labels, Cluster.maxPosterior(:), [K 1], @median, NaN);
    bar(1:K, medConf, 'FaceColor',[0.25 0.25 0.25]);
    ylim([0 1.02]);
    xlabel('Cluster'); ylabel('Median posterior');
    title('Confidence by cluster');
    box off
else
    axis off
    text(0.1, 0.5, 'maxPosterior unavailable', 'FontSize', 12);
end

nexttile;
if isfield(Cluster, 'Xpca') && size(Cluster.Xpca,2) >= 2
    hold on
    for k = 1:K
        idx = labels == k;
        if any(idx)
            scatter(Cluster.Xpca(idx,1), Cluster.Xpca(idx,2), 8, ...
                'MarkerFaceColor', ethColors(k,:), ...
                'MarkerFaceAlpha', 0.35, ...
                'MarkerEdgeColor', 'none');
        end
    end
    xlabel('PC1'); ylabel('PC2');
    title('Integrated motif embedding');
    box off; grid on
else
    axis off
end

nexttile;
if isfield(Cluster, 'PCA') && isfield(Cluster.PCA, 'explained')
    ex = Cluster.PCA.explained(:);
    nShow = min(20, numel(ex));
    yyaxis left
    plot(1:nShow, ex(1:nShow), 'o-', 'LineWidth',1.5, 'MarkerSize',4);
    ylabel('Explained variance (%)');
    yyaxis right
    plot(1:nShow, cumsum(ex(1:nShow)), 's--', 'LineWidth',1.5, 'MarkerSize',4);
    ylabel('Cumulative explained (%)');
    xlabel('PC');
    title('PCA variance explained');
    box off
else
    axis off
end

nexttile;
if isfield(Cluster, 'stability') && isfield(Cluster.stability, 'ARI')
    histogram(Cluster.stability.ARI, 'BinMethod','sturges', 'FaceColor',[0.3 0.3 0.3]);
    xlabel('Bootstrap ARI'); ylabel('Count');
    title('Stability distribution');
    box off
else
    axis off
end

sgtitle('Motif discovery overview', 'FontWeight','bold', 'FontSize',16);

%% Figure 2: ethogram-style raster
Fig.ethogram = figure('Color','w', 'Name','Motif ethogram', ...
    'Position',[70 70 1700 700]);

sess = P.ExampleSession;
idxSess = find(sessionIdx == sess);
assert(~isempty(idxSess), 'Requested ExampleSession not found in Cluster.Data.anchorTable.');

[~, order] = sort(anchorFrame(idxSess));
idxSess = idxSess(order);

labsSess = labels(idxSess);
timeSess = anchorTime(idxSess);
frameSess = anchorFrame(idxSess);

labsSmooth = local_smooth_labels(labsSess, P.EthogramSmoothFrames);

tiledlayout(3,1, 'TileSpacing','compact', 'Padding','compact');

nexttile;
imagesc(timeSess(:)', 1, labsSess(:)');
set(gca, 'YDir', 'normal');
colormap(gca, ethColors);
clim([1 K]);
ylabel('Raw');
title(sprintf('Session %d motif ethogram', sess));
set(gca,'YTick',[]);
box off

nexttile;
imagesc(timeSess(:)', 1, labsSmooth(:)');
set(gca, 'YDir', 'normal');
colormap(gca, ethColors);
clim([1 K]);
ylabel('Smoothed');
title(sprintf('Session %d bout-smoothed ethogram', sess));
set(gca,'YTick',[]);
box off

nexttile;
hold on
for k = 1:K
    idx = labsSmooth == k;
    if any(idx)
        plot(timeSess(idx), k * ones(nnz(idx),1), '.', 'Color', ethColors(k,:), 'MarkerSize', 8);
    end
end
ylim([0.5 K+0.5]);
xlabel('Time (s)');
ylabel('Motif');
title('Smoothed motif assignments');
box off
grid on

%% Figure 3: raw feature traces for session (if ChunkSet available)
Fig.raw_session = figure('Color','w', 'Name','Raw feature traces by motif', ...
    'Position',[90 90 1700 1000]);

if isempty(P.ChunkSet)
    axis off
    text(0.05, 0.5, ['ChunkSet not provided. Re-run with ' ...
        '''ChunkSet'', ChunkSet to show raw multiscale traces.'], 'FontSize', 12);
else
    nFeat = numel(rawFeatureNames);
    tiledlayout(nFeat + 1, 1, 'TileSpacing','compact', 'Padding','compact');

    nexttile;
    imagesc(timeSess(:)', 1, labsSmooth(:)');
    set(gca, 'YDir', 'normal');
    colormap(gca, ethColors);
    clim([1 K]);
    ylabel('Motif');
    title(sprintf('Session %d ethogram with requested raw features', sess));
    set(gca,'YTick',[]);
    box off

    % Use reference scale if available
    refScaleIdx = local_choose_reference_scale(P.ChunkSet, Cluster);
    scaleObj = P.ChunkSet.scale(refScaleIdx);

    for f = 1:nFeat
        nexttile;
        hold on
        featName = rawFeatureNames(f);

        ok = false;
        [tTrace, xTrace] = local_extract_session_feature_trace(scaleObj, sess, featName);
        if ~isempty(tTrace) && ~isempty(xTrace)
            ok = true;
            local_add_ethogram_background(gca, timeSess, labsSmooth, ethColors);
            plot(tTrace, xTrace, 'k-', 'LineWidth', 1.2);
        end

        ylabel(strrep(featName,'_','\_'), 'Interpreter','tex');
        if f == nFeat
            xlabel('Time (s)');
        else
            set(gca,'XTickLabel',[]);
        end
        grid on
        box off

        if ~ok
            text(0.05, 0.5, sprintf('Could not recover feature %s from ChunkSet.', featName), ...
                'Units','normalized', 'FontSize',11);
        end
    end
end

%% Figure 4: medoid-centered multiscale raw exemplars
Fig.medoids = figure('Color','w', 'Name','Medoid-centered multiscale exemplars', ...
    'Position',[110 110 1800 1100]);

if isempty(P.ChunkSet)
    axis off
    text(0.05, 0.5, ['ChunkSet not provided. Re-run with ' ...
        '''ChunkSet'', ChunkSet to show medoid-centered raw exemplars.'], 'FontSize', 12);
else
    clustersToShow = 1:min(K, P.MaxClustersToShow);
    nFeat = min(numel(rawFeatureNames), 5);
    tiledlayout(numel(clustersToShow), nFeat, 'TileSpacing','compact', 'Padding','compact');

    scaleIdxAll = local_get_selected_scale_indices(P.ChunkSet, Cluster);

    for r = 1:numel(clustersToShow)
        k = clustersToShow(r);

        if isnan(Cluster.medoidIdx(k))
            for c = 1:nFeat
                nexttile; axis off
            end
            continue
        end

        rowIdx = Cluster.medoidIdx(k);
        sessK = sessionIdx(rowIdx);
        frameK = anchorFrame(rowIdx);

        for c = 1:nFeat
            nexttile;
            hold on
            featName = rawFeatureNames(c);
            plotted = false;

            for s = 1:numel(scaleIdxAll)
                sc = P.ChunkSet.scale(scaleIdxAll(s));
                [tt, xx] = local_extract_anchor_feature_trace(sc, sessK, frameK, featName);
                if ~isempty(tt) && ~isempty(xx)
                    plot(tt - tt(round(end/2)), xx, 'LineWidth', 1.0);
                    plotted = true;
                end
            end

            xline(0, 'k--', 'LineWidth', 0.8);

            if r == 1
                title(strrep(featName,'_','\_'), 'Interpreter','tex');
            end
            if c == 1
                ylabel(sprintf('Cluster %d', k));
            end
            if r == numel(clustersToShow)
                xlabel('Time relative to medoid anchor (s)');
            end

            grid on
            box off

            if ~plotted
                text(0.05, 0.5, 'Trace unavailable', 'Units','normalized', 'FontSize', 10);
            end
        end
    end
end

sgtitle('Medoid-centered multiscale raw feature exemplars', 'FontWeight','bold', 'FontSize',16);

end

function cols = local_ethogram_colors(K)
base = [ ...
    0.1216 0.4667 0.7059
    0.8392 0.1529 0.1569
    0.1725 0.6275 0.1725
    0.5804 0.4039 0.7412
    1.0000 0.4980 0.0549
    0.0902 0.7451 0.8118
    0.8902 0.4667 0.7608
    0.5490 0.3373 0.2941
    0.4980 0.4980 0.4980
    0.7373 0.7412 0.1333];
if K <= size(base,1)
    cols = base(1:K,:);
else
    cols = lines(K);
end
end

function y = local_smooth_labels(x, win)
x = x(:);
win = max(1, round(win));
if win <= 1
    y = x;
    return
end
pad = floor(win/2);
xp = [repmat(x(1), pad, 1); x; repmat(x(end), pad, 1)];
y = x;
for i = 1:numel(x)
    seg = xp(i:(i+win-1));
    y(i) = mode(seg);
end
end

function v = local_get_median_max_post(Cluster)
if isfield(Cluster, 'maxPosterior') && ~isempty(Cluster.maxPosterior)
    v = median(Cluster.maxPosterior, 'omitnan');
else
    v = NaN;
end
end

function v = local_get_stability_mean_ari(Cluster)
if isfield(Cluster, 'stability') && isfield(Cluster.stability, 'meanARI')
    v = Cluster.stability.meanARI;
else
    v = NaN;
end
end

function v = local_get_npcs_used(Cluster)
if isfield(Cluster, 'PCA') && isfield(Cluster.PCA, 'nPCsUsed')
    v = Cluster.PCA.nPCsUsed;
elseif isfield(Cluster, 'Xpca')
    v = size(Cluster.Xpca,2);
else
    v = NaN;
end
end

function local_add_ethogram_background(ax, t, labs, cmap)
axes(ax); %#ok<LAXES>
yl = ylim(ax);
change = [1; find(diff(labs(:)) ~= 0) + 1; numel(labs)+1];
for j = 1:numel(change)-1
    a = change(j);
    b = change(j+1)-1;
    patch(ax, [t(a) t(b) t(b) t(a)], [yl(1) yl(1) yl(2) yl(2)], ...
        cmap(labs(a),:), 'FaceAlpha', 0.08, 'EdgeColor', 'none');
end
uistack(findobj(ax,'Type','line'),'top');
end

function sIdx = local_choose_reference_scale(ChunkSet, Cluster)
if isfield(Cluster.Data, 'selectedScaleIdx') && ~isempty(Cluster.Data.selectedScaleIdx)
    sIdx = Cluster.Data.selectedScaleIdx(1);
elseif isfield(ChunkSet, 'scale') && ~isempty(ChunkSet.scale)
    sIdx = 1;
else
    error('Could not determine reference scale.');
end
end

function idx = local_get_selected_scale_indices(ChunkSet, Cluster)
if isfield(Cluster.Data, 'selectedScaleIdx') && ~isempty(Cluster.Data.selectedScaleIdx)
    idx = Cluster.Data.selectedScaleIdx(:)';
else
    idx = 1:numel(ChunkSet.scale);
end
idx = idx(isfinite(idx) & idx >= 1 & idx <= numel(ChunkSet.scale));
end

function [t, x] = local_extract_session_feature_trace(scaleObj, sess, featName)
t = [];
x = [];

if ~isfield(scaleObj, 'meta') || ~isfield(scaleObj, 'Xraw')
    return
end
M = scaleObj.meta;

sessCol = local_pick_var(M.Properties.VariableNames, {'session_index','sessionIdx'});
timeCol = local_pick_var(M.Properties.VariableNames, {'anchor_time_s','timeSec'});
featCol = [];

if isempty(sessCol) || isempty(timeCol)
    return
end

if isfield(scaleObj, 'featureNames')
    featCol = find(strcmp(string(scaleObj.featureNames), string(featName)), 1);
elseif isfield(scaleObj, 'channelMeta') && ismember('BaseFeature', scaleObj.channelMeta.Properties.VariableNames)
    featCol = find(strcmp(string(scaleObj.channelMeta.BaseFeature), string(featName)), 1);
end

if isempty(featCol)
    return
end

idx = find(M.(sessCol) == sess);
if isempty(idx)
    return
end

[tt, ord] = sort(M.(timeCol)(idx));
idx = idx(ord);

Xr = scaleObj.Xraw;
if ndims(Xr) ~= 3
    return
end

% use center frame of chunk
mid = round(size(Xr,2) / 2);
x = squeeze(Xr(idx, mid, featCol));
t = tt(:);
x = x(:);
end

function [t, x] = local_extract_anchor_feature_trace(scaleObj, sess, frameAnchor, featName)
t = [];
x = [];

if ~isfield(scaleObj, 'meta') || ~isfield(scaleObj, 'Xraw')
    return
end

M = scaleObj.meta;
sessCol = local_pick_var(M.Properties.VariableNames, {'session_index','sessionIdx'});
frameCol = local_pick_var(M.Properties.VariableNames, {'anchor_frame','anchorFrame'});

if isempty(sessCol) || isempty(frameCol)
    return
end

if isfield(scaleObj, 'featureNames')
    featCol = find(strcmp(string(scaleObj.featureNames), string(featName)), 1);
elseif isfield(scaleObj, 'channelMeta') && ismember('BaseFeature', scaleObj.channelMeta.Properties.VariableNames)
    featCol = find(strcmp(string(scaleObj.channelMeta.BaseFeature), string(featName)), 1);
else
    featCol = [];
end

if isempty(featCol)
    return
end

idxSess = find(M.(sessCol) == sess);
if isempty(idxSess)
    return
end

[~, ix] = min(abs(M.(frameCol)(idxSess) - frameAnchor));
row = idxSess(ix);

Xr = scaleObj.Xraw;
if ndims(Xr) ~= 3
    return
end

x = squeeze(Xr(row,:,featCol));
x = x(:);

fps = 1;
if isfield(scaleObj, 'fps') && ~isempty(scaleObj.fps)
    fps = scaleObj.fps;
elseif isfield(M, 'fps')
    fps = M.fps(row);
end

t = ((1:numel(x))' - ceil(numel(x)/2)) ./ fps;
end

function vn = local_pick_var(allVars, candidates)
vn = '';
for i = 1:numel(candidates)
    if ismember(candidates{i}, allVars)
        vn = candidates{i};
        return
    end
end
end