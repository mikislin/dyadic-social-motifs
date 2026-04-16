function Fig = plot_motif_medoid_exemplars(Cluster, varargin)
%PLOT_MOTIF_MEDOID_EXEMPLARS Inspect medoid anchors in raw dyad feature space.

p = inputParser;
p.addParameter('MotifsToShow', [], @(x)isempty(x) || isnumeric(x));
p.addParameter('RawFeatureNames', {'centroid_dist','mutual_facing','radial_speed_12','tangential_speed_12'}, @(x)iscell(x) || isstring(x));
p.addParameter('Palette', [], @(x)isnumeric(x) || isempty(x));
p.parse(varargin{:});
P = p.Results;

K = Cluster.NumClusters;
if isempty(P.MotifsToShow)
    motifs = 1:min(K,6);
else
    motifs = P.MotifsToShow(:)';
end
if isempty(P.Palette)
    P.Palette = lines(max(motifs));
end
featNames = string(P.RawFeatureNames);
selSec = Cluster.Data.selectedScaleSec(:)';

nRows = numel(motifs);
nCols = numel(featNames);
Fig = figure('Color','w', 'Name','Integrated medoid exemplars', 'Position',[140 140 1800 max(900, 280*nRows)]);
tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact');

for r = 1:nRows
    k = motifs(r);
    med = Cluster.medoids(k);
    for c = 1:nCols
        nexttile;
        hold on
        feat = featNames(c);
        plotted = false;
        for s = 1:numel(selSec)
            [tt, x] = local_extract_medoid_trace(Cluster, med.session_index, med.anchor_frame, Cluster.Data.selectedScaleIdx(s), feat);
            if isempty(x)
                continue
            end
            plot(tt, x, 'LineWidth', 1.2, 'Color', P.Palette(k,:) * (0.65 + 0.35*s/numel(selSec)));
            plotted = true;
        end
        if ~plotted
            axis off
        else
            yline(0, '-', 'Color', [0.75 0.75 0.75]);
            box off
            grid on
        end
        if r == 1
            title(strrep(feat, '_', '\_'), 'Interpreter','tex');
        end
        if c == 1
            ylabel(sprintf('Motif %d\n sess %d\n frame %d', k, med.session_index, med.anchor_frame));
        end
        if r == nRows
            xlabel('Relative frame in chunk');
        end
    end
end
sgtitle('Integrated medoid signatures across selected scales', 'FontWeight','bold', 'FontSize', 16);
end

function [tt, x] = local_extract_medoid_trace(Cluster, sessionIdx, anchorFrame, scaleIndex, featName)
tt = [];
x = [];
sc = Cluster.ChunkSet.scale(scaleIndex);
if ~isfield(sc, 'meta') || ~isfield(sc, 'Xraw')
    return
end
meta = sc.meta;
if ~istable(meta)
    meta = struct2table(meta);
end
vars = string(meta.Properties.VariableNames);
if ismember('session_index', vars)
    sess = meta.session_index;
elseif ismember('sessionIdx', vars)
    sess = meta.sessionIdx;
else
    return
end
if ismember('anchor_frame', vars)
    fr = meta.anchor_frame;
elseif ismember('anchorFrame', vars)
    fr = meta.anchorFrame;
else
    return
end
row = find(double(sess) == double(sessionIdx) & double(fr) == double(anchorFrame), 1);
if isempty(row)
    return
end

col = local_find_feature_column(sc, featName);
if isempty(col)
    return
end
chunk = squeeze(sc.Xraw(row,:,:));
if isempty(chunk)
    return
end
x = double(chunk(:,col));
tt = (1:numel(x))';
end

function col = local_find_feature_column(sc, featName)
col = [];
if isfield(sc, 'channelMeta') && istable(sc.channelMeta)
    cm = sc.channelMeta;
    if ismember('BaseFeature', cm.Properties.VariableNames)
        col = find(strcmp(string(cm.BaseFeature), string(featName)), 1);
    elseif ismember('ObsName', cm.Properties.VariableNames)
        col = find(strcmp(string(cm.ObsName), string(featName)), 1);
    end
elseif isfield(sc, 'featureNames')
    col = find(strcmp(string(sc.featureNames), string(featName)), 1);
end
end
