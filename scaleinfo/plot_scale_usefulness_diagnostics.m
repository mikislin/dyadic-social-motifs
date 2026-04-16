function Fig = plot_scale_usefulness_diagnostics(ScaleScore, ChunkSet, varargin)
%PLOT_SCALE_USEFULNESS_DIAGNOSTICS Journal-style diagnostics for scale scoring.
%
% Inputs
%   ScaleScore : output of score_multiscale_chunk_bank
%   ChunkSet   : output of build_multiscale_chunk_dataset
%
% Name-value pairs
%   'featureNames'      : features to show in exemplar trajectories
%   'nExampleChunks'    : number of chunk exemplars per selected scale
%   'invertLog1p'       : true/false, invert log1p for distance-like features
%   'sessionIndex'      : session to use for exemplars
%   'rngSeed'           : RNG seed for exemplar sampling
%
% Output
%   Fig struct with figure handles

p = inputParser;
p.addParameter('featureNames', ["centroid_dist","body2body_dist","mutual_facing","heading_diff_deg"], ...
    @(x)isstring(x) || iscell(x));
p.addParameter('nExampleChunks', 5, @(x)isscalar(x) && x >= 1);
p.addParameter('invertLog1p', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('sessionIndex', 1, @(x)isscalar(x) && x >= 1);
p.addParameter('rngSeed', 1, @isscalar);
p.parse(varargin{:});
P = p.Results;

rng(P.rngSeed);
featNames = string(P.featureNames);

ST = ScaleScore.scaleTable;
Sel = ScaleScore.selectedTable;

Fig = struct();

%% Figure 1: score overview
Fig.overview = figure('Color','w', 'Name','Scale-usefulness scoring overview', ...
    'Position',[60 60 1700 950]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

% Representation
nexttile;
hold on
plot(log10(ST.chunk_sec), ST.pc1_explained, 'o-', 'LineWidth', 1.5);
plot(log10(ST.chunk_sec), ST.cum5_explained, 's-', 'LineWidth', 1.5);
xlabel('log10 scale (s)');
ylabel('Explained variance (%)');
title('Representation quality across scales');
legend({'PC1','First 5 PCs'}, 'Location','best');
box off; grid on

% Behavioral coherence
nexttile;
hold on
leg = strings(0,1);
if ismember('silhouette_like', ST.Properties.VariableNames)
    plot(log10(ST.chunk_sec), ST.silhouette_like, 'o-', 'LineWidth', 1.5);
    leg(end+1) = "Silhouette-like"; %#ok<AGROW>
end
if ismember('between_cluster_separation', ST.Properties.VariableNames)
    plot(log10(ST.chunk_sec), ST.between_cluster_separation, 's-', 'LineWidth', 1.5);
    leg(end+1) = "Between-cluster separation"; %#ok<AGROW>
end
xlabel('log10 scale (s)');
ylabel('Coherence metric');
title('Behavioral coherence across scales');
if ~isempty(leg)
    legend(cellstr(leg), 'Location','best');
end
box off; grid on

% Temporal persistence
nexttile;
hold on
plot(log10(ST.chunk_sec), ST.lag1_embedding_corr, 'o-', 'LineWidth', 1.5);
plot(log10(ST.chunk_sec), ST.label_run_frames, 's-', 'LineWidth', 1.5);
xlabel('log10 scale (s)');
ylabel('Persistence metric');
title('Temporal persistence across scales');
legend({'Lag-1 embedding corr','Run length (anchors)'}, 'Location','best');
box off; grid on

% Predictive value and redundancy
nexttile;
hold on
leg = strings(0,1);
if ismember('predict_short_r2', ST.Properties.VariableNames)
    plot(log10(ST.chunk_sec), ST.predict_short_r2, 'o-', 'LineWidth', 1.5);
    leg(end+1) = "Predict short"; %#ok<AGROW>
end
if ismember('predict_long_r2', ST.Properties.VariableNames)
    plot(log10(ST.chunk_sec), ST.predict_long_r2, 's-', 'LineWidth', 1.5);
    leg(end+1) = "Predict long"; %#ok<AGROW>
end
if ismember('redundancy_short', ST.Properties.VariableNames)
    plot(log10(ST.chunk_sec), ST.redundancy_short, '^-', 'LineWidth', 1.2);
    leg(end+1) = "Redundant short"; %#ok<AGROW>
end
if ismember('redundancy_long', ST.Properties.VariableNames)
    plot(log10(ST.chunk_sec), ST.redundancy_long, 'd-', 'LineWidth', 1.2);
    leg(end+1) = "Redundant long"; %#ok<AGROW>
end
xlabel('log10 scale (s)');
ylabel('Cross-scale relation');
title('Predictive value and redundancy');
if ~isempty(leg)
    legend(cellstr(leg), 'Location','best');
end
box off; grid on

% Composite scores
nexttile;
hold on
plot(log10(ST.chunk_sec), ST.composite_micro, 'o-', 'LineWidth', 1.5);
plot(log10(ST.chunk_sec), ST.composite_motif, 's-', 'LineWidth', 1.5);
plot(log10(ST.chunk_sec), ST.composite_context, 'd-', 'LineWidth', 1.5);
if ismember('composite_global', ST.Properties.VariableNames)
    plot(log10(ST.chunk_sec), ST.composite_global, 'k-', 'LineWidth', 2);
end
xlabel('log10 scale (s)');
ylabel('Composite score');
title('Role-specific and global composite scores');
legend({'Micro','Motif','Context','Global'}, 'Location','best');
box off; grid on

% Text panel
nexttile;
axis off
nMicro = sum(Sel.hierarchical_role == "micro");
nMotif = sum(Sel.hierarchical_role == "motif");
nContext = sum(Sel.hierarchical_role == "context");
txt = {
    sprintf('Total scales: %d', height(ST))
    sprintf('Selected scales: %d', height(Sel))
    sprintf('Micro selected: %d', nMicro)
    sprintf('Motif selected: %d', nMotif)
    sprintf('Context selected: %d', nContext)
    };
text(0.02, 0.95, txt, 'VerticalAlignment','top', 'FontSize', 13);

sgtitle('Scale-usefulness scoring overview', 'FontWeight','bold', 'FontSize', 18);

%% Figure 2: z-scored metric heatmap
Fig.heatmap = figure('Color','w', 'Name','Scale metric heatmap', ...
    'Position',[80 80 1500 700]);

metricVars = ["z_pc1","z_cum5","z_sil","z_between","z_persist","z_run","z_pred_short","z_pred_long"];
metricVars = metricVars(ismember(metricVars, string(ST.Properties.VariableNames)));

if ~isempty(metricVars)
    M = zeros(numel(metricVars), height(ST));
    for i = 1:numel(metricVars)
        M(i,:) = ST.(metricVars(i))';
    end
    imagesc(M);
    colormap(parula);
    colorbar
    yticks(1:numel(metricVars));
    yticklabels(cellstr(metricVars));
    xticks(1:height(ST));
    xticklabels(compose('%.2f', ST.chunk_sec));
    xtickangle(90)
    xlabel('Scale index');
    ylabel('Metric');
    title('Z-scored scale metrics');
else
    axis off
    text(0.1, 0.5, 'No z-scored metrics found.', 'FontSize', 12);
end

%% Figure 3: selected scales
Fig.selected = figure('Color','w', 'Name','Selected operational scales', ...
    'Position',[100 100 1500 700]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
hold on
if ismember('composite_global', ST.Properties.VariableNames)
    plot(log10(ST.chunk_sec), ST.composite_global, 'k.-', 'LineWidth', 1.5, 'MarkerSize', 14);
end

roleColors = containers.Map( ...
    {'micro','motif','context'}, ...
    {[0.85 0.33 0.10],[0.30 0.60 0.20],[0.50 0.20 0.70]});

for i = 1:height(Sel)
    role = char(Sel.hierarchical_role(i));
    c = roleColors(role);
    x = log10(Sel.chunk_sec(i));
    y = Sel.composite_global(i);
    scatter(x, y, 80, 'filled', 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k');
    text(x, y, sprintf('  %.2fs', Sel.chunk_sec(i)), 'Color', c, 'FontSize', 11);
end
xlabel('log10 scale (s)');
ylabel('Global composite');
title('Selected operational scales');
box off; grid on

nexttile;
G = groupsummary(Sel, 'hierarchical_role', 'mean', 'chunk_sec');
bar(categorical(cellstr(string(G.hierarchical_role))), G.mean_chunk_sec);
ylabel('Mean selected scale (s)');
title('Hierarchical temporal summary');
box off

%% Figure 4: exemplar trajectories using raw chunk features
Fig.exemplars = figure('Color','w', 'Name','Selected-scale exemplar trajectories', ...
    'Position',[120 120 1900 1200]);

nSel = height(Sel);
nFeat = numel(featNames);
tiledlayout(nSel, nFeat, 'TileSpacing','compact', 'Padding','compact');

for r = 1:nSel
    scaleVal = Sel.chunk_sec(r);
    [~, s] = min(abs([ChunkSet.scale.chunkSec] - scaleVal));

    meta = ChunkSet.scale(s).meta;
    idxSess = find(meta.session_index == P.sessionIndex);

    if isempty(idxSess)
        useIdx = [];
    else
        useIdx = idxSess(randperm(numel(idxSess), min(P.nExampleChunks, numel(idxSess))));
    end

    for c = 1:nFeat
        nexttile;
        hold on

        feat = featNames(c);
        allY = [];
        for j = 1:numel(useIdx)
            x = extract_chunk_feature_trace(ChunkSet, s, useIdx(j), feat, ...
                'invertLog1p', P.invertLog1p);
            tt = (0:numel(x)-1)' ./ ChunkSet.sessions{P.sessionIndex}.fps;
            plot(tt, x, 'LineWidth', 1.0);
            allY = [allY; x(:)]; %#ok<AGROW>
        end

        if r == 1
            title(strrep(feat, '_', '\_'), 'Interpreter','tex');
        end
        if c == 1
            ylabel(sprintf('%s\n%.2fs', char(Sel.hierarchical_role(r)), scaleVal));
        end
        if r == nSel
            xlabel('Time within chunk (s)');
        end

        if ~isempty(allY)
            ylo = prctile(allY, 1);
            yhi = prctile(allY, 99);
            if isfinite(ylo) && isfinite(yhi) && yhi > ylo
                ylim([ylo yhi]);
            end
        end

        box off
        grid on
    end
end

sgtitle('Selected-scale exemplar trajectories', 'FontWeight','bold', 'FontSize', 18);
end
