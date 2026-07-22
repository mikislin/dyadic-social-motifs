function figureManifest = make_run08_graph_qc_figures(outRoot, params)
%MAKE_RUN08_GRAPH_QC_FIGURES Generate CSV-backed run-08 graph QC figures.

outRoot = string(outRoot);
figDir = fullfile(outRoot, 'figures');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

style = local_style(params);
rows = table();

rows = local_add_figure(rows, local_plot_topology_overview(outRoot, style), ...
    string(fullfile(outRoot, 'graph_node_manifest.csv')) + ";" + ...
    string(fullfile(outRoot, 'graph_degree_audit.csv')), ...
    "condition-blind graph overview on run_07 global PC audit axes");
rows = local_add_figure(rows, local_plot_degree_radius(outRoot, style), ...
    string(fullfile(outRoot, 'graph_degree_audit.csv')) + ";" + ...
    string(fullfile(outRoot, 'graph_neighbor_composition_audit.csv')), ...
    "graph degree and nearest-neighbor radius diagnostics by primary scale");
rows = local_add_figure(rows, local_plot_components(outRoot, style), ...
    fullfile(outRoot, 'graph_component_audit.csv'), ...
    "connected-component size and within-component support audit");
rows = local_add_figure(rows, local_plot_neighbor_composition(outRoot, style), ...
    fullfile(outRoot, 'graph_neighbor_composition_audit.csv'), ...
    "post-fit neighbor composition audit by scale");
rows = local_add_figure(rows, local_plot_k_sensitivity(outRoot, style), ...
    fullfile(outRoot, 'graph_k_sensitivity_audit.csv'), ...
    "condition-blind graph topology sensitivity across k values");
rows = local_add_figure(rows, local_plot_event_coverage(outRoot, style), ...
    fullfile(outRoot, 'graph_rare_event_coverage_audit.csv'), ...
    "post-fit rare and social event coverage audit");
rows = local_add_figure(rows, local_plot_scale_mixing(outRoot, style), ...
    fullfile(outRoot, 'graph_scale_mixing_matrix_audit.csv'), ...
    "complete directed scale-mixing matrix with random-mixing normalization");
rows = local_add_figure(rows, local_plot_global_pca_density(outRoot, style), ...
    string(fullfile(outRoot, 'graph_node_manifest.csv')) + ";" + ...
    string(fullfile(outRoot, 'graph_global_pca_cumulative_variance_audit.csv')), ...
    "global PCA density, scale-dependent dispersion, and cumulative variance audit");
rows = local_add_figure(rows, local_plot_global_pca_3d(outRoot, style), ...
    fullfile(outRoot, 'graph_node_manifest.csv'), ...
    "three-dimensional global PCA visualization audit; not a motif definition");
rows = local_add_figure(rows, local_plot_event_prevalence_fold(outRoot, style), ...
    fullfile(outRoot, 'graph_event_prevalence_fold_audit.csv'), ...
    "baseline-versus-enriched event prevalence and fold-enrichment audit");
rows = local_add_figure(rows, local_plot_umap_2d(outRoot, style), ...
    fullfile(outRoot, 'graph_umap_embedding_audit.csv'), ...
    "deterministic condition-blind UMAP 2D visualization audit only");
rows = local_add_figure(rows, local_plot_umap_3d(outRoot, style), ...
    fullfile(outRoot, 'graph_umap_embedding_audit.csv'), ...
    "deterministic condition-blind UMAP 3D visualization audit only");
rows = local_add_figure(rows, local_plot_sensitivity_overview(outRoot, style), ...
    string(fullfile(outRoot, 'graph_anchor_stage_sensitivity_audit.csv')) + ";" + ...
    string(fullfile(outRoot, 'graph_session_excluded_sensitivity_audit.csv')) + ";" + ...
    string(fullfile(outRoot, 'graph_neighborhood_resampling_audit.csv')), ...
    "anchor-stage, session-exclusion, and neighborhood-resampling sensitivity audit");

figureManifest = rows;
end

function style = local_style(params)
style = struct();
style.fontName = char(params.figure_font_name);
style.fontSize = params.figure_font_size;
style.titleFontSize = params.figure_font_size + 2;
style.dpi = params.figure_dpi;
style.exportPng = logical(params.figure_export_png);
style.exportPdf = logical(params.figure_export_pdf);
style.blue = [0.000 0.447 0.698];
style.gold = [0.902 0.624 0.000];
style.green = [0.000 0.620 0.451];
style.red = [0.800 0.200 0.200];
style.purple = [0.494 0.184 0.556];
style.gray = [0.45 0.45 0.45];
style.lightGray = [0.86 0.86 0.86];
style.grid = [0.86 0.86 0.86];
end

function files = local_plot_topology_overview(outRoot, style)
N = local_read_csv(fullfile(outRoot, 'graph_node_manifest.csv'));
D = local_read_csv(fullfile(outRoot, 'graph_degree_audit.csv'));
if ~all(ismember(["graph_plot_pc1", "graph_plot_pc2"], string(N.Properties.VariableNames)))
    files = strings(0, 1);
    return
end

[tf, loc] = ismember(N.graph_node_id, D.graph_node_id);
assert(all(tf), 'make_run08_graph_qc_figures:NodeDegreeMismatch', ...
    'graph_node_manifest rows must match graph_degree_audit rows by graph_node_id.');
D = D(loc, :);

rng(108);
n = height(N);
if n > 9000
    idx = sort(randperm(n, 9000));
else
    idx = 1:n;
end

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1500 720]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
scatter(N.graph_plot_pc1(idx), N.graph_plot_pc2(idx), 10, double(N.chunk_sec(idx)), 'filled');
xlabel('Global PC1', 'Interpreter', 'none');
ylabel('Global PC2', 'Interpreter', 'none');
title('Graph input space colored by scale', 'Interpreter', 'none');
colorbar;
local_format_axis(gca, style);

nexttile;
scatter(N.graph_plot_pc1(idx), N.graph_plot_pc2(idx), 10, D.undirected_degree(idx), 'filled');
xlabel('Global PC1', 'Interpreter', 'none');
ylabel('Global PC2', 'Interpreter', 'none');
title('Graph input space colored by degree', 'Interpreter', 'none');
colorbar;
local_format_axis(gca, style);

files = local_export(fig, fullfile(outRoot, 'figures'), ...
    'graph_embedding_topology_overview_audit', style);
close(fig);
end

function files = local_plot_degree_radius(outRoot, style)
D = local_read_csv(fullfile(outRoot, 'graph_degree_audit.csv'));
S = local_read_csv(fullfile(outRoot, 'graph_neighbor_composition_audit.csv'));

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1500 720]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
histogram(D.undirected_degree, 'BinMethod', 'integers', ...
    'FaceColor', style.blue, 'EdgeColor', 'none');
xlabel('Undirected degree', 'Interpreter', 'none');
ylabel('Node count', 'Interpreter', 'none');
title('Graph degree distribution', 'Interpreter', 'none');
local_format_axis(gca, style);

nexttile;
hold on;
plot(S.chunk_sec, S.mean_neighbor_distance, 'o-', ...
    'Color', style.green, 'MarkerFaceColor', style.green, 'LineWidth', 1.5);
plot(S.chunk_sec, local_scale_median(D, 'knn_radius'), 's-', ...
    'Color', style.gold, 'MarkerFaceColor', style.gold, 'LineWidth', 1.5);
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Distance', 'Interpreter', 'none');
title('Neighbor distance and kNN radius by scale', 'Interpreter', 'none');
legend({'Mean neighbor distance','Median kNN radius'}, ...
    'Location', 'best', 'Box', 'off');
local_format_axis(gca, style);

files = local_export(fig, fullfile(outRoot, 'figures'), ...
    'graph_degree_radius_audit', style);
close(fig);
end

function files = local_plot_components(outRoot, style)
C = local_read_csv(fullfile(outRoot, 'graph_component_audit.csv'));
if isempty(C)
    files = strings(0, 1);
    return
end
nShow = min(height(C), 25);
Cshow = C(1:nShow, :);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1450 720]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
bar(Cshow.node_fraction, 'FaceColor', style.purple, 'EdgeColor', 'none');
xlabel('Component rank by size', 'Interpreter', 'none');
ylabel('Node fraction', 'Interpreter', 'none');
title('Connected-component size audit', 'Interpreter', 'none');
local_format_axis(gca, style);

nexttile;
hold on;
plot(Cshow.n_scales, 'o-', 'Color', style.green, ...
    'MarkerFaceColor', style.green, 'LineWidth', 1.5);
plot(Cshow.n_sessions, 's-', 'Color', style.gold, ...
    'MarkerFaceColor', style.gold, 'LineWidth', 1.5);
xlabel('Component rank by size', 'Interpreter', 'none');
ylabel('Support count', 'Interpreter', 'none');
title('Scale and session support in components', 'Interpreter', 'none');
legend({'Scales','Sessions'}, 'Location', 'best', 'Box', 'off');
local_format_axis(gca, style);

files = local_export(fig, fullfile(outRoot, 'figures'), ...
    'graph_component_summary_audit', style);
close(fig);
end

function files = local_plot_neighbor_composition(outRoot, style)
T = local_read_csv(fullfile(outRoot, 'graph_neighbor_composition_audit.csv'));

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1450 720]);
hold on;
plot(T.chunk_sec, T.mean_same_scale_neighbor_fraction, 'o-', ...
    'Color', style.blue, 'MarkerFaceColor', style.blue, 'LineWidth', 1.5);
plot(T.chunk_sec, T.mean_same_session_neighbor_fraction, 's-', ...
    'Color', style.green, 'MarkerFaceColor', style.green, 'LineWidth', 1.5);
plot(T.chunk_sec, T.mean_same_arena_neighbor_fraction_posthoc, 'd-', ...
    'Color', style.gold, 'MarkerFaceColor', style.gold, 'LineWidth', 1.5);
set(gca, 'XScale', 'log');
ylim([0 1.02]);
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Neighbor fraction', 'Interpreter', 'none');
title('Neighbor composition after graph fitting', 'Interpreter', 'none');
legend({'Same scale','Same session','Same arena (post-fit)'}, ...
    'Location', 'best', 'Box', 'off');
local_format_axis(gca, style);

files = local_export(fig, fullfile(outRoot, 'figures'), ...
    'graph_neighbor_composition_audit', style);
close(fig);
end

function files = local_plot_k_sensitivity(outRoot, style)
T = local_read_csv(fullfile(outRoot, 'graph_k_sensitivity_audit.csv'));

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1450 720]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(T.k_neighbors, T.edge_jaccard_to_primary_k, 'o-', ...
    'Color', style.green, 'MarkerFaceColor', style.green, 'LineWidth', 1.6);
ylim([0 1.02]);
xlabel('k neighbors', 'Interpreter', 'none');
ylabel('Edge Jaccard to primary k', 'Interpreter', 'none');
title('Edge-set sensitivity', 'Interpreter', 'none');
local_format_axis(gca, style);

nexttile;
yyaxis left;
plot(T.k_neighbors, T.n_components, 's-', ...
    'Color', style.blue, 'MarkerFaceColor', style.blue, 'LineWidth', 1.6);
ylabel('Components', 'Interpreter', 'none');
yyaxis right;
plot(T.k_neighbors, T.largest_component_fraction, 'd-', ...
    'Color', style.gold, 'MarkerFaceColor', style.gold, 'LineWidth', 1.6);
ylabel('Largest component fraction', 'Interpreter', 'none');
xlabel('k neighbors', 'Interpreter', 'none');
title('Connectivity sensitivity', 'Interpreter', 'none');
local_format_axis(gca, style);

files = local_export(fig, fullfile(outRoot, 'figures'), ...
    'graph_k_sensitivity_audit', style);
close(fig);
end

function files = local_plot_event_coverage(outRoot, style)
T = local_read_csv(fullfile(outRoot, 'graph_rare_event_coverage_audit.csv'));
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1550 780]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

if isempty(T)
    nexttile([1 2]);
    axis off;
    text(0.02, 0.55, 'Rare/social event coverage audit was unavailable.', ...
        'FontName', style.fontName, 'FontSize', style.titleFontSize, ...
        'Interpreter', 'none');
    files = local_export(fig, fullfile(outRoot, 'figures'), ...
        'graph_rare_event_coverage_audit', style);
    close(fig);
    return
end

events = unique(string(T.event_id), 'stable');
scales = unique(double(T.chunk_sec), 'stable')';
M = NaN(numel(events), numel(scales));
A = NaN(numel(events), numel(scales));
for i = 1:numel(events)
    for j = 1:numel(scales)
        idx = string(T.event_id) == events(i) & double(T.chunk_sec) == scales(j);
        if any(idx)
            M(i, j) = mean(T.selected_event_coverage_fraction(idx), 'omitnan');
            A(i, j) = mean(T.available_event_fraction(idx), 'omitnan');
        end
    end
end

nexttile;
imagesc(M);
set(gca, 'YTick', 1:numel(events), 'YTickLabel', local_clean_label(events), ...
    'XTick', 1:numel(scales), 'XTickLabel', compose('%.3g', scales), ...
    'TickLabelInterpreter', 'none');
xtickangle(45);
caxis([0 1]);
colorbar;
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Event audit', 'Interpreter', 'none');
title('Graph coverage of available events', 'Interpreter', 'none');
local_format_axis(gca, style);

nexttile;
imagesc(A);
set(gca, 'YTick', 1:numel(events), 'YTickLabel', local_clean_label(events), ...
    'XTick', 1:numel(scales), 'XTickLabel', compose('%.3g', scales), ...
    'TickLabelInterpreter', 'none');
xtickangle(45);
caxis([0 1]);
colorbar;
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Event audit', 'Interpreter', 'none');
title('Event prevalence in primary bank', 'Interpreter', 'none');
local_format_axis(gca, style);

files = local_export(fig, fullfile(outRoot, 'figures'), ...
    'graph_rare_event_coverage_audit', style);
close(fig);
end

function files = local_plot_scale_mixing(outRoot, style)
pathText = fullfile(outRoot, 'graph_scale_mixing_matrix_audit.csv');
if ~isfile(pathText)
    files = strings(0, 1);
    return
end
T = local_read_csv(pathText);
scales = unique(double(T.source_chunk_sec), 'stable');
M = nan(numel(scales));
for i = 1:numel(scales)
    for j = 1:numel(scales)
        idx = double(T.source_chunk_sec) == scales(i) & double(T.target_chunk_sec) == scales(j);
        if any(idx)
            M(i, j) = log2(max(double(T.observed_over_random_ratio(find(idx, 1))), eps));
        end
    end
end
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 980 820]);
imagesc(M);
axis square;
set(gca, 'YTick', 1:numel(scales), 'YTickLabel', compose('%.3g', scales), ...
    'XTick', 1:numel(scales), 'XTickLabel', compose('%.3g', scales));
xtickangle(45);
xlabel('Target chunk scale (s)', 'Interpreter', 'none');
ylabel('Source chunk scale (s)', 'Interpreter', 'none');
title('Scale mixing: log_2(observed / random expectation)', 'Interpreter', 'tex');
cb = colorbar;
cb.Label.String = 'log_2 enrichment';
local_format_axis(gca, style);
files = local_export(fig, fullfile(outRoot, 'figures'), ...
    'graph_scale_mixing_matrix_audit', style);
close(fig);
end

function files = local_plot_global_pca_density(outRoot, style)
nodePath = fullfile(outRoot, 'graph_node_manifest.csv');
variancePath = fullfile(outRoot, 'graph_global_pca_cumulative_variance_audit.csv');
if ~(isfile(nodePath) && isfile(variancePath))
    files = strings(0, 1);
    return
end
N = local_read_csv(nodePath);
V = local_read_csv(variancePath);
if ~all(ismember({'graph_plot_pc1','graph_plot_pc2'}, N.Properties.VariableNames))
    files = strings(0, 1);
    return
end
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [40 80 1800 620]);
tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
histogram2(double(N.graph_plot_pc1), double(N.graph_plot_pc2), [110 110], ...
    'DisplayStyle', 'tile', 'ShowEmptyBins', 'off');
view(2);
set(gca, 'ColorScale', 'log');
colorbar;
xlabel('Global PC1'); ylabel('Global PC2');
title('Point density (all graph nodes)');
local_format_axis(gca, style);

nexttile;
scales = unique(double(N.chunk_sec), 'stable');
medRadius = nan(size(scales));
p90Radius = nan(size(scales));
r = hypot(double(N.graph_plot_pc1), double(N.graph_plot_pc2));
for i = 1:numel(scales)
    idx = double(N.chunk_sec) == scales(i);
    medRadius(i) = median(r(idx), 'omitnan');
    p90Radius(i) = quantile(r(idx), 0.90);
end
semilogx(scales, medRadius, '-o', 'Color', style.blue, 'LineWidth', 1.5, ...
    'MarkerFaceColor', style.blue); hold on;
semilogx(scales, p90Radius, '-s', 'Color', style.gold, 'LineWidth', 1.5, ...
    'MarkerFaceColor', style.gold);
xlabel('Chunk scale (s)'); ylabel('Radius in displayed PC1-PC2 plane');
title('Scale-dependent projected dispersion');
legend({'Median radius','90th percentile'}, 'Location', 'best', 'Box', 'off');
local_format_axis(gca, style);

nexttile;
plot(double(V.global_pc_index), double(V.cumulative_explained_variance_percent), ...
    '-', 'Color', style.purple, 'LineWidth', 1.8); hold on;
xline(2, ':', 'PC2');
lastGraph = max(double(V.global_pc_index(logical(V.selected_for_run08_graph))));
xline(lastGraph, '--', sprintf('Graph PC%d', lastGraph));
yline(100, ':');
xlabel('Number of global PCs'); ylabel('Cumulative explained variance (%)');
title('Displayed axes versus graph input');
xlim([1 max(double(V.global_pc_index))]); ylim([0 101]);
local_format_axis(gca, style);

files = local_export(fig, fullfile(outRoot, 'figures'), ...
    'graph_global_pca_density_variance_audit', style);
close(fig);
end

function files = local_plot_global_pca_3d(outRoot, style)
pathText = fullfile(outRoot, 'graph_node_manifest.csv');
if ~isfile(pathText)
    files = strings(0, 1);
    return
end
N = local_read_csv(pathText);
if ~all(ismember({'graph_plot_pc1','graph_plot_pc2','graph_plot_pc3'}, N.Properties.VariableNames))
    files = strings(0, 1);
    return
end
idx = local_plot_sample(height(N), 12000, 108);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1100 860]);
scatter3(double(N.graph_plot_pc1(idx)), double(N.graph_plot_pc2(idx)), ...
    double(N.graph_plot_pc3(idx)), 8, log10(double(N.chunk_sec(idx))), 'filled');
xlabel('Global PC1'); ylabel('Global PC2'); zlabel('Global PC3');
title('Global PCA 3D audit (visualization only)');
cb = colorbar; cb.Label.String = 'log_{10} chunk scale (s)';
view(42, 24);
axis vis3d;
local_format_axis(gca, style);
files = local_export(fig, fullfile(outRoot, 'figures'), ...
    'graph_global_pca_3d_audit', style);
close(fig);
end

function files = local_plot_event_prevalence_fold(outRoot, style)
pathText = fullfile(outRoot, 'graph_event_prevalence_fold_audit.csv');
if ~isfile(pathText)
    files = strings(0, 1);
    return
end
T = local_read_csv(pathText);
T = T(string(T.aggregation_scope) == "all_scales", :);
if isempty(T)
    files = strings(0, 1);
    return
end
events = local_clean_label(string(T.event_id));
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [40 80 1700 760]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
barh(1:height(T), [double(T.baseline_prevalence), double(T.current_prevalence)]);
set(gca, 'XScale', 'log', 'YTick', 1:height(T), 'YTickLabel', events, ...
    'TickLabelInterpreter', 'none', 'YDir', 'reverse');
xlabel('Event prevalence'); ylabel('Event audit');
title('Baseline and rare-enriched bank prevalence');
legend({'Baseline primary','Rare enriched'}, 'Location', 'southoutside', ...
    'Orientation', 'horizontal', 'Box', 'off');
local_format_axis(gca, style);
nexttile;
barh(1:height(T), double(T.log2_prevalence_fold), 'FaceColor', style.green);
xline(0, ':');
set(gca, 'YTick', 1:height(T), 'YTickLabel', events, ...
    'TickLabelInterpreter', 'none', 'YDir', 'reverse');
xlabel('log_2 enriched / baseline prevalence'); ylabel('Event audit');
title('Prevalence fold enrichment (half-count stabilized)');
local_format_axis(gca, style);
files = local_export(fig, fullfile(outRoot, 'figures'), ...
    'graph_event_prevalence_fold_audit', style);
close(fig);
end

function files = local_plot_umap_2d(outRoot, style)
pathText = fullfile(outRoot, 'graph_umap_embedding_audit.csv');
if ~isfile(pathText)
    files = strings(0, 1);
    return
end
U = local_read_csv(pathText);
if isempty(U)
    files = strings(0, 1);
    return
end
idx = local_plot_sample(height(U), 14000, 108);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [40 80 1500 680]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
histogram2(double(U.umap2_x), double(U.umap2_y), [120 120], ...
    'DisplayStyle', 'tile', 'ShowEmptyBins', 'off');
view(2); set(gca, 'ColorScale', 'log'); colorbar;
xlabel('UMAP 1'); ylabel('UMAP 2'); title('UMAP 2D density');
local_format_axis(gca, style);
nexttile;
scatter(double(U.umap2_x(idx)), double(U.umap2_y(idx)), 8, ...
    log10(double(U.chunk_sec(idx))), 'filled');
xlabel('UMAP 1'); ylabel('UMAP 2'); title('UMAP 2D colored by chunk scale');
cb = colorbar; cb.Label.String = 'log_{10} chunk scale (s)';
local_format_axis(gca, style);
sgtitle('Condition-blind UMAP audit: visualization only, not motif evidence', ...
    'FontName', style.fontName, 'FontSize', style.titleFontSize);
files = local_export(fig, fullfile(outRoot, 'figures'), ...
    'graph_umap_2d_visualization_audit', style);
close(fig);
end

function files = local_plot_umap_3d(outRoot, style)
pathText = fullfile(outRoot, 'graph_umap_embedding_audit.csv');
if ~isfile(pathText)
    files = strings(0, 1);
    return
end
U = local_read_csv(pathText);
if isempty(U)
    files = strings(0, 1);
    return
end
idx = local_plot_sample(height(U), 14000, 108);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1100 860]);
scatter3(double(U.umap3_x(idx)), double(U.umap3_y(idx)), double(U.umap3_z(idx)), ...
    8, log10(double(U.chunk_sec(idx))), 'filled');
xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3');
title('Condition-blind UMAP 3D audit (visualization only)');
cb = colorbar; cb.Label.String = 'log_{10} chunk scale (s)';
view(42, 24); axis vis3d;
local_format_axis(gca, style);
files = local_export(fig, fullfile(outRoot, 'figures'), ...
    'graph_umap_3d_visualization_audit', style);
close(fig);
end

function files = local_plot_sensitivity_overview(outRoot, style)
stagePath = fullfile(outRoot, 'graph_anchor_stage_sensitivity_audit.csv');
sessionPath = fullfile(outRoot, 'graph_session_excluded_sensitivity_audit.csv');
resamplePath = fullfile(outRoot, 'graph_neighborhood_resampling_audit.csv');
if ~(isfile(stagePath) && isfile(sessionPath) && isfile(resamplePath))
    files = strings(0, 1);
    return
end
S = local_read_csv(stagePath);
E = local_read_csv(sessionPath);
R = local_read_csv(resamplePath);
S = S(string(S.status) == "completed", :);
E = E(string(E.source_scope) == "source_scale", :);
R = R(string(R.resampling_type) ~= "reference", :);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [30 80 1800 680]);
tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
variantLabels = local_clean_label(string(S.graph_variant));
bar(categorical(variantLabels, variantLabels, variantLabels), ...
    [double(S.mean_same_scale_neighbor_fraction), double(S.mean_same_session_neighbor_fraction)]);
set(gca, 'TickLabelInterpreter', 'none');
ylabel('Neighbor fraction'); title('Anchor-stage graph sensitivities');
legend({'Same scale','Same session'}, 'Location', 'northoutside', ...
    'Orientation', 'horizontal', 'Box', 'off');
local_format_axis(gca, style);
nexttile;
semilogx(double(E.chunk_sec), double(E.mean_primary_neighbor_retention_fraction), ...
    '-o', 'Color', style.blue, 'MarkerFaceColor', style.blue, 'LineWidth', 1.5); hold on;
semilogx(double(E.chunk_sec), double(E.median_neighbor_jaccard_to_primary), ...
    '-s', 'Color', style.gold, 'MarkerFaceColor', style.gold, 'LineWidth', 1.5);
xlabel('Chunk scale (s)'); ylabel('Retention / Jaccard');
title('Session-excluded sensitivity');
legend({'Primary neighbor retention','Median Jaccard'}, 'Location', 'best', 'Box', 'off');
local_format_axis(gca, style);
nexttile;
types = unique(string(R.resampling_type), 'stable');
data = cell(numel(types), 1);
for i = 1:numel(types)
    data{i} = double(R.median_neighbor_jaccard_to_panel_reference(string(R.resampling_type) == types(i)));
end
maxN = max(cellfun(@numel, data));
M = nan(maxN, numel(types));
for i = 1:numel(types), M(1:numel(data{i}), i) = data{i}; end
typeLabels = local_clean_label(types);
typeCategories = categorical(types, types, typeLabels);
boxchart(repelem(typeCategories, maxN), M(:));
set(gca, 'TickLabelInterpreter', 'none');
ylabel('Median node-neighborhood Jaccard');
title('Resampling-based neighborhood retention');
local_format_axis(gca, style);
files = local_export(fig, fullfile(outRoot, 'figures'), ...
    'graph_sensitivity_overview_audit', style);
close(fig);
end

function idx = local_plot_sample(n, maxN, seed)
rng(seed, 'twister');
if n > maxN
    idx = sort(randperm(n, maxN));
else
    idx = 1:n;
end
end

function values = local_scale_median(T, varName)
scales = unique(T.scale_index, 'stable')';
values = NaN(numel(scales), 1);
for i = 1:numel(scales)
    idx = T.scale_index == scales(i);
    values(i) = median(T.(varName)(idx), 'omitnan');
end
end

function rows = local_add_figure(rows, files, sourceCsv, description)
for i = 1:numel(files)
    one = table();
    [~, stem, ext] = fileparts(files(i));
    one.figure_id = string(stem);
    one.figure_file = string(files(i));
    one.source_csv = string(sourceCsv);
    one.description = string(description);
    one.file_type = erase(string(ext), ".");
    rows = [rows; one]; %#ok<AGROW>
end
end

function files = local_export(fig, figDir, stem, style)
files = strings(0, 1);
if style.exportPng
    pngPath = fullfile(figDir, stem + ".png");
    exportgraphics(fig, pngPath, 'Resolution', style.dpi);
    files(end + 1, 1) = string(pngPath);
end
if style.exportPdf
    pdfPath = fullfile(figDir, stem + ".pdf");
    % Dense scatter/density figures can crash MATLAB's Windows vector
    % hardcopy renderer. A high-resolution raster embedded in PDF preserves
    % the reviewer-facing artifact without changing any analytical values.
    exportgraphics(fig, pdfPath, 'ContentType', 'image', 'Resolution', style.dpi);
    files(end + 1, 1) = string(pdfPath);
end
end

function labels = local_clean_label(values)
labels = string(values);
labels = replace(labels, "_", " ");
labels = replace(labels, "|", " | ");
end

function local_format_axis(ax, style)
set(ax, 'Box', 'off', 'FontName', style.fontName, 'FontSize', style.fontSize, ...
    'LineWidth', 0.8, 'TickDir', 'out');
grid(ax, 'on');
ax.GridColor = style.grid;
ax.GridAlpha = 0.35;
end

function T = local_read_csv(pathText)
opts = detectImportOptions(pathText, 'FileType', 'text', ...
    'Delimiter', ',', 'TextType', 'string');
T = readtable(pathText, opts);
end
