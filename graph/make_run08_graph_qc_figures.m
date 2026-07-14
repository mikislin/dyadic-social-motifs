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
    exportgraphics(fig, pdfPath, 'ContentType', 'vector');
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
