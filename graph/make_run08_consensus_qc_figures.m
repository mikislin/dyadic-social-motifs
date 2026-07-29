function Manifest = make_run08_consensus_qc_figures(outRoot, params)
%MAKE_RUN08_CONSENSUS_QC_FIGURES Reviewer-facing run-08 consensus audits.

outRoot = string(outRoot);
figDir = fullfile(outRoot, 'figures');
if ~isfolder(figDir)
    mkdir(figDir);
end
style = local_style(params);
Manifest = table('Size', [0 5], ...
    'VariableTypes', {'string','string','string','string','string'}, ...
    'VariableNames', {'figure_id','figure_file','source_csv','description','file_type'});

thresholdPath = fullfile(outRoot, 'graph_consensus_threshold_sensitivity.csv');
edgePath = fullfile(outRoot, 'graph_resampled_edge_support_audit.csv');
nodePath = fullfile(outRoot, 'graph_consensus_node_stability_audit.csv');
replicatePath = fullfile(outRoot, 'graph_replicate_manifest.csv');

files = local_plot_threshold_gates(thresholdPath, figDir, style, params);
Manifest = local_add(Manifest, files, thresholdPath, ...
    "predefined consensus-support gates across k and recurrence thresholds");
files = local_plot_edge_support(edgePath, figDir, style, params);
Manifest = local_add(Manifest, files, edgePath, ...
    "co-inclusion-normalized sparse edge recurrence at primary k");
files = local_plot_node_stability(nodePath, figDir, style, params);
Manifest = local_add(Manifest, files, nodePath, ...
    "selected consensus degree and stable-node coverage by scale and anchor stage");
files = local_plot_replicates(replicatePath, figDir, style);
Manifest = local_add(Manifest, files, replicatePath, ...
    "condition-blind dimension and stage-balanced replicate diagnostics");
end

function style = local_style(params)
style.exportPng = logical(params.figure_export_png);
style.exportPdf = logical(params.figure_export_pdf);
style.dpi = params.figure_dpi;
style.fontName = char(params.figure_font_name);
style.fontSize = params.figure_font_size;
style.blue = [0.10 0.36 0.64];
style.orange = [0.90 0.47 0.13];
style.green = [0.18 0.58 0.39];
style.red = [0.76 0.18 0.20];
style.gray = [0.45 0.48 0.52];
style.grid = [0.85 0.87 0.90];
end

function files = local_plot_threshold_gates(pathText, figDir, style, params)
T = local_read(pathText);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [40 80 1750 700]);
tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
kValues = unique(T.k_neighbors, 'stable')';
colors = lines(numel(kValues));

nexttile;
hold on;
for i = 1:numel(kValues)
    idx = T.k_neighbors == kValues(i);
    plot(T.support_threshold(idx), T.stable_node_fraction(idx), '-o', ...
        'LineWidth', 1.6, 'Color', colors(i, :), 'DisplayName', "k=" + string(kValues(i)));
end
yline(params.consensus_gate_min_stable_node_fraction, '--', 'predefined gate', ...
    'Color', style.gray, 'HandleVisibility', 'off');
xlabel('Conditional support threshold'); ylabel('Stable-node fraction');
title('Node coverage'); legend('Location', 'southwest'); local_format(gca, style);

nexttile;
hold on;
for i = 1:numel(kValues)
    idx = T.k_neighbors == kValues(i);
    plot(T.support_threshold(idx), T.largest_component_fraction(idx), '-o', ...
        'LineWidth', 1.6, 'Color', colors(i, :), 'DisplayName', "k=" + string(kValues(i)));
end
yline(params.consensus_gate_min_largest_component_fraction, '--', 'predefined gate', ...
    'Color', style.gray, 'HandleVisibility', 'off');
xlabel('Conditional support threshold'); ylabel('Largest-component fraction');
title('Consensus connectivity'); local_format(gca, style);

nexttile;
hold on;
for i = 1:numel(kValues)
    idx = T.k_neighbors == kValues(i);
    plot(T.support_threshold(idx), T.rare_stage_stable_fraction(idx), '-o', ...
        'LineWidth', 1.6, 'Color', colors(i, :), 'DisplayName', "k=" + string(kValues(i)));
end
yline(params.consensus_gate_min_rare_stable_fraction, '--', 'predefined gate', ...
    'Color', style.gray, 'HandleVisibility', 'off');
selected = logical(T.selected_for_frozen_handoff);
if any(selected)
    scatter(T.support_threshold(selected), T.rare_stage_stable_fraction(selected), ...
        90, style.red, 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'selected');
end
xlabel('Conditional support threshold'); ylabel('Rare-stage degree coverage');
title('Enrichment-stage coverage'); local_format(gca, style);
sgtitle('Run 08 consensus threshold sensitivity (no motif or cluster definitions)', ...
    'FontName', style.fontName, 'FontWeight', 'bold');
files = local_export(fig, figDir, 'graph_consensus_threshold_sensitivity_audit', style);
close(fig);
end

function files = local_plot_edge_support(pathText, figDir, style, params)
T = local_read(pathText);
varName = sprintf('conditional_neighbor_support_k%d', params.k_neighbors);
support = double(T.(varName));
eligible = logical(T.eligible_for_support);
support = support(eligible);
co = double(T.co_inclusion_replicates(eligible));
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1450 680]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
histogram(support, 40, 'Normalization', 'probability', ...
    'FaceColor', style.blue, 'EdgeColor', 'none');
for x = params.consensus_support_thresholds
    xline(x, '--', compose('%.1f', x), 'Color', style.red);
end
xlabel('Conditional any-direction recurrence'); ylabel('Fraction of union edges');
title(sprintf('Primary k=%d edge support', params.k_neighbors)); local_format(gca, style);
nexttile;
boxchart(categorical(co), support, 'BoxFaceColor', style.orange, 'MarkerStyle', '.');
xlabel('Endpoint co-inclusion replicates'); ylabel('Conditional support');
title('Support denominator audit'); local_format(gca, style);
sgtitle('Sparse consensus edge support (numeric-PC graphs only)', ...
    'FontName', style.fontName, 'FontWeight', 'bold');
files = local_export(fig, figDir, 'graph_consensus_edge_support_audit', style);
close(fig);
end

function files = local_plot_node_stability(pathText, figDir, style, params)
T = local_read(pathText);
scales = unique(double(T.scale_index), 'stable')';
stableFraction = zeros(numel(scales), 1);
medianDegree = zeros(numel(scales), 1);
chunkSec = zeros(numel(scales), 1);
for i = 1:numel(scales)
    idx = double(T.scale_index) == scales(i);
    stableFraction(i) = mean(logical(T.consensus_stable_for_run09(idx)));
    medianDegree(i) = median(double(T.consensus_degree(idx)), 'omitnan');
    chunkSec(i) = median(double(T.chunk_sec(idx)), 'omitnan');
end
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [40 80 1650 680]);
tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
bar(1:numel(scales), stableFraction, 0.75, 'FaceColor', style.green, 'EdgeColor', 'none');
yline(params.consensus_gate_min_stable_node_fraction, '--', 'all-node gate', 'Color', style.gray);
xticks(1:numel(scales)); xticklabels(compose('%.3g', chunkSec)); xtickangle(45);
xlabel('Chunk duration (s)'); ylabel('Stable-node fraction'); title('Coverage by scale'); local_format(gca, style);
nexttile;
bar(1:numel(scales), medianDegree, 0.75, 'FaceColor', style.blue, 'EdgeColor', 'none');
yline(params.consensus_stable_min_degree, '--', 'stable degree', 'Color', style.gray);
xticks(1:numel(scales)); xticklabels(compose('%.3g', chunkSec)); xtickangle(45);
xlabel('Chunk duration (s)'); ylabel('Median consensus degree'); title('Degree by scale'); local_format(gca, style);
nexttile;
stages = unique(string(T.anchor_stage), 'stable');
values = zeros(numel(stages), 1);
for i = 1:numel(stages)
    idx = string(T.anchor_stage) == stages(i);
    values(i) = mean(logical(T.consensus_stable_for_run09(idx)));
end
bar(categorical(replace(stages, "_", " ")), values, 0.65, ...
    'FaceColor', style.orange, 'EdgeColor', 'none');
ylabel('Stable-node fraction'); title('Coverage by sampling stage'); local_format(gca, style);
sgtitle('Selected run 08 consensus graph stability', ...
    'FontName', style.fontName, 'FontWeight', 'bold');
files = local_export(fig, figDir, 'graph_consensus_node_stability_audit', style);
close(fig);
end

function files = local_plot_replicates(pathText, figDir, style)
T = local_read(pathText);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [60 80 1500 680]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
bar(double(T.replicate_id), double(T.mean_same_scale_candidate_fraction), ...
    0.75, 'FaceColor', style.blue, 'EdgeColor', 'none');
xlabel('Replicate'); ylabel('Same-scale candidate fraction');
title('Replicate scale composition'); local_format(gca, style);
nexttile;
hold on;
designs = unique(string(T.dimension_design), 'stable');
colors = lines(numel(designs));
for i = 1:numel(designs)
    idx = string(T.dimension_design) == designs(i);
    scatter(double(T.n_graph_dimensions(idx)), double(T.elapsed_seconds(idx)), ...
        70, colors(i, :), 'filled', 'DisplayName', replace(designs(i), "_", " "));
end
xlabel('Numeric graph dimensions'); ylabel('Elapsed seconds');
title('Exact graph replicate cost'); legend('Location', 'best'); local_format(gca, style);
sgtitle('Dimension/stage-balanced graph ensemble audit', ...
    'FontName', style.fontName, 'FontWeight', 'bold');
files = local_export(fig, figDir, 'graph_consensus_replicate_audit', style);
close(fig);
end

function T = local_add(T, files, sourceCsv, description)
for i = 1:numel(files)
    [~, stem, ext] = fileparts(files(i));
    one = table(string(stem), string(files(i)), string(sourceCsv), ...
        string(description), erase(string(ext), "."), ...
        'VariableNames', {'figure_id','figure_file','source_csv','description','file_type'});
    T = [T; one]; %#ok<AGROW>
end
end

function files = local_export(fig, figDir, stem, style)
files = strings(0, 1);
if style.exportPng
    pathText = fullfile(figDir, stem + ".png");
    exportgraphics(fig, pathText, 'Resolution', style.dpi);
    files(end + 1, 1) = string(pathText);
end
if style.exportPdf
    pathText = fullfile(figDir, stem + ".pdf");
    exportgraphics(fig, pathText, 'ContentType', 'image', 'Resolution', style.dpi);
    files(end + 1, 1) = string(pathText);
end
end

function local_format(ax, style)
set(ax, 'Box', 'off', 'FontName', style.fontName, 'FontSize', style.fontSize, ...
    'LineWidth', 0.8, 'TickDir', 'out');
grid(ax, 'on');
ax.GridColor = style.grid;
ax.GridAlpha = 0.35;
end

function T = local_read(pathText)
opts = detectImportOptions(pathText, 'FileType', 'text', ...
    'Delimiter', ',', 'TextType', 'string');
T = readtable(pathText, opts);
end
