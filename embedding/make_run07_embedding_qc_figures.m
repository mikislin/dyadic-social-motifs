function figureManifest = make_run07_embedding_qc_figures(outRoot, params)
%MAKE_RUN07_EMBEDDING_QC_FIGURES Generate CSV-backed run-07 QC figures.

outRoot = string(outRoot);
figDir = fullfile(outRoot, 'figures');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

style = local_style(params);
rows = table();

rows = local_add_figure(rows, local_plot_anchor_coverage(outRoot, style), ...
    fullfile(outRoot, 'embedding_input_anchor_audit.csv'), ...
    "primary scale-specific anchor coverage used by run_07 embedding");
rows = local_add_figure(rows, local_plot_pca_dimensions(outRoot, style), ...
    fullfile(outRoot, 'embedding_pca_by_scale.csv'), ...
    "per-scale PCA dimensionality from run_06 guidance and run_07 caps");
rows = local_add_figure(rows, local_plot_preprocess_dimension_audit(outRoot, style), ...
    fullfile(outRoot, 'embedding_preprocess_dimension_audit.csv'), ...
    "condition-blind sparse-feature scale safeguard and standardized-tail audit");
rows = local_add_figure(rows, local_plot_anchor_stage_sensitivity(outRoot, style), ...
    fullfile(outRoot, 'embedding_anchor_stage_pca_sensitivity_audit.csv'), ...
    "audit-only base rare-enriched and combined PCA subspace sensitivity");
rows = local_add_figure(rows, local_plot_matrix_manifest(outRoot, style), ...
    fullfile(outRoot, 'embedding_matrix_manifest.csv'), ...
    "embedding matrix finite-value and artifact audit");
rows = local_add_figure(rows, local_plot_scale_weights(outRoot, style), ...
    fullfile(outRoot, 'embedding_scale_weight_audit.csv'), ...
    "auditable predefined scale weights before global PCA");
rows = local_add_figure(rows, local_plot_stability(outRoot, style), ...
    fullfile(outRoot, 'embedding_stability_audit.csv'), ...
    "condition-blind split-half PCA stability by primary scale");
rows = local_add_figure(rows, local_plot_arena_sensitivity(outRoot, style), ...
    fullfile(outRoot, 'embedding_arena_sensitivity_audit.csv'), ...
    "post-fit arena sensitivity audit only");
rows = local_add_figure(rows, local_plot_global_embedding(outRoot, style), ...
    string(fullfile(outRoot, 'embedding_global_scores.csv')) + ";" + ...
    string(fullfile(outRoot, 'embedding_row_manifest.csv')), ...
    "global condition-blind embedding overview colored by scale and arena for audit");

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
style.grid = [0.86 0.86 0.86];
end

function files = local_plot_anchor_coverage(outRoot, style)
T = local_read_csv(fullfile(outRoot, 'embedding_input_anchor_audit.csv'));
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1450 720]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
plot(T.chunk_sec, T.n_primary_anchor_rows_available, 'o-', 'Color', style.gray, 'LineWidth', 1.5);
plot(T.chunk_sec, T.n_run07_anchor_rows_used, 's-', 'Color', style.green, 'LineWidth', 1.7);
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Anchor rows', 'Interpreter', 'none');
title('Primary bank anchors used by embedding', 'Interpreter', 'none');
legend({'Primary bank available','Run 07 used'}, 'Location', 'best', 'Box', 'off');
local_format_axis(gca, style);

nexttile;
hold on;
plot(T.chunk_sec, T.n_sessions_available, 'o-', 'Color', style.gray, 'LineWidth', 1.5);
plot(T.chunk_sec, T.n_sessions_used, 's-', 'Color', style.green, 'LineWidth', 1.7);
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Sessions', 'Interpreter', 'none');
title('Session support by primary scale', 'Interpreter', 'none');
legend({'Available','Used'}, 'Location', 'best', 'Box', 'off');
local_format_axis(gca, style);

files = local_export(fig, fullfile(outRoot, 'figures'), 'embedding_anchor_coverage_by_scale', style);
close(fig);
end

function files = local_plot_pca_dimensions(outRoot, style)
T = local_read_csv(fullfile(outRoot, 'embedding_pca_by_scale.csv'));
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1500 760]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
plot(T.chunk_sec, T.run06_recommended_pcs_for_next_embedding, 'o-', 'Color', style.gray, 'LineWidth', 1.5);
plot(T.chunk_sec, T.run07_role_pc_cap, 'd--', 'Color', style.purple, 'LineWidth', 1.3);
plot(T.chunk_sec, T.n_pcs_selected, 's-', 'Color', style.green, 'LineWidth', 1.7);
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('PC count', 'Interpreter', 'none');
title('Scale-specific PCA dimensions', 'Interpreter', 'none');
legend({'Run 06 recommendation','Run 07 cap','Selected PCs'}, 'Location', 'best', 'Box', 'off');
local_format_axis(gca, style);

nexttile;
hold on;
plot(T.chunk_sec, T.pc1_explained, 'o-', 'Color', style.blue, 'LineWidth', 1.5);
plot(T.chunk_sec, T.cum5_explained, 's-', 'Color', style.gold, 'LineWidth', 1.5);
plot(T.chunk_sec, T.cum_selected_explained, 'd-', 'Color', style.green, 'LineWidth', 1.5);
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Explained variance (%)', 'Interpreter', 'none');
title('Per-scale PCA variance audit', 'Interpreter', 'none');
legend({'PC1','First 5 PCs','Selected PCs'}, 'Location', 'best', 'Box', 'off');
local_format_axis(gca, style);

files = local_export(fig, fullfile(outRoot, 'figures'), 'embedding_pca_dimension_and_variance_by_scale', style);
close(fig);
end

function files = local_plot_matrix_manifest(outRoot, style)
T = local_read_csv(fullfile(outRoot, 'embedding_matrix_manifest.csv'));
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1450 760]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
bar(T.finite_value_fraction, 'FaceColor', style.blue, 'EdgeColor', 'none');
ylim([0 1.02]);
set(gca, 'XTick', 1:height(T), 'XTickLabel', local_clean_label(T.matrix_name), ...
    'TickLabelInterpreter', 'none');
xtickangle(45);
ylabel('Finite value fraction', 'Interpreter', 'none');
title('Matrix finite-value audit', 'Interpreter', 'none');
local_format_axis(gca, style);

nexttile;
scatter(T.n_rows, T.n_columns, 80, style.purple, 'filled');
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Rows', 'Interpreter', 'none');
ylabel('Columns', 'Interpreter', 'none');
title('Matrix size manifest', 'Interpreter', 'none');
local_format_axis(gca, style);

files = local_export(fig, fullfile(outRoot, 'figures'), 'embedding_matrix_manifest_finite_audit', style);
close(fig);
end

function files = local_plot_preprocess_dimension_audit(outRoot, style)
T = local_read_csv(fullfile(outRoot, 'embedding_preprocess_dimension_audit.csv'));
scales = unique(T.chunk_sec, 'stable');
nGuard = zeros(numel(scales), 1);
maxStd = nan(numel(scales), 1);
maxTailFraction = nan(numel(scales), 1);
for i = 1:numel(scales)
    idx = T.chunk_sec == scales(i);
    nGuard(i) = nnz(logical(double(T.sparse_scale_safeguard_triggered(idx))));
    maxStd(i) = max(T.standardized_std(idx), [], 'omitnan');
    maxTailFraction(i) = max(T.fraction_abs_gt_tail_threshold(idx), [], 'omitnan');
end

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1500 760]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(scales, nGuard, 'o-', 'Color', style.red, 'MarkerFaceColor', style.red, 'LineWidth', 1.6);
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Dimensions using sparse-scale guard', 'Interpreter', 'none');
title('Configurable IQR-to-SD safeguard', 'Interpreter', 'none');
local_format_axis(gca, style);

nexttile;
yyaxis left;
plot(scales, maxStd, 's-', 'Color', style.blue, 'LineWidth', 1.5);
ylabel('Maximum standardized SD', 'Interpreter', 'none');
yyaxis right;
plot(scales, maxTailFraction, 'd-', 'Color', style.gold, 'LineWidth', 1.5);
ylabel('Maximum fraction beyond tail threshold', 'Interpreter', 'none');
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
title('Post-standardization tail diagnostics', 'Interpreter', 'none');
local_format_axis(gca, style);

files = local_export(fig, fullfile(outRoot, 'figures'), ...
    'embedding_preprocess_sparse_scale_audit', style);
close(fig);
end

function files = local_plot_anchor_stage_sensitivity(outRoot, style)
T = local_read_csv(fullfile(outRoot, 'embedding_anchor_stage_pca_sensitivity_audit.csv'));
complete = string(T.audit_status) == "complete";
if ~any(complete)
    files = strings(0, 1);
    return
end
T = T(complete, :);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1450 720]);
hold on;
plot(T.chunk_sec, T.combined_vs_base_subspace_similarity, 'o-', ...
    'Color', style.blue, 'LineWidth', 1.5);
plot(T.chunk_sec, T.combined_vs_rare_enriched_subspace_similarity, 's-', ...
    'Color', style.red, 'LineWidth', 1.5);
plot(T.chunk_sec, T.base_vs_rare_enriched_subspace_similarity, 'd-', ...
    'Color', style.gold, 'LineWidth', 1.5);
set(gca, 'XScale', 'log');
ylim([0 1.02]);
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Mean principal-angle cosine', 'Interpreter', 'none');
title('Anchor-stage PCA sensitivity (audit only)', 'Interpreter', 'none');
legend({'Combined vs base','Combined vs rare enriched','Base vs rare enriched'}, ...
    'Location', 'best', 'Box', 'off');
local_format_axis(gca, style);
files = local_export(fig, fullfile(outRoot, 'figures'), ...
    'embedding_anchor_stage_pca_sensitivity', style);
close(fig);
end

function files = local_plot_scale_weights(outRoot, style)
T = local_read_csv(fullfile(outRoot, 'embedding_scale_weight_audit.csv'));
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1350 700]);
yyaxis left;
plot(T.chunk_sec, T.n_pcs_weighted, 'o-', 'Color', style.blue, 'LineWidth', 1.6);
ylabel('Weighted PCs', 'Interpreter', 'none');
yyaxis right;
plot(T.chunk_sec, T.scale_weight, 's-', 'Color', style.gold, 'LineWidth', 1.6);
ylabel('Scale weight', 'Interpreter', 'none');
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
title('Predefined scale weights before global PCA', 'Interpreter', 'none');
local_format_axis(gca, style);
files = local_export(fig, fullfile(outRoot, 'figures'), 'embedding_scale_weight_audit', style);
close(fig);
end

function files = local_plot_stability(outRoot, style)
T = local_read_csv(fullfile(outRoot, 'embedding_stability_audit.csv'));
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1450 720]);
hold on;
roles = unique(T.hierarchical_role, 'stable');
colors = [style.red; style.green; style.purple; style.blue];
for i = 1:numel(roles)
    idx = T.hierarchical_role == roles(i);
    c = colors(1 + mod(i - 1, size(colors, 1)), :);
    plot(T.chunk_sec(idx), T.median_subspace_similarity(idx), 'o-', ...
        'Color', c, 'MarkerFaceColor', c, 'LineWidth', 1.4);
end
set(gca, 'XScale', 'log');
ylim([0 1.02]);
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Median split-half subspace similarity', 'Interpreter', 'none');
title('Condition-blind PCA stability audit', 'Interpreter', 'none');
legend(cellstr(local_clean_label(roles)), 'Location', 'best', 'Box', 'off', 'Interpreter', 'none');
local_format_axis(gca, style);
files = local_export(fig, fullfile(outRoot, 'figures'), 'embedding_split_half_pca_stability', style);
close(fig);
end

function files = local_plot_arena_sensitivity(outRoot, style)
T = local_read_csv(fullfile(outRoot, 'embedding_arena_sensitivity_audit.csv'));
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1450 720]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(T.chunk_sec, T.embedding_arena_shift_iqr_units, 'o-', ...
    'Color', style.gold, 'MarkerFaceColor', style.gold, 'LineWidth', 1.6);
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Arena shift (IQR units)', 'Interpreter', 'none');
title('Arena sensitivity after fitting (audit only)', 'Interpreter', 'none');
local_format_axis(gca, style);

nexttile;
bar(1:height(T), T.top_base_feature_loading_fraction, ...
    'FaceColor', style.purple, 'EdgeColor', 'none');
set(gca, 'XTick', 1:height(T), ...
    'XTickLabel', local_clean_label(T.top_base_feature_by_loading), ...
    'TickLabelInterpreter', 'none');
ylabel('Top feature loading fraction', 'Interpreter', 'none');
title('Top loading feature by scale', 'Interpreter', 'none');
xtickangle(45);
local_format_axis(gca, style);

files = local_export(fig, fullfile(outRoot, 'figures'), 'embedding_arena_sensitivity_audit_only', style);
close(fig);
end

function files = local_plot_global_embedding(outRoot, style)
T = local_read_csv(fullfile(outRoot, 'embedding_global_scores.csv'));
if ~all(ismember(["global_pc01", "global_pc02"], string(T.Properties.VariableNames)))
    files = strings(0, 1);
    return
end
rng(107);
n = height(T);
if n > 7000
    idx = sort(randperm(n, 7000));
else
    idx = 1:n;
end
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1500 720]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
scatter(T.global_pc01(idx), T.global_pc02(idx), 12, double(T.chunk_sec(idx)), 'filled');
xlabel('Global PC1', 'Interpreter', 'none');
ylabel('Global PC2', 'Interpreter', 'none');
title('Global embedding colored by scale', 'Interpreter', 'none');
colorbar;
local_format_axis(gca, style);

nexttile;
arenaLabel = local_arena_labels_for_scores(outRoot, T);
arena = categorical(arenaLabel(idx));
gscatter(T.global_pc01(idx), T.global_pc02(idx), arena);
xlabel('Global PC1', 'Interpreter', 'none');
ylabel('Global PC2', 'Interpreter', 'none');
title('Arena overlay after fitting (audit only)', 'Interpreter', 'none');
local_format_axis(gca, style);

files = local_export(fig, fullfile(outRoot, 'figures'), 'embedding_global_pc_overview_audit', style);
close(fig);
end

function arenaLabel = local_arena_labels_for_scores(outRoot, scoreTable)
if ismember('arena_label', scoreTable.Properties.VariableNames)
    arenaLabel = string(scoreTable.arena_label);
    return
end
manifestPath = fullfile(outRoot, 'embedding_row_manifest.csv');
R = local_read_csv(manifestPath);
assert(all(ismember(["embedding_row_id", "arena_label"], string(R.Properties.VariableNames))), ...
    'make_run07_embedding_qc_figures:MissingArenaAuditManifest', ...
    'Global arena audit figure requires embedding_row_manifest.csv with embedding_row_id and arena_label.');
[tf, loc] = ismember(scoreTable.embedding_row_id, R.embedding_row_id);
assert(all(tf), 'make_run07_embedding_qc_figures:ScoreRowsMissingManifest', ...
    'Global score rows could not be matched to the arena audit row manifest.');
arenaLabel = string(R.arena_label(loc));
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
