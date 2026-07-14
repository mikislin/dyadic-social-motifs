function figureManifest = make_run06_chunk_qc_figures(outRoot, params)
%MAKE_RUN06_CHUNK_QC_FIGURES Generate CSV-backed run-06 QC figures.

outRoot = string(outRoot);
figDir = fullfile(outRoot, 'figures');
sourceDir = fullfile(outRoot, 'figure_sources');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

style = local_style(params);
rows = table();

rows = local_add_figure(rows, local_plot_coverage(outRoot, style), ...
    fullfile(outRoot, 'chunk_bank_summary.csv'), ...
    "chunk coverage by temporal scale");
rows = local_add_figure(rows, local_plot_scale_anchor_survey(outRoot, style), ...
    fullfile(outRoot, 'scale_anchor_coverage_audit.csv'), ...
    "scale-specific anchor coverage and all-scale retention audit");
rows = local_add_figure(rows, local_plot_session_scale_counts(sourceDir, style), ...
    fullfile(sourceDir, 'chunk_valid_counts_by_session_scale.csv'), ...
    "valid chunk counts by session and scale");
rows = local_add_figure(rows, local_plot_arena_qc_counts(sourceDir, style), ...
    fullfile(sourceDir, 'chunk_valid_counts_by_arena_qc.csv'), ...
    "valid chunk counts stratified by arena and QC set for audit only");
rows = local_add_figure(rows, local_plot_scale_selection(outRoot, style), ...
    fullfile(outRoot, 'scale_selection_audit.csv'), ...
    "scale usefulness criteria and selected operational scales");
rows = local_add_figure(rows, local_plot_primary_scale_decision(outRoot, style), ...
    fullfile(outRoot, 'primary_operational_scales.csv'), ...
    "primary operational scales after stability and dimensionality guards");
rows = local_add_figure(rows, local_plot_motif_band_dimension_stability(outRoot, style), ...
    string(fullfile(outRoot, 'selected_operational_scales.csv')) + "; " + ...
    string(fullfile(outRoot, 'embedding_dimension_audit.csv')), ...
    "motif-band PC1 PC2 effective dimension stability and primary status audit");
rows = local_add_figure(rows, local_plot_embedding_dimension_audit(outRoot, style), ...
    fullfile(outRoot, 'embedding_dimension_audit.csv'), ...
    "multiresolution representation compression and PCA dimensionality audit");
rows = local_add_figure(rows, local_plot_scale_selection_stability(outRoot, style), ...
    fullfile(outRoot, 'scale_selection_stability.csv'), ...
    "condition-blind bootstrap stability of selected operational scales");
rows = local_add_figure(rows, local_plot_pca_loading_stability(outRoot, style), ...
    fullfile(outRoot, 'pca_loading_stability.csv'), ...
    "condition-blind split-half PCA loading stability by scale");
rows = local_add_figure(rows, local_plot_primary_scale_specific_coverage(outRoot, style), ...
    fullfile(outRoot, 'primary_scale_specific_chunk_bank_summary.csv'), ...
    "scale-specific anchor coverage for primary operational scales");
rows = local_add_figure(rows, local_plot_primary_event_summary(sourceDir, style), ...
    fullfile(sourceDir, 'primary_chunk_event_summary_source.csv'), ...
    "within-chunk event summaries for primary operational scales");
rows = local_add_figure(rows, local_plot_validity(outRoot, style), ...
    fullfile(outRoot, 'chunk_validity_summary.csv'), ...
    "chunk validity and frame-mask propagation by scale");
rows = local_add_figure(rows, local_plot_examples(sourceDir, style), ...
    fullfile(sourceDir, 'chunk_example_trace_source.csv'), ...
    "example chunk traces across selected scales");
rows = local_add_figure(rows, local_plot_transform_audit(outRoot, style), ...
    fullfile(outRoot, 'chunk_feature_transform_audit.csv'), ...
    "feature and channel transform audit");
rows = local_add_figure(rows, local_plot_arena_shift(sourceDir, style), ...
    fullfile(sourceDir, 'chunk_arena_shift_after_transform.csv'), ...
    "arena shift diagnostics after chunk transform and scaling for audit only");
rows = local_add_figure(rows, local_plot_arena_embedding_sensitivity(sourceDir, style), ...
    fullfile(sourceDir, 'arena_embedding_sensitivity_source.csv'), ...
    "arena sensitivity and feature-family loading audit after condition-blind fitting");

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

function files = local_plot_coverage(outRoot, style)
T = readtable(fullfile(outRoot, 'chunk_bank_summary.csv'), 'TextType', 'string');
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1350 720]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
if ismember('n_scale_specific_candidate_anchors', T.Properties.VariableNames)
    plot(T.chunk_sec, T.n_scale_specific_candidate_anchors, 'o-', 'Color', style.purple, 'LineWidth', 1.6);
end
plot(T.chunk_sec, T.n_candidate_anchors_before_cap, 'o-', 'Color', style.gray, 'LineWidth', 1.6);
plot(T.chunk_sec, T.n_chunks, 'o-', 'Color', style.blue, 'LineWidth', 1.8);
plot(T.chunk_sec, T.n_valid_chunks, 'o-', 'Color', style.green, 'LineWidth', 1.8);
hold off;
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Chunk or anchor count', 'Interpreter', 'none');
title('Coverage by temporal scale', 'Interpreter', 'none');
if ismember('n_scale_specific_candidate_anchors', T.Properties.VariableNames)
    legend({'Scale-specific candidates', 'All-scale candidates', 'Materialized chunks', 'Valid chunks'}, ...
        'Location', 'best', 'Interpreter', 'none', 'Box', 'off');
    local_expand_count_axis(gca, [T.n_scale_specific_candidate_anchors; T.n_candidate_anchors_before_cap; T.n_chunks; T.n_valid_chunks]);
else
    legend({'All-scale candidates', 'Materialized chunks', 'Valid chunks'}, ...
        'Location', 'best', 'Interpreter', 'none', 'Box', 'off');
    local_expand_count_axis(gca, [T.n_candidate_anchors_before_cap; T.n_chunks; T.n_valid_chunks]);
end
local_format_axis(gca, style);

nexttile;
plot(T.chunk_sec, T.n_sessions_with_valid_chunks, 'o-', 'Color', style.green, 'LineWidth', 1.8);
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Sessions with valid chunks', 'Interpreter', 'none');
title('Session support by scale', 'Interpreter', 'none');
local_expand_count_axis(gca, T.n_sessions_with_valid_chunks);
local_format_axis(gca, style);

files = local_export(fig, fullfile(outRoot, 'figures'), 'chunk_coverage_by_temporal_scale', style);
close(fig);
end

function files = local_plot_scale_anchor_survey(outRoot, style)
path = fullfile(outRoot, 'scale_anchor_coverage_audit.csv');
if ~isfile(path)
    files = strings(0, 1);
    return
end
T = readtable(path, 'TextType', 'string');
if isempty(T)
    files = strings(0, 1);
    return
end
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1450 760]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
plot(T.chunk_sec, T.n_scale_specific_candidate_anchors, 'o-', 'Color', style.purple, 'LineWidth', 1.7);
plot(T.chunk_sec, T.n_common_all_scale_candidate_anchors, 'o-', 'Color', style.gray, 'LineWidth', 1.7);
hold off;
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Candidate anchors', 'Interpreter', 'none');
title('Scale-specific versus all-scale anchors', 'Interpreter', 'none');
legend({'Scale-specific candidates', 'All-scale common candidates'}, ...
    'Location', 'best', 'Interpreter', 'none', 'Box', 'off');
local_expand_count_axis(gca, [T.n_scale_specific_candidate_anchors; T.n_common_all_scale_candidate_anchors]);
local_format_axis(gca, style);

nexttile;
hold on;
plot(T.chunk_sec, T.n_sessions_with_scale_specific_anchors, 'o-', 'Color', style.purple, 'LineWidth', 1.7);
plot(T.chunk_sec, T.n_sessions_with_common_all_scale_anchors, 'o-', 'Color', style.gray, 'LineWidth', 1.7);
yyaxis right;
plot(T.chunk_sec, T.common_anchor_retention_fraction, 's--', 'Color', style.red, 'LineWidth', 1.4);
ylabel('Common-anchor retention', 'Interpreter', 'none');
ylim([0 1.02]);
yyaxis left;
hold off;
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Sessions with anchors', 'Interpreter', 'none');
title('Long-context attrition audit', 'Interpreter', 'none');
legend({'Scale-specific sessions', 'All-scale common sessions', 'Retention fraction'}, ...
    'Location', 'best', 'Interpreter', 'none', 'Box', 'off');
local_format_axis(gca, style);

files = local_export(fig, fullfile(outRoot, 'figures'), 'scale_specific_anchor_coverage_audit', style);
close(fig);
end

function files = local_plot_session_scale_counts(sourceDir, style)
T = readtable(fullfile(sourceDir, 'chunk_valid_counts_by_session_scale.csv'), 'TextType', 'string');
scales = unique(T.scale_index, 'stable');
sessions = unique(T.raw_index, 'stable');
M = zeros(numel(sessions), numel(scales));
for i = 1:numel(sessions)
    for s = 1:numel(scales)
        idx = T.raw_index == sessions(i) & T.scale_index == scales(s);
        if any(idx)
            M(i, s) = sum(T.n_valid_chunks(idx));
        end
    end
end

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1500 780]);
imagesc(M);
colormap(parula);
colorbar;
set(gca, 'XTick', 1:numel(scales), 'XTickLabel', string(scales), ...
    'YTick', 1:numel(sessions), 'YTickLabel', string(sessions), ...
    'TickLabelInterpreter', 'none');
xlabel('Scale index', 'Interpreter', 'none');
ylabel('Raw session index', 'Interpreter', 'none');
title('Valid chunk counts by session and scale', 'Interpreter', 'none');
local_format_axis(gca, style);
files = local_export(fig, fullfile(fileparts(sourceDir), 'figures'), 'valid_chunk_counts_by_session_and_scale', style);
close(fig);
end

function files = local_plot_arena_qc_counts(sourceDir, style)
T = readtable(fullfile(sourceDir, 'chunk_valid_counts_by_arena_qc.csv'), 'TextType', 'string');
T.group_label = local_clean_label(T.arena_label + " | " + T.qc_set);
groups = unique(T.group_label, 'stable');
scales = unique(T.scale_index, 'stable');
Y = zeros(numel(scales), numel(groups));
for s = 1:numel(scales)
    for g = 1:numel(groups)
        idx = T.scale_index == scales(s) & T.group_label == groups(g);
        Y(s, g) = sum(T.n_valid_chunks(idx));
    end
end

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1450 760]);
bar(Y, 'stacked', 'EdgeColor', 'none');
set(gca, 'XTick', 1:numel(scales), 'XTickLabel', string(scales), ...
    'TickLabelInterpreter', 'none');
xlabel('Scale index', 'Interpreter', 'none');
ylabel('Valid chunks', 'Interpreter', 'none');
title('Arena and QC-set coverage audit only', 'Interpreter', 'none');
legend(cellstr(groups), 'Location', 'eastoutside', 'Interpreter', 'none', 'Box', 'off');
local_format_axis(gca, style);
files = local_export(fig, fullfile(fileparts(sourceDir), 'figures'), 'valid_chunk_counts_by_arena_and_qc_audit_only', style);
close(fig);
end

function files = local_plot_scale_selection(outRoot, style)
T = readtable(fullfile(outRoot, 'scale_selection_audit.csv'), 'TextType', 'string');
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1500 820]);
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
plot(T.chunk_sec, T.composite_micro, 'o-', 'Color', style.red, 'LineWidth', 1.5);
plot(T.chunk_sec, T.composite_motif, 's-', 'Color', style.green, 'LineWidth', 1.5);
plot(T.chunk_sec, T.composite_context, 'd-', 'Color', style.purple, 'LineWidth', 1.5);
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Composite score', 'Interpreter', 'none');
title('Role-specific usefulness', 'Interpreter', 'none');
legend({'Micro','Motif','Context'}, 'Location', 'best', 'Box', 'off');
local_format_axis(gca, style);

nexttile;
hold on;
plot(T.chunk_sec, T.composite_global, 'k.-', 'LineWidth', 1.6, 'MarkerSize', 12);
sel = logical(T.selected_operational_scale);
scatter(T.chunk_sec(sel), T.composite_global(sel), 80, style.gold, 'filled', 'MarkerEdgeColor', 'k');
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Global score', 'Interpreter', 'none');
title('Selected operational scales', 'Interpreter', 'none');
local_format_axis(gca, style);

nexttile;
plot(T.chunk_sec, T.lag1_embedding_corr, 'o-', 'LineWidth', 1.4);
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Lag-1 embedding correlation', 'Interpreter', 'none');
title('Temporal persistence criterion', 'Interpreter', 'none');
local_format_axis(gca, style);

nexttile;
if ismember('cum_embedding_pcs_explained', T.Properties.VariableNames)
    yyaxis left;
    plot(T.chunk_sec, T.cum_embedding_pcs_explained, 'o-', 'LineWidth', 1.4);
    ylabel('Scoring PCs explained (%)', 'Interpreter', 'none');
    yyaxis right;
    if ismember('n_pcs_90pct', T.Properties.VariableNames)
        plot(T.chunk_sec, T.n_pcs_90pct, 's--', 'LineWidth', 1.2);
        ylabel('PCs for 90% variance', 'Interpreter', 'none');
    else
        plot(T.chunk_sec, T.effective_dim, 's--', 'LineWidth', 1.2);
        ylabel('Effective dimensionality', 'Interpreter', 'none');
    end
else
    plot(T.chunk_sec, T.cum5_explained, 'o-', 'LineWidth', 1.4);
    ylabel('First 5 PCs explained (%)', 'Interpreter', 'none');
end
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
title('Representation and dimensionality audit', 'Interpreter', 'none');
local_format_axis(gca, style);

files = local_export(fig, fullfile(outRoot, 'figures'), 'scale_usefulness_and_selected_operational_scales', style);
close(fig);
end

function files = local_plot_primary_scale_decision(outRoot, style)
selectedPath = fullfile(outRoot, 'selected_operational_scales.csv');
primaryPath = fullfile(outRoot, 'primary_operational_scales.csv');
if ~isfile(selectedPath) || ~isfile(primaryPath)
    files = strings(0, 1);
    return
end
S = readtable(selectedPath, 'TextType', 'string');
P = readtable(primaryPath, 'TextType', 'string');
if isempty(S)
    files = strings(0, 1);
    return
end
isPrimary = ismember(S.scale_index, P.scale_index);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1450 760]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
scatter(S.chunk_sec, S.bootstrap_selection_frequency, 70, style.gray, 'filled');
scatter(S.chunk_sec(isPrimary), S.bootstrap_selection_frequency(isPrimary), 105, style.green, ...
    'filled', 'MarkerEdgeColor', 'k');
if ismember('passes_stability_threshold', S.Properties.VariableNames)
    yline(0.7, 'k--', 'LineWidth', 1.1);
end
set(gca, 'XScale', 'log');
ylim([0 1.02]);
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Bootstrap selection frequency', 'Interpreter', 'none');
title('Primary scales pass selection stability', 'Interpreter', 'none');
legend({'Score-selected','Primary scale','Stability threshold'}, ...
    'Location', 'best', 'Interpreter', 'none', 'Box', 'off');
local_format_axis(gca, style);

nexttile;
hold on;
scatter(S.effective_dim, S.pc1_explained, 70, style.gray, 'filled');
scatter(S.effective_dim(isPrimary), S.pc1_explained(isPrimary), 105, style.green, ...
    'filled', 'MarkerEdgeColor', 'k');
xline(S.dimension_guard_min_effective_dim(1), 'k--', 'LineWidth', 1.1);
yline(S.dimension_guard_max_pc1_explained_pct(1), 'k--', 'LineWidth', 1.1);
xlabel('Effective PCA dimensionality', 'Interpreter', 'none');
ylabel('PC1 explained (%)', 'Interpreter', 'none');
title('Primary scales pass dimensionality guard', 'Interpreter', 'none');
legend({'Score-selected','Primary scale','Guard threshold'}, ...
    'Location', 'best', 'Interpreter', 'none', 'Box', 'off');
local_format_axis(gca, style);

files = local_export(fig, fullfile(outRoot, 'figures'), 'primary_operational_scale_decision', style);
close(fig);
end

function files = local_plot_motif_band_dimension_stability(outRoot, style)
path = fullfile(outRoot, 'selected_operational_scales.csv');
if ~isfile(path)
    files = strings(0, 1);
    return
end
T = readtable(path, 'TextType', 'string');
if isempty(T) || ~ismember('initial_band', T.Properties.VariableNames)
    files = strings(0, 1);
    return
end
M = T(string(T.initial_band) == "motif", :);
if isempty(M)
    files = strings(0, 1);
    return
end
M = sortrows(M, 'chunk_sec');
pc2 = nan(height(M), 1);
if ismember('pc2_explained', M.Properties.VariableNames)
    pc2 = M.pc2_explained;
end
primaryStatus = false(height(M), 1);
if ismember('stable_and_dimension_supported', M.Properties.VariableNames)
    primaryStatus = logical(M.stable_and_dimension_supported);
end
stability = nan(height(M), 1);
if ismember('bootstrap_selection_frequency', M.Properties.VariableNames)
    stability = M.bootstrap_selection_frequency;
end
rankWithin = nan(height(M), 1);
if ismember('rank_within_role', M.Properties.VariableNames)
    rankWithin = M.rank_within_role;
end

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1450 900]);
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
plot(M.chunk_sec, M.pc1_explained, 'o-', 'Color', style.blue, 'LineWidth', 1.6);
plot(M.chunk_sec, pc2, 's-', 'Color', style.gold, 'LineWidth', 1.6);
scatter(M.chunk_sec(primaryStatus), M.pc1_explained(primaryStatus), 80, ...
    style.green, 'filled', 'MarkerEdgeColor', 'k');
hold off;
set(gca, 'XScale', 'log');
xlabel('Motif-band scale (s)', 'Interpreter', 'none');
ylabel('Variance explained (%)', 'Interpreter', 'none');
title('Motif PCA concentration', 'Interpreter', 'none');
legend({'PC1', 'PC2', 'Primary scale'}, 'Location', 'best', 'Box', 'off');
local_format_axis(gca, style);

nexttile;
plot(M.chunk_sec, M.effective_dim, 'o-', 'Color', style.purple, 'LineWidth', 1.7);
hold on;
if ismember('dimension_guard_min_effective_dim', M.Properties.VariableNames)
    yline(M.dimension_guard_min_effective_dim(1), 'k--', 'LineWidth', 1.1);
end
scatter(M.chunk_sec(primaryStatus), M.effective_dim(primaryStatus), 80, ...
    style.green, 'filled', 'MarkerEdgeColor', 'k');
hold off;
set(gca, 'XScale', 'log');
xlabel('Motif-band scale (s)', 'Interpreter', 'none');
ylabel('Effective PCA dimension', 'Interpreter', 'none');
title('Motif dimensional support', 'Interpreter', 'none');
legend({'Effective dimension', 'Dimension guard', 'Primary scale'}, ...
    'Location', 'best', 'Box', 'off');
local_format_axis(gca, style);

nexttile;
plot(M.chunk_sec, stability, 'o-', 'Color', style.green, 'LineWidth', 1.7);
hold on;
if ismember('stability_selection_frequency_threshold', M.Properties.VariableNames)
    yline(M.stability_selection_frequency_threshold(1), 'k--', 'LineWidth', 1.1);
elseif ismember('passes_stability_threshold', M.Properties.VariableNames)
    yline(0.7, 'k--', 'LineWidth', 1.1);
end
scatter(M.chunk_sec(primaryStatus), stability(primaryStatus), 80, ...
    style.green, 'filled', 'MarkerEdgeColor', 'k');
hold off;
set(gca, 'XScale', 'log');
ylim([0 1.05]);
xlabel('Motif-band scale (s)', 'Interpreter', 'none');
ylabel('Bootstrap selection frequency', 'Interpreter', 'none');
title('Motif selection stability', 'Interpreter', 'none');
legend({'Selection frequency', 'Stability threshold', 'Primary scale'}, ...
    'Location', 'best', 'Box', 'off');
local_format_axis(gca, style);

nexttile;
bar(M.chunk_sec, double(primaryStatus), 0.6, 'FaceColor', style.gray, 'EdgeColor', 'none');
hold on;
scatter(M.chunk_sec, double(primaryStatus), 70, rankWithin, 'filled', 'MarkerEdgeColor', 'k');
hold off;
set(gca, 'XScale', 'log');
ylim([-0.05 1.2]);
yticks([0 1]);
yticklabels({'survey only','primary'});
xlabel('Motif-band scale (s)', 'Interpreter', 'none');
ylabel('Primary status', 'Interpreter', 'none');
title('Motif primary-promotion status', 'Interpreter', 'none');
cb = colorbar;
cb.Label.String = 'Rank within motif role';
local_format_axis(gca, style);

files = local_export(fig, fullfile(outRoot, 'figures'), 'motif_band_dimension_stability_audit', style);
close(fig);
end

function files = local_plot_embedding_dimension_audit(outRoot, style)
path = fullfile(outRoot, 'embedding_dimension_audit.csv');
if ~isfile(path)
    files = strings(0, 1);
    return
end
T = readtable(path, 'TextType', 'string');
if isempty(T)
    files = strings(0, 1);
    return
end

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1500 780]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
plot(T.chunk_sec, T.n_raw_flattened_dims, 'o-', 'Color', style.gray, 'LineWidth', 1.6);
plot(T.chunk_sec, T.n_summary_dims, 's-', 'Color', style.blue, 'LineWidth', 1.8);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Representation dimensions', 'Interpreter', 'none');
title('Raw flattened versus multiresolution dimensions', 'Interpreter', 'none');
legend({'Raw frame x channel dimensions','Summary dimensions'}, ...
    'Location', 'best', 'Interpreter', 'none', 'Box', 'off');
local_format_axis(gca, style);

nexttile;
hold on;
plot(T.chunk_sec, T.n_pcs_90pct, 'o-', 'Color', style.purple, 'LineWidth', 1.6);
plot(T.chunk_sec, T.n_score_pcs_retained, 's--', 'Color', style.gray, 'LineWidth', 1.3);
if ismember('recommended_pcs_for_next_embedding', T.Properties.VariableNames)
    plot(T.chunk_sec, T.recommended_pcs_for_next_embedding, 'd-', 'Color', style.green, 'LineWidth', 1.6);
end
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('PCA components', 'Interpreter', 'none');
title('PCA dimensionality on compact summaries', 'Interpreter', 'none');
legend({'PCs for 90% variance','Scoring PCs retained','Recommended cap for next embedding'}, ...
    'Location', 'best', 'Interpreter', 'none', 'Box', 'off');
local_format_axis(gca, style);

files = local_export(fig, fullfile(outRoot, 'figures'), 'embedding_dimension_audit_multiresolution', style);
close(fig);
end

function files = local_plot_pca_loading_stability(outRoot, style)
path = fullfile(outRoot, 'pca_loading_stability.csv');
if ~isfile(path)
    files = strings(0, 1);
    return
end
T = readtable(path, 'TextType', 'string');
if isempty(T)
    files = strings(0, 1);
    return
end

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1450 720]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
roles = unique(T.initial_band, 'stable');
colors = [style.red; style.green; style.purple; style.blue];
for i = 1:numel(roles)
    idx = T.initial_band == roles(i);
    c = colors(1 + mod(i - 1, size(colors, 1)), :);
    plot(T.chunk_sec(idx), T.median_subspace_similarity(idx), 'o-', ...
        'Color', c, 'MarkerFaceColor', c, 'LineWidth', 1.4);
end
yline(T.loading_stability_threshold(1), 'k--', 'LineWidth', 1.1);
set(gca, 'XScale', 'log');
ylim([0 1.02]);
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Median subspace similarity', 'Interpreter', 'none');
title('Split-half PCA loading stability', 'Interpreter', 'none');
legend([cellstr(local_clean_label(roles)); {'Threshold'}], ...
    'Location', 'best', 'Interpreter', 'none', 'Box', 'off');
local_format_axis(gca, style);

nexttile;
hold on;
plot(T.chunk_sec, T.median_pc1_abs_correlation, 'o-', 'Color', style.blue, 'LineWidth', 1.5);
plot(T.chunk_sec, T.pc1_abs_correlation_p10, ':', 'Color', style.gray, 'LineWidth', 1.1);
plot(T.chunk_sec, T.pc1_abs_correlation_p90, ':', 'Color', style.gray, 'LineWidth', 1.1);
set(gca, 'XScale', 'log');
ylim([0 1.02]);
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('PC1 absolute loading correlation', 'Interpreter', 'none');
title('PC1 loading agreement interval', 'Interpreter', 'none');
legend({'Median','P10/P90 interval'}, 'Location', 'best', 'Interpreter', 'none', 'Box', 'off');
local_format_axis(gca, style);

files = local_export(fig, fullfile(outRoot, 'figures'), 'pca_loading_stability', style);
close(fig);
end

function files = local_plot_primary_scale_specific_coverage(outRoot, style)
path = fullfile(outRoot, 'primary_scale_specific_chunk_bank_summary.csv');
if ~isfile(path)
    files = strings(0, 1);
    return
end
T = readtable(path, 'TextType', 'string');
if isempty(T)
    files = strings(0, 1);
    return
end

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1450 720]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
plot(T.chunk_sec, T.n_scale_specific_candidate_anchors, 'o-', 'Color', style.purple, 'LineWidth', 1.5);
plot(T.chunk_sec, T.n_primary_scale_specific_anchors, 's-', 'Color', style.green, 'LineWidth', 1.7);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Anchor count', 'Interpreter', 'none');
title('Primary scale-specific anchor bank (log count)', 'Interpreter', 'none');
legend({'Candidate anchors','Retained primary anchors'}, ...
    'Location', 'best', 'Interpreter', 'none', 'Box', 'off');
local_expand_count_axis(gca, [T.n_scale_specific_candidate_anchors; T.n_primary_scale_specific_anchors]);
local_format_axis(gca, style);

nexttile;
hold on;
plot(T.chunk_sec, T.n_sessions_with_candidates, 'o-', 'Color', style.gray, 'LineWidth', 1.5);
plot(T.chunk_sec, T.n_sessions_with_primary_anchors, 's-', 'Color', style.green, 'LineWidth', 1.7);
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Sessions', 'Interpreter', 'none');
title('Session support for primary anchors', 'Interpreter', 'none');
legend({'Candidate support','Retained-anchor support'}, ...
    'Location', 'best', 'Interpreter', 'none', 'Box', 'off');
local_expand_count_axis(gca, [T.n_sessions_with_candidates; T.n_sessions_with_primary_anchors]);
local_format_axis(gca, style);

files = local_export(fig, fullfile(outRoot, 'figures'), 'primary_scale_specific_anchor_coverage', style);
close(fig);
end

function files = local_plot_scale_selection_stability(outRoot, style)
path = fullfile(outRoot, 'scale_selection_stability.csv');
if ~isfile(path)
    files = strings(0, 1);
    return
end
T = readtable(path, 'TextType', 'string');
if isempty(T)
    files = strings(0, 1);
    return
end

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1500 780]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
roles = unique(T.initial_band, 'stable');
colors = [style.red; style.green; style.purple; style.blue];
for i = 1:numel(roles)
    idx = T.initial_band == roles(i);
    c = colors(1 + mod(i - 1, size(colors, 1)), :);
    plot(T.chunk_sec(idx), T.selection_frequency(idx), 'o-', ...
        'Color', c, 'MarkerFaceColor', c, 'LineWidth', 1.4);
end
if ismember('stability_threshold', T.Properties.VariableNames)
    yline(T.stability_threshold(1), 'k--', 'LineWidth', 1.2);
end
sel = logical(T.selected_in_original_run);
scatter(T.chunk_sec(sel), T.selection_frequency(sel), 90, style.gold, ...
    'filled', 'MarkerEdgeColor', 'k');
set(gca, 'XScale', 'log');
ylim([0 1.02]);
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Bootstrap selection frequency', 'Interpreter', 'none');
title('Condition-blind scale-selection stability', 'Interpreter', 'none');
legend([cellstr(local_clean_label(roles)); {'Threshold'; 'Original selected'}], ...
    'Location', 'best', 'Interpreter', 'none', 'Box', 'off');
local_format_axis(gca, style);

nexttile;
plot(T.chunk_sec, T.median_role_rank, 'o-', 'Color', style.blue, 'LineWidth', 1.5);
hold on;
plot(T.chunk_sec, T.role_rank_p10, ':', 'Color', style.gray, 'LineWidth', 1.1);
plot(T.chunk_sec, T.role_rank_p90, ':', 'Color', style.gray, 'LineWidth', 1.1);
set(gca, 'XScale', 'log', 'YDir', 'reverse');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Role rank', 'Interpreter', 'none');
title('Bootstrap role-rank interval', 'Interpreter', 'none');
legend({'Median rank','P10/P90 interval'}, ...
    'Location', 'best', 'Interpreter', 'none', 'Box', 'off');
local_format_axis(gca, style);

files = local_export(fig, fullfile(outRoot, 'figures'), 'scale_selection_bootstrap_stability', style);
close(fig);
end

function files = local_plot_primary_event_summary(sourceDir, style)
path = fullfile(sourceDir, 'primary_chunk_event_summary_source.csv');
if ~isfile(path)
    files = strings(0, 1);
    return
end
T = readtable(path, 'TextType', 'string');
if isempty(T)
    files = strings(0, 1);
    return
end
G = groupsummary(T, 'chunk_sec', {'median','mean'}, ...
    {'contact_dwell_fraction','contact_transition_count', ...
    'mutual_approach_dwell_fraction','withdrawal_dwell_fraction', ...
    'role_asymmetry_bias_mean','radial_speed_sign_change_count'});

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1450 820]);
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
plot(G.chunk_sec, G.median_contact_dwell_fraction, 'o-', 'Color', style.blue, 'LineWidth', 1.5);
plot(G.chunk_sec, G.median_mutual_approach_dwell_fraction, 's-', 'Color', style.green, 'LineWidth', 1.5);
plot(G.chunk_sec, G.median_withdrawal_dwell_fraction, 'd-', 'Color', style.red, 'LineWidth', 1.5);
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Median dwell fraction', 'Interpreter', 'none');
title('Event dwell summaries by primary scale', 'Interpreter', 'none');
legend({'Contact','Mutual approach','Withdrawal'}, 'Location', 'best', 'Interpreter', 'none', 'Box', 'off');
local_format_axis(gca, style);

nexttile;
plot(G.chunk_sec, G.median_contact_transition_count, 'o-', 'Color', style.purple, 'LineWidth', 1.5);
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Median transitions per chunk', 'Interpreter', 'none');
title('Contact transition counts', 'Interpreter', 'none');
local_format_axis(gca, style);

nexttile;
plot(G.chunk_sec, G.median_role_asymmetry_bias_mean, 'o-', 'Color', style.gold, 'LineWidth', 1.5);
yline(0, 'k--', 'LineWidth', 1.0);
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Median role asymmetry bias', 'Interpreter', 'none');
title('Role-asymmetry event summary', 'Interpreter', 'none');
local_format_axis(gca, style);

nexttile;
plot(G.chunk_sec, G.median_radial_speed_sign_change_count, 'o-', 'Color', style.gray, 'LineWidth', 1.5);
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Median sign changes', 'Interpreter', 'none');
title('Direction-change summary', 'Interpreter', 'none');
local_format_axis(gca, style);

files = local_export(fig, fullfile(fileparts(sourceDir), 'figures'), 'primary_chunk_event_summary_audit', style);
close(fig);
end

function files = local_plot_validity(outRoot, style)
T = readtable(fullfile(outRoot, 'chunk_validity_summary.csv'), 'TextType', 'string');
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1350 720]);
hold on;
plot(T.chunk_sec, T.mean_valid_frac, 'o-', 'Color', style.blue, 'LineWidth', 1.7);
plot(T.chunk_sec, T.min_valid_frac, 's--', 'Color', style.gray, 'LineWidth', 1.3);
plot(T.chunk_sec, T.mean_feature_finite_frac, 'd-', 'Color', style.green, 'LineWidth', 1.7);
set(gca, 'XScale', 'log');
ylim([0 1.02]);
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Fraction', 'Interpreter', 'none');
title('Frame-mask propagation and feature availability by scale', 'Interpreter', 'none');
legend({'Mean frame-valid fraction','Minimum frame-valid fraction','Mean finite feature fraction'}, ...
    'Location', 'southwest', 'Box', 'off');
local_format_axis(gca, style);
files = local_export(fig, fullfile(outRoot, 'figures'), 'chunk_validity_frame_mask_propagation_by_scale', style);
close(fig);
end

function files = local_plot_examples(sourceDir, style)
T = readtable(fullfile(sourceDir, 'chunk_example_trace_source.csv'), 'TextType', 'string');
if isempty(T)
    files = strings(0, 1);
    return
end
features = unique(T.feature_name, 'stable');
scales = unique(T.scale_index, 'stable');
nRows = numel(scales);
nCols = min(numel(features), 4);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1680 260 * max(nRows, 1)]);
tiledlayout(fig, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
for s = 1:nRows
    for f = 1:nCols
        nexttile;
        hold on;
        idx = T.scale_index == scales(s) & T.feature_name == features(f);
        S = T(idx, :);
        examples = unique(S.example_rank, 'stable');
        for e = examples(:)'
            E = S(S.example_rank == e, :);
            plot(E.time_within_chunk_s, E.value, 'LineWidth', 0.9);
        end
        if s == 1
            title(local_clean_label(features(f)), 'Interpreter', 'none');
        end
        if f == 1
            sec = T.chunk_sec(find(T.scale_index == scales(s), 1));
            ylabel(sprintf('Scale %d (%.4g s)', scales(s), sec), 'Interpreter', 'none');
        end
        if s == nRows
            xlabel('Time within chunk (s)', 'Interpreter', 'none');
        end
        local_format_axis(gca, style);
    end
end
files = local_export(fig, fullfile(fileparts(sourceDir), 'figures'), 'example_chunk_traces_across_selected_scales', style);
close(fig);
end

function files = local_plot_arena_embedding_sensitivity(sourceDir, style)
path = fullfile(sourceDir, 'arena_embedding_sensitivity_source.csv');
if ~isfile(path)
    files = strings(0, 1);
    return
end
T = readtable(path, 'TextType', 'string');
if isempty(T)
    files = strings(0, 1);
    return
end

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1450 760]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
plot(T.chunk_sec, T.distance_family_loading_fraction, 'o-', 'Color', style.purple, 'LineWidth', 1.5);
plot(T.chunk_sec, T.kinematics_family_loading_fraction, 's-', 'Color', style.blue, 'LineWidth', 1.5);
plot(T.chunk_sec, T.contact_family_loading_fraction, 'd-', 'Color', style.green, 'LineWidth', 1.5);
set(gca, 'XScale', 'log');
ylim([0 1]);
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('PCA loading energy fraction', 'Interpreter', 'none');
title('Feature-family loading audit', 'Interpreter', 'none');
legend({'Distance','Kinematics','Contact'}, 'Location', 'best', 'Interpreter', 'none', 'Box', 'off');
local_format_axis(gca, style);

nexttile;
yyaxis left;
plot(T.chunk_sec, T.embedding_arena_shift_iqr_units, 'o-', 'Color', style.gold, 'LineWidth', 1.6);
ylabel('Arena shift in score space (IQR units)', 'Interpreter', 'none');
yyaxis right;
plot(T.chunk_sec, T.distance_family_loading_fraction, 's--', 'Color', style.purple, 'LineWidth', 1.2);
ylabel('Distance loading fraction', 'Interpreter', 'none');
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
title('Arena sensitivity audit only', 'Interpreter', 'none');
local_format_axis(gca, style);

files = local_export(fig, fullfile(fileparts(sourceDir), 'figures'), 'arena_embedding_sensitivity_audit_only', style);
close(fig);
end

function files = local_plot_transform_audit(outRoot, style)
T = readtable(fullfile(outRoot, 'chunk_feature_transform_audit.csv'), 'TextType', 'string');
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1450 760]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
G = groupsummary(T, 'ChannelType');
bar(categorical(local_clean_label(G.ChannelType)), G.GroupCount, 'FaceColor', style.blue, 'EdgeColor', 'none');
ylabel('Channels', 'Interpreter', 'none');
title('Transformed channel types', 'Interpreter', 'none');
local_format_axis(gca, style);

nexttile;
plot(T.channel_index, T.robust_scale, 'o', 'Color', style.purple, 'MarkerFaceColor', style.purple);
xlabel('Channel index', 'Interpreter', 'none');
ylabel('Robust scale', 'Interpreter', 'none');
title('Condition-blind robust scaling factors', 'Interpreter', 'none');
local_format_axis(gca, style);

files = local_export(fig, fullfile(outRoot, 'figures'), 'feature_channel_transform_audit', style);
close(fig);
end

function files = local_plot_arena_shift(sourceDir, style)
T = readtable(fullfile(sourceDir, 'chunk_arena_shift_after_transform.csv'), 'TextType', 'string');
if isempty(T)
    files = strings(0, 1);
    return
end
[~, ord] = sort(abs(T.robust_median_difference_iqr), 'ascend');
T = T(ord, :);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1300 940]);
barh(T.robust_median_difference_iqr, 'FaceColor', style.gold, 'EdgeColor', 'none');
set(gca, 'YTick', 1:height(T), 'YTickLabel', local_clean_label(T.obs_name), ...
    'TickLabelInterpreter', 'none');
xlabel('Arena median shift after transform/scaling (IQR units)', 'Interpreter', 'none');
title('Arena shift diagnostic after condition-blind transform', 'Interpreter', 'none');
local_format_axis(gca, style);
files = local_export(fig, fullfile(fileparts(sourceDir), 'figures'), 'arena_shift_after_chunk_transform_audit_only', style);
close(fig);
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

function local_expand_count_axis(ax, values)
values = values(isfinite(values));
if isempty(values)
    return
end
upper = max(values);
if upper <= 0
    upper = 1;
end
if strcmpi(ax.YScale, 'log')
    positive = values(values > 0);
    if isempty(positive)
        ylim(ax, [1, upper * 1.08]);
    else
        lower = max(min(positive) * 0.8, eps);
        ylim(ax, [lower, upper * 1.2]);
    end
else
    ylim(ax, [0, upper * 1.08]);
end
end

function local_format_axis(ax, style)
set(ax, 'Box', 'off', 'FontName', style.fontName, 'FontSize', style.fontSize, ...
    'LineWidth', 0.8, 'TickDir', 'out');
grid(ax, 'on');
ax.GridColor = style.grid;
ax.GridAlpha = 0.35;
end
