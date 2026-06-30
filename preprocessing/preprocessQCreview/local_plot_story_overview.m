function files = local_plot_story_overview(T, styles, cfg, figDir)
ok = isfinite(T.badframe_fraction);
if ~any(ok)
    files = strings(0,1);
    return
end

fig = figure('Visible','off', 'Color','w', 'Position',[80 80 1500 950]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');

nexttile(tl);
metricNames = ["prediction_issue_fraction","repaired_prediction_issue_fraction", ...
    "unresolved_prediction_issue_fraction","badframe_fraction"];
metricLabels = ["Prediction issue","Repaired","Unresolved","Final bad"];
stageColors = [ ...
    local_cfg_rgb(cfg, 'story_color_prediction_issue'); ...
    local_cfg_rgb(cfg, 'story_color_repaired'); ...
    local_cfg_rgb(cfg, 'story_color_unresolved'); ...
    local_cfg_rgb(cfg, 'story_color_badframe')];
values = nan(1, numel(metricNames));
lo = values;
hi = values;
for i = 1:numel(metricNames)
    x = T.(metricNames(i));
    x = x(isfinite(x)) .* 100;
    values(i) = median(x);
    lo(i) = prctile(x, 25);
    hi(i) = prctile(x, 75);
end
b = bar(values, 'FaceColor','flat');
b.CData = stageColors;
hold on;
errorbar(1:numel(values), values, values-lo, hi-values, 'k.', 'LineWidth', 0.9);
set(gca, 'XTick', 1:numel(metricLabels), 'XTickLabel', metricLabels);
ylabel('% animal-frames');
title('A. What preprocessing encountered and retained');
local_style_axes(gca, cfg);

nexttile(tl);
hold on;
conds = unique(string(T.condition_id(ok)), 'stable');
for i = 1:numel(conds)
    idx = ok & string(T.condition_id) == conds(i);
    scatter(T.prediction_issue_fraction(idx).*100, T.badframe_fraction(idx).*100, ...
        cfg.figures.marker_size, local_condition_rgb(styles, conds(i)), 'filled');
end
plot([0 100], [0 100], '--', 'Color', local_cfg_rgb(cfg, 'grid_color'), 'LineWidth', 1);
xlabel('% animal-frames with prediction issue');
ylabel('% final bad animal-frames');
title('B. Flagged prediction issues are audited separately from final exclusions');
xlim([0 max(65, max(T.prediction_issue_fraction(ok).*100) * 1.05)]);
ylim([0 max(60, max(T.badframe_fraction(ok).*100) * 1.05)]);
local_style_axes(gca, cfg);

nexttile(tl);
repair = T.prediction_issue_repair_rate(isfinite(T.prediction_issue_repair_rate));
usable = T.prediction_issue_usable_rate(isfinite(T.prediction_issue_usable_rate));
unresolved = T.prediction_issue_unresolved_rate(isfinite(T.prediction_issue_unresolved_rate));
usableNoInterp = max(median(usable) - median(repair), 0);
stackVals = 100 .* [median(repair), usableNoInterp, median(unresolved)];
b = bar(1, stackVals, 'stacked');
b(1).FaceColor = local_cfg_rgb(cfg, 'story_color_repaired');
b(2).FaceColor = local_cfg_rgb(cfg, 'story_color_interpolated');
b(3).FaceColor = local_cfg_rgb(cfg, 'story_color_unresolved');
set(gca, 'XTick', 1, 'XTickLabel', "Prediction-issue frames");
ylabel('% of prediction-issue animal-frames');
legend({'repaired and usable','usable without interpolation','unresolved'}, ...
    'Location','southoutside', 'Orientation','horizontal', 'Box','off');
title('C. Fate of prediction-issue frames');
ylim([0 100]);
local_style_axes(gca, cfg);

nexttile(tl);
motifMask = select_analysis_cohort(T, "motif_discovery");
reviewMask = T.qc_warn_motif_discovery;
passCleanMask = T.qc_pass_motif_discovery & ~reviewMask;
failMask = motifMask & T.preprocess_success & ~T.qc_pass_motif_discovery;
counts = [nnz(passCleanMask), nnz(reviewMask), nnz(failMask)];
b = bar(counts, 'FaceColor','flat');
b.CData = [local_cfg_rgb(cfg, 'story_color_repaired'); ...
    local_cfg_rgb(cfg, 'story_color_interpolated'); ...
    local_cfg_rgb(cfg, 'story_color_badframe')];
set(gca, 'XTick', 1:3, 'XTickLabel', ["pass","review","fail"]);
text(1:3, counts + max(counts) * 0.03, string(counts), ...
    'HorizontalAlignment','center', 'FontName', char(cfg.figures.font_name), ...
    'FontSize', cfg.figures.font_size);
ylabel('Motif-discovery sessions');
title('D. Motif-discovery QC decision');
ylim([0 max(counts) * 1.15]);
local_style_axes(gca, cfg);

files = local_export_figure(fig, figDir, 'preprocess_qc_story_overview', cfg);
close(fig);
end
