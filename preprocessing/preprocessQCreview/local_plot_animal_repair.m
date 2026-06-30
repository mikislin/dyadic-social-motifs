function files = local_plot_animal_repair(animalTable, styles, cfg, figDir)
if isempty(animalTable) || height(animalTable) == 0
    files = strings(0,1);
    return
end

ok = isfinite(animalTable.pctPredictionIssueRepaired);
if ~any(ok)
    files = strings(0,1);
    return
end

G = groupsummary(animalTable(ok,:), {'analysis_group_id','plot_group_label','condition_id'}, ...
    'median', 'pctPredictionIssueRepaired');
G = local_sort_group_summary(G, styles);
[~, loc] = ismember(string(animalTable.plot_group_label), string(G.plot_group_label));
valid = ok & loc > 0;

fig = figure('Visible','off', 'Color','w', 'Position',[80 80 1200 1000]);
hold on;
for i = 1:height(G)
    idx = valid & loc == i;
    x = animalTable.pctPredictionIssueRepaired(idx);
    y = i + local_jitter(nnz(idx), cfg);
    c = local_condition_rgb(styles, G.condition_id(i));
    scatter(x, y, cfg.figures.marker_size, c, 'filled', ...
        'MarkerFaceAlpha', 0.65, 'MarkerEdgeColor', 'none');
    plot([G.median_pctPredictionIssueRepaired(i) G.median_pctPredictionIssueRepaired(i)], ...
        [i-0.32 i+0.32], 'Color', c, 'LineWidth', 2.2);
end
set(gca, 'YTick', 1:height(G), 'YTickLabel', G.plot_group_label);
xlabel('% prediction-issue animal-frames repaired and usable');
ylabel('Analysis group (cohort and arena preserved)');
title('Animal-level repair audit');
xMax = max(animalTable.pctPredictionIssueRepaired(valid));
if isfinite(xMax)
    xlim([0 min(100, max(55, ceil(xMax * 1.1 / 5) * 5))]);
else
    xlim([0 100]);
end
local_style_axes(gca, cfg);

files = local_export_figure(fig, figDir, 'preprocess_qc_animal_repair_by_group', cfg);
close(fig);
end
