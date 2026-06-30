function files = local_plot_condition_arena_summary(T, styles, cfg, figDir)
ok = T.preprocess_success & isfinite(T.badframe_fraction);
if ~any(ok)
    files = strings(0,1);
    return
end

G = groupsummary(T(ok,:), {'arena_condition','arena_condition_label','condition_id'}, ...
    'median', {'badframe_fraction','prediction_issue_repair_rate'});
G = local_sort_group_summary(G, styles);

fig = figure('Visible','off', 'Color','w', 'Position',[80 80 1250 1000]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing','compact', 'Padding','compact');

nexttile(tl);
hold on;
for i = 1:height(G)
    barh(i, G.median_badframe_fraction(i) * 100, ...
        'FaceColor', local_condition_rgb(styles, G.condition_id(i)), 'EdgeColor', 'none');
end
set(gca, 'YTick', 1:height(G), 'YTickLabel', G.arena_condition_label);
xlabel('Median % bad animal-frames');
title('Final exclusion burden');
local_style_axes(gca, cfg);

nexttile(tl);
hold on;
for i = 1:height(G)
    barh(i, G.median_prediction_issue_repair_rate(i) * 100, ...
        'FaceColor', local_condition_rgb(styles, G.condition_id(i)), 'EdgeColor', 'none');
end
set(gca, 'YTick', 1:height(G), 'YTickLabel', []);
xlabel('Median % prediction issues repaired');
title('Repair yield');
xlim([0 100]);
local_style_axes(gca, cfg);

files = local_export_figure(fig, figDir, 'preprocess_qc_condition_arena_summary', cfg);
close(fig);
end
