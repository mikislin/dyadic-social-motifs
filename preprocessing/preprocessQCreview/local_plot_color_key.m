function files = local_plot_color_key(styles, cfg, figDir)
fig = figure('Visible','off', 'Color','w', 'Position',[80 80 800 850]);
hold on;
for i = 1:height(styles.table)
    barh(i, 1, 'FaceColor', styles.colors(i,:), 'EdgeColor', 'none');
end
set(gca, 'YTick', 1:height(styles.table), 'YTickLabel', styles.table.display_label, ...
    'XTick', []);
xlim([0 1]);
ylabel('Experimental group');
title('Stable condition color key');
local_style_axes(gca, cfg);
files = local_export_figure(fig, figDir, 'preprocess_qc_experimental_group_color_key', cfg);
close(fig);
end