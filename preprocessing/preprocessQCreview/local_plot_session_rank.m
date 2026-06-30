function files = local_plot_session_rank(T, styles, cfg, figDir)
ok = T.preprocess_success & isfinite(T.badframe_fraction);
if ~any(ok)
    files = strings(0,1);
    return
end

R = sortrows(T(ok,:), 'badframe_fraction', 'descend');
fig = figure('Visible','off', 'Color','w', 'Position',[80 80 1500 850]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing','compact', 'Padding','compact');

nexttile(tl);
hold on;
for i = 1:height(R)
    c = local_condition_rgb(styles, R.condition_id(i));
    bar(i, R.badframe_fraction(i) * 100, 'FaceColor', c, 'EdgeColor', 'none');
end
yline(cfg.qc.motif_warn_badframe_fraction * 100, '--', 'review', ...
    'Color', local_cfg_rgb(cfg, 'story_color_interpolated'), 'LineWidth', 1);
yline(cfg.qc.motif_max_badframe_fraction * 100, '--', 'fail', ...
    'Color', local_cfg_rgb(cfg, 'story_color_badframe'), 'LineWidth', 1);
ylabel('% bad animal-frames');
title('Dyadic/session-level bad-frame ranking');
xlim([0.5 height(R)+0.5]);
local_style_axes(gca, cfg);

nexttile(tl);
hold on;
for i = 1:height(R)
    c = local_condition_rgb(styles, R.condition_id(i));
    scatter(i, R.prediction_issue_repair_rate(i) * 100, cfg.figures.marker_size, c, 'filled');
end
ylabel('% prediction issues repaired');
xlabel('Session rank by bad-frame fraction');
xlim([0.5 height(R)+0.5]);
ylim([0 100]);
local_style_axes(gca, cfg);

files = local_export_figure(fig, figDir, 'preprocess_qc_dyadic_session_rank', cfg);
close(fig);
end