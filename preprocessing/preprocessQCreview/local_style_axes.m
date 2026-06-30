function local_style_axes(ax, cfg)
set(ax, 'FontName', char(cfg.figures.font_name), 'FontSize', cfg.figures.font_size, ...
    'LineWidth', 0.8, 'Box','off', 'TickDir','out');
grid(ax, 'on');
ax.GridColor = local_cfg_rgb(cfg, 'grid_color');
ax.GridAlpha = 0.35;
title(ax, ax.Title.String, 'FontSize', cfg.figures.title_font_size, ...
    'FontWeight','bold');
end