function y = local_jitter(n, cfg)
if n == 0
    y = [];
else
    y = linspace(-cfg.figures.jitter_width, cfg.figures.jitter_width, n)';
end
end