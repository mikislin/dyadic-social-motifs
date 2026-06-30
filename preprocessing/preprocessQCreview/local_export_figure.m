function files = local_export_figure(fig, figDir, baseName, cfg)
files = strings(0,1);
figureSubdir = char(cfg.paths.preprocess_qc_figure_subdir);
if cfg.figures.export_png
    pngFile = fullfile(figDir, [baseName '.png']);
    exportgraphics(fig, pngFile, 'Resolution', cfg.figures.dpi);
    files(end+1,1) = replace(string(fullfile(figureSubdir, [baseName '.png'])), '\', '/');
end
if cfg.figures.export_pdf
    pdfFile = fullfile(figDir, [baseName '.pdf']);
    exportgraphics(fig, pdfFile, 'ContentType','vector');
    files(end+1,1) = replace(string(fullfile(figureSubdir, [baseName '.pdf'])), '\', '/');
end
end