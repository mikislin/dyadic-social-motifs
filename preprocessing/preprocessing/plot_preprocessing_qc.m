function fig = plot_preprocessing_qc(sessionPreproc, outFile)
%PLOT_PREPROCESSING_QC Create a single high-resolution audit figure.

if nargin < 2
    outFile = '';
end

tracksClean = sessionPreproc.clean.tracks;
if ndims(tracksClean) == 3
    tracksClean = reshape(tracksClean, size(tracksClean,1), size(tracksClean,2), size(tracksClean,3), 1);
end
nAnimals = size(tracksClean,4);
anchorNode = min(sessionPreproc.params.qc.arena_anchor_node, size(tracksClean,2));

hasRaw = ~isempty(sessionPreproc.raw) && isfield(sessionPreproc.raw, 'SLEAPtracks') && ~isempty(sessionPreproc.raw.SLEAPtracks);
if hasRaw
    tracksRaw = sessionPreproc.raw.SLEAPtracks;
    if ndims(tracksRaw) == 3
        tracksRaw = reshape(tracksRaw, size(tracksRaw,1), size(tracksRaw,2), size(tracksRaw,3), 1);
    end
else
    tracksRaw = [];
end

fig = figure('Visible','off','Color','w','Position',[100 100 1600 900]);
tiledlayout(3, nAnimals + 1, 'TileSpacing','compact', 'Padding','compact');

for m = 1:nAnimals
    xyClean = squeeze(tracksClean(:,anchorNode,:,m));

    nexttile;
    if hasRaw
        xyRaw = squeeze(tracksRaw(:,anchorNode,:,m));
        plot(xyRaw(:,1), xyRaw(:,2), '-', 'LineWidth', 0.4); hold on;
        plot(xyClean(:,1), xyClean(:,2), '-', 'LineWidth', 0.6);
        legend({'raw','clean'}, 'Location','best');
    else
        plot(xyClean(:,1), xyClean(:,2), '-', 'LineWidth', 0.6);
        legend({'clean'}, 'Location','best');
    end
    axis equal tight;
    title(sprintf('Animal %d anchor trajectory', m));

    nexttile;
    dClean = sqrt(sum(diff(xyClean,1,1).^2,2));
    if hasRaw
        dRaw = sqrt(sum(diff(xyRaw,1,1).^2,2));
        histogram(dRaw, 120, 'Normalization','probability'); hold on;
        histogram(dClean, 120, 'Normalization','probability');
        legend({'raw','clean'}, 'Location','best');
    else
        histogram(dClean, 120, 'Normalization','probability');
        legend({'clean'}, 'Location','best');
    end
    title(sprintf('Animal %d anchor displacement', m));

    nexttile;
    bodyLen = sessionPreproc.qc.frames.bodyLength(:,m);
    yBad = sessionPreproc.qc.badframes(:,m) .* max(bodyLen,[],'omitnan');
    plot(bodyLen, 'k-', 'LineWidth', 0.6); hold on;
    plot(yBad, 'r-');
    title(sprintf('Animal %d body length + badframes', m));
end

nexttile;
imagesc([sessionPreproc.qc.frames.fracJump, ...
         sessionPreproc.qc.frames.fracInterp, ...
         sessionPreproc.qc.frames.fracGeom, ...
         sessionPreproc.qc.frames.fracArena]');
colormap(parula);
colorbar;
yticks(1:4*nAnimals);
yticklabels(make_labels(nAnimals));
title('Frame-level QC fractions');
xlabel('Frame');

sgtitle(sprintf('Preprocessing audit: %s', sessionPreproc.session_id), 'Interpreter','none');

if ~isempty(outFile)
    dpi = 220;
    if isfield(sessionPreproc.params, 'output') && isfield(sessionPreproc.params.output, 'audit_plot_dpi')
        dpi = sessionPreproc.params.output.audit_plot_dpi;
    end
    exportgraphics(fig, outFile, 'Resolution', dpi);
end
end

function labels = make_labels(nAnimals)
labels = cell(4*nAnimals,1);
k = 1;
for m = 1:nAnimals
    labels{k} = sprintf('jump a%d', m); k = k + 1;
    labels{k} = sprintf('interp a%d', m); k = k + 1;
    labels{k} = sprintf('geom a%d', m); k = k + 1;
    labels{k} = sprintf('arena a%d', m); k = k + 1;
end
end
