function Fig = plot_multiscale_chunk_bank_diagnostics(ChunkSet, Report, varargin)
%PLOT_MULTISCALE_CHUNK_BANK_DIAGNOSTICS Journal-style scale-bank diagnostics.

p = inputParser;
p.addParameter('sessionIndex', 1, @(x)isscalar(x) && x >= 1);
p.addParameter('maxScalesForExamples', 9, @(x)isscalar(x) && x >= 1);
p.addParameter('exampleFeatures', ["centroid_dist","mutual_facing","radial_speed_12","in_contact"], @(x)isstring(x) || iscell(x));
p.addParameter('nExamplesPerScale', 5, @(x)isscalar(x) && x >= 1);
p.parse(varargin{:});
P = p.Results;
featNames = string(P.exampleFeatures);

Fig = struct();
R = Report.scaleSummary;

Fig.overview = figure('Color','w','Name','Dense scale-bank overview','Position',[60 60 1700 980]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile;
plot(R.chunk_sec, R.n_chunks, 'o-', 'LineWidth', 1.7); set(gca,'XScale','log');
xlabel('Chunk scale (s)'); ylabel('Number of chunks'); title('Coverage across log-spaced scales'); box off; grid on

nexttile;
plot(R.chunk_sec, R.mean_valid_frac, 'o-', 'LineWidth', 1.7); hold on
plot(R.chunk_sec, R.min_valid_frac, 's--', 'LineWidth', 1.2);
set(gca,'XScale','log'); ylim([0.85 1.01]);
xlabel('Chunk scale (s)'); ylabel('Valid fraction'); title('Mean and minimum validity by scale'); legend({'Mean','Minimum'}, 'Location','southwest'); box off; grid on

nexttile;
G = groupsummary(R, 'temporal_band', 'mean', 'chunk_sec');
G = sortrows(G, 'mean_chunk_sec');
bar(categorical(cellstr(string(G.temporal_band))), G.mean_chunk_sec);
ylabel('Mean scale (s)'); title('Temporal-band organization'); box off

nexttile;
scatter(log10(R.chunk_sec), R.chunk_density_per_sec, 55, double(categorical(R.temporal_band)), 'filled');
xlabel('log10 scale (s)'); ylabel('Chunk density per second'); title('Sampling density vs scale'); box off; grid on

nexttile;
plot(log10(R.chunk_sec(2:end)), R.log10_scale_gap(2:end), 'o-', 'LineWidth', 1.5);
xlabel('log10 scale (s)'); ylabel('Adjacent log-scale gap'); title('Scale spacing regularity'); box off; grid on

nexttile;
axis off
text(0.02, 0.95, { ...
    sprintf('Scales: %d', height(R)), ...
    sprintf('Range: %.3f to %.3f s', min(R.chunk_sec), max(R.chunk_sec)), ...
    sprintf('Median stride: %.1f frames', median(R.median_anchor_stride_frames, 'omitnan')), ...
    sprintf('Min validity observed: %.3f', min(R.min_valid_frac)), ...
    sprintf('Mean chunks / scale: %.1f', mean(R.n_chunks)), ...
    sprintf('Validation passed: %d', Report.isValid)}, 'VerticalAlignment','top','FontSize',12);
sgtitle('Dense log-spaced scale bank diagnostics','FontWeight','bold','FontSize',16);

% Example traces at representative scales.
idxShow = unique(round(linspace(1, numel(ChunkSet.scale), min(P.maxScalesForExamples, numel(ChunkSet.scale)))));
Fig.examples = figure('Color','w','Name','Representative scales: example chunk traces','Position',[90 90 1800 1100]);
tiledlayout(numel(idxShow), min(4,numel(featNames)), 'TileSpacing','compact','Padding','compact');
rng(1);
for rr = 1:numel(idxShow)
    s = idxShow(rr); Sc = ChunkSet.scale(s); meta = Sc.meta;
    idxSess = find(meta.session_index == P.sessionIndex);
    if isempty(idxSess)
        useIdx = [];
    else
        useIdx = idxSess(randperm(numel(idxSess), min(P.nExamplesPerScale, numel(idxSess))));
    end
    for c = 1:min(4, numel(featNames))
        nexttile; hold on
        fIdx = find(strcmp(ChunkSet.featureNames, featNames(c)), 1);
        if isempty(fIdx)
            axis off; continue
        end
        for jj = 1:numel(useIdx)
            x = squeeze(Sc.Xraw(useIdx(jj), :, fIdx));
            v = Sc.valid(useIdx(jj), :);
            x(~v) = NaN;
            tt = (0:numel(x)-1) ./ ChunkSet.sessions{P.sessionIndex}.fps;
            plot(tt, x, 'LineWidth', 0.8);
        end
        if rr == 1, title(strrep(featNames(c), '_', '\_'),'Interpreter','tex'); end
        if c == 1, ylabel(sprintf('%.2fs', Sc.chunkSec)); end
        if rr == numel(idxShow), xlabel('Time within chunk (s)'); end
        box off; grid on
    end
end

% Coverage ribbons for representative scales.
Fig.coverage = figure('Color','w','Name','Representative scales: session coverage','Position',[120 120 1700 1000]);
seq = ChunkSet.sessions{P.sessionIndex};
tiledlayout(numel(idxShow)+1,1,'TileSpacing','compact','Padding','compact');
nexttile; plot(seq.time, double(seq.validMask), 'k-'); ylim([-0.1 1.1]); ylabel('Valid'); title(sprintf('Session %d valid-mask timeline', P.sessionIndex)); box off
for rr = 1:numel(idxShow)
    nexttile; hold on; Sc = ChunkSet.scale(idxShow(rr)); meta = Sc.meta; idx = find(meta.session_index == P.sessionIndex);
    idx = idx(1:min(numel(idx), 600));
    for k = 1:numel(idx)
        plot([meta.start_time_s(idx(k)), meta.stop_time_s(idx(k))], [k k], '-', 'LineWidth', 0.8);
    end
    ylabel(sprintf('%.2fs', Sc.chunkSec)); box off
end
xlabel('Time (s)');
end
