function qc = review_preprocessing_outputs(preprocDir, outDir)
%REVIEW_PREPROCESSING_OUTPUTS Review preprocessing QC across session output files.
%
% Usage:
%   qc = review_preprocessing_outputs( ...
%       '/projects/m/mkislin/DyadicInteractions/preproc_outputs', ...
%       '/projects/m/mkislin/DyadicInteractions/preproc_qc_review');
%
% Inputs
%   preprocDir : folder with session_*_preproc.mat files
%   outDir     : folder to save QC tables and figures
%
% Output
%   qc : struct with fields
%       .summaryTable
%       .animalTable
%       .files
%
% Notes
%   - expects each MAT file to contain variable `out` or `sessionPreproc`
%   - works with 1-animal and 2-animal sessions
%   - saves a CSV summary and several PNG/PDF figures

if nargin < 1 || isempty(preprocDir)
    error('Please provide preprocDir');
end
if nargin < 2 || isempty(outDir)
    outDir = fullfile(preprocDir, 'qc_review');
end
if ~isfolder(outDir)
    mkdir(outDir);
end

files = dir(fullfile(preprocDir, '*_preproc.mat'));
assert(~isempty(files), 'No *_preproc.mat files found in %s', preprocDir);

[~, order] = sort({files.name});
files = files(order);

sessionRows = [];
animalRows = [];

fprintf('Reviewing %d preprocessing output files...\n', numel(files));

for i = 1:numel(files)
    inFile = fullfile(files(i).folder, files(i).name);
    S = load(inFile);

    if isfield(S, 'out')
        P = S.out;
    elseif isfield(S, 'sessionPreproc')
        P = S.sessionPreproc;
    else
        warning('Skipping %s: no `out` or `sessionPreproc` variable found.', files(i).name);
        continue
    end

    nAnimals = size(P.clean.tracks, 4);
    nFrames = size(P.clean.tracks, 1);
    nNodes = size(P.clean.tracks, 2);

    sessionId = erase(files(i).name, '_preproc.mat');

    % Session-level aggregates
    badframes = P.qc.badframes;
    if isvector(badframes)
        badframes = badframes(:);
    end

    sessionPctBad_any = mean(any(logical(badframes), 2), 'omitnan') * 100;
    sessionPctBad_all = mean(all(logical(badframes), 2), 'omitnan') * 100;

    row = table;
    row.session_id = string(sessionId);
    row.file_name = string(files(i).name);
    row.nFrames = nFrames;
    row.nNodes = nNodes;
    row.nAnimals = nAnimals;
    row.pctBad_anyAnimal = sessionPctBad_any;
    row.pctBad_allAnimals = sessionPctBad_all;

    if isfield(P, 'debug') && isfield(P.debug, 'frameHasAnyNaN')
        row.pctFrameHasAnyNaN_raw = mean(logical(P.debug.frameHasAnyNaN)) * 100;
    else
        row.pctFrameHasAnyNaN_raw = NaN;
    end

    sessionRows = [sessionRows; row]; %#ok<AGROW>

    % Animal-level stats
    if isfield(P.qc, 'sessionStats') && isfield(P.qc.sessionStats, 'animal')
        statsAnimal = P.qc.sessionStats.animal;
    else
        statsAnimal = struct([]);
    end

    for m = 1:nAnimals
        arow = table;
        arow.session_id = string(sessionId);
        arow.file_name = string(files(i).name);
        arow.animal = m;
        arow.nFrames = nFrames;

        if ~isempty(statsAnimal) && numel(statsAnimal) >= m
            st = statsAnimal(m);

            arow.pctBadframes = getfield_safe(st, 'pctBadframes');
            arow.pctInterpFrames = getfield_safe(st, 'pctInterpFrames');
            arow.pctJumpFrames = getfield_safe(st, 'pctJumpFrames');
            arow.pctGeomFrames = getfield_safe(st, 'pctGeomFrames');
            arow.pctArenaFrames = getfield_safe(st, 'pctArenaFrames');
            arow.pctLowConfFrames = getfield_safe(st, 'pctLowConfFrames');
            arow.pctJumpSamples = getfield_safe(st, 'pctJumpSamples');
            arow.pctInterpSamples = getfield_safe(st, 'pctInterpSamples');
            arow.medianBodyLength = getfield_safe(st, 'medianBodyLength');
            arow.medianBodyLengthRaw = getfield_safe(st, 'medianBodyLengthRaw');
            arow.medianDistortionAnchorDispRatio = getfield_safe(st, 'medianDistortionAnchorDispRatio');
            arow.medianCentroidSpeedRaw = getfield_safe(st, 'medianCentroidSpeedRaw');
            arow.medianCentroidSpeedClean = getfield_safe(st, 'medianCentroidSpeedClean');
        else
            arow.pctBadframes = NaN;
            arow.pctInterpFrames = NaN;
            arow.pctJumpFrames = NaN;
            arow.pctGeomFrames = NaN;
            arow.pctArenaFrames = NaN;
            arow.pctLowConfFrames = NaN;
            arow.pctJumpSamples = NaN;
            arow.pctInterpSamples = NaN;
            arow.medianBodyLength = NaN;
            arow.medianBodyLengthRaw = NaN;
            arow.medianDistortionAnchorDispRatio = NaN;
            arow.medianCentroidSpeedRaw = NaN;
            arow.medianCentroidSpeedClean = NaN;
        end

        % Optional direct mask-derived counts for cross-checking
        if isfield(P.qc, 'animals') && numel(P.qc.animals) >= m
            qa = P.qc.animals(m);

            if isfield(qa, 'jumpMask')
                jm = full(qa.jumpMask);
                arow.maskPctJumpFrames_anyNode = mean(any(jm > 0, 2)) * 100;
            else
                arow.maskPctJumpFrames_anyNode = NaN;
            end

            if isfield(qa, 'interpMask')
                im = full(qa.interpMask);
                arow.maskPctInterpFrames_anyNode = mean(any(im > 0, 2)) * 100;
            else
                arow.maskPctInterpFrames_anyNode = NaN;
            end

            if isfield(qa, 'finalNanMask')
                fnm = full(qa.finalNanMask);
                arow.maskPctFinalNaNFrames_anyNode = mean(any(fnm > 0, 2)) * 100;
            else
                arow.maskPctFinalNaNFrames_anyNode = NaN;
            end
        else
            arow.maskPctJumpFrames_anyNode = NaN;
            arow.maskPctInterpFrames_anyNode = NaN;
            arow.maskPctFinalNaNFrames_anyNode = NaN;
        end

        animalRows = [animalRows; arow]; %#ok<AGROW>
    end
end

assert(~isempty(sessionRows), 'No valid preprocessing outputs were parsed.');

summaryTable = sessionRows;
animalTable = animalRows;

% Save tables
writetable(summaryTable, fullfile(outDir, 'preprocessing_qc_session_summary.csv'));
writetable(animalTable, fullfile(outDir, 'preprocessing_qc_animal_summary.csv'));

% Save MAT bundle
qc = struct();
qc.summaryTable = summaryTable;
qc.animalTable = animalTable;
qc.files = files;
save(fullfile(outDir, 'preprocessing_qc_review.mat'), 'qc', '-v7.3');

% Print quick text summary
print_summary_to_console(summaryTable, animalTable);

% Make figures
make_overview_histograms(animalTable, outDir);
make_ranked_session_plot(animalTable, outDir);
make_metric_heatmap(animalTable, outDir);
make_scatter_diagnostics(animalTable, outDir);
make_flag_report(animalTable, outDir);

fprintf('QC review saved to: %s\n', outDir);
end


function v = getfield_safe(s, f)
if isstruct(s) && isfield(s, f)
    v = s.(f);
else
    v = NaN;
end
end


function print_summary_to_console(summaryTable, animalTable)
fprintf('\n===== Preprocessing QC summary =====\n');
fprintf('Sessions reviewed: %d\n', height(summaryTable));
fprintf('Animal rows:       %d\n\n', height(animalTable));

metrics = { ...
    'pctBadframes', ...
    'pctInterpFrames', ...
    'pctJumpFrames', ...
    'pctGeomFrames', ...
    'pctArenaFrames', ...
    'pctLowConfFrames', ...
    'medianDistortionAnchorDispRatio'};

for i = 1:numel(metrics)
    x = animalTable.(metrics{i});
    x = x(isfinite(x));
    if isempty(x)
        continue
    end
    fprintf('%-32s median=%7.3f   p90=%7.3f   max=%7.3f\n', ...
        metrics{i}, median(x), prctile(x,90), max(x));
end
fprintf('\n');
end


function make_overview_histograms(T, outDir)
metrics = { ...
    'pctBadframes', ...
    'pctInterpFrames', ...
    'pctJumpFrames', ...
    'pctGeomFrames', ...
    'pctArenaFrames', ...
    'pctLowConfFrames', ...
    'pctJumpSamples', ...
    'pctInterpSamples', ...
    'medianBodyLength', ...
    'medianDistortionAnchorDispRatio'};

fig = figure('Visible','off','Position',[100 100 1600 900]);
tl = tiledlayout(fig, 2, 5, 'TileSpacing','compact', 'Padding','compact');

for i = 1:numel(metrics)
    nexttile(tl);
    x = T.(metrics{i});
    x = x(isfinite(x));
    if isempty(x)
        title(metrics{i}, 'Interpreter','none');
        axis off;
        continue
    end
    histogram(x, 25);
    title(strrep(metrics{i}, '_', '\_'), 'Interpreter','tex');
    box off;
end

exportgraphics(fig, fullfile(outDir, 'qc_histograms.png'), 'Resolution', 220);
exportgraphics(fig, fullfile(outDir, 'qc_histograms.pdf'));
close(fig);
end


function make_ranked_session_plot(T, outDir)
% Rank sessions by badframes and overlay main error sources
[groups, sessionNames] = findgroups(T.session_id);

pctBad = splitapply(@median, T.pctBadframes, groups);
pctInterp = splitapply(@median, T.pctInterpFrames, groups);
pctJump = splitapply(@median, T.pctJumpFrames, groups);
pctGeom = splitapply(@median, T.pctGeomFrames, groups);
pctArena = splitapply(@median, T.pctArenaFrames, groups);
pctLowConf = splitapply(@median, T.pctLowConfFrames, groups);

M = table(sessionNames, pctBad, pctInterp, pctJump, pctGeom, pctArena, pctLowConf);
M = sortrows(M, 'pctBad', 'descend');

fig = figure('Visible','off','Position',[100 100 1700 800]);
tiledlayout(fig, 2, 1, 'TileSpacing','compact', 'Padding','compact');

nexttile;
bar(M.pctBad);
ylabel('% bad frames');
title('Sessions ranked by median % bad frames across animals');
xlim([0.5 height(M)+0.5]);
box off;

nexttile;
plot(M.pctInterp, '-','LineWidth',1); hold on;
plot(M.pctJump, '-','LineWidth',1);
plot(M.pctGeom, '-','LineWidth',1);
plot(M.pctArena, '-','LineWidth',1);
plot(M.pctLowConf, '-','LineWidth',1);
ylabel('% frames');
xlabel('Session rank');
legend({'interp','jump','geom','arena','lowconf'}, 'Location','northeastoutside');
title('Median QC components across ranked sessions');
xlim([1 height(M)]);
box off;

exportgraphics(fig, fullfile(outDir, 'qc_ranked_sessions.png'), 'Resolution', 220);
exportgraphics(fig, fullfile(outDir, 'qc_ranked_sessions.pdf'));
close(fig);

writetable(M, fullfile(outDir, 'qc_ranked_sessions.csv'));
end


function make_metric_heatmap(T, outDir)
[groups, sessionNames] = findgroups(T.session_id);

metrics = { ...
    'pctBadframes', ...
    'pctInterpFrames', ...
    'pctJumpFrames', ...
    'pctGeomFrames', ...
    'pctArenaFrames', ...
    'pctLowConfFrames', ...
    'pctJumpSamples', ...
    'pctInterpSamples', ...
    'medianDistortionAnchorDispRatio'};

X = nan(numel(sessionNames), numel(metrics));
for j = 1:numel(metrics)
    X(:,j) = splitapply(@median, T.(metrics{j}), groups);
end

% z-score each column for comparability
Xz = X;
for j = 1:size(Xz,2)
    col = Xz(:,j);
    mu = mean(col, 'omitnan');
    sd = std(col, 'omitnan');
    if isfinite(sd) && sd > 0
        Xz(:,j) = (col - mu) ./ sd;
    else
        Xz(:,j) = 0;
    end
end

fig = figure('Visible','off','Position',[100 100 1500 900]);
imagesc(Xz);
colorbar;
colormap(parula);
xlabel('Metric');
ylabel('Session');
title('Session-by-metric QC heatmap (column z-score)');
set(gca, 'XTick', 1:numel(metrics), 'XTickLabel', metrics, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:numel(sessionNames), 'YTickLabel', sessionNames);

exportgraphics(fig, fullfile(outDir, 'qc_heatmap.png'), 'Resolution', 220);
exportgraphics(fig, fullfile(outDir, 'qc_heatmap.pdf'));
close(fig);
end


function make_scatter_diagnostics(T, outDir)
fig = figure('Visible','off','Position',[100 100 1400 1000]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');

nexttile(tl);
scatter(T.pctLowConfFrames, T.pctBadframes, 30, 'filled');
xlabel('% low-conf frames');
ylabel('% bad frames');
title('Badframes vs low confidence');
box off;

nexttile(tl);
scatter(T.pctInterpFrames, T.pctBadframes, 30, 'filled');
xlabel('% interp frames');
ylabel('% bad frames');
title('Badframes vs interpolation');
box off;

nexttile(tl);
scatter(T.pctJumpFrames, T.pctBadframes, 30, 'filled');
xlabel('% jump frames');
ylabel('% bad frames');
title('Badframes vs jump corrections');
box off;

nexttile(tl);
scatter(T.medianDistortionAnchorDispRatio, T.pctBadframes, 30, 'filled');
xlabel('Median anchor distortion ratio');
ylabel('% bad frames');
title('Badframes vs distortion');
box off;

exportgraphics(fig, fullfile(outDir, 'qc_scatter_diagnostics.png'), 'Resolution', 220);
exportgraphics(fig, fullfile(outDir, 'qc_scatter_diagnostics.pdf'));
close(fig);
end


function make_flag_report(T, outDir)
% Flag suspicious sessions for manual review
[groups, sessionNames] = findgroups(T.session_id);

pctBad = splitapply(@median, T.pctBadframes, groups);
pctInterp = splitapply(@median, T.pctInterpFrames, groups);
pctJump = splitapply(@median, T.pctJumpFrames, groups);
pctGeom = splitapply(@median, T.pctGeomFrames, groups);
pctArena = splitapply(@median, T.pctArenaFrames, groups);
pctLowConf = splitapply(@median, T.pctLowConfFrames, groups);
distort = splitapply(@median, T.medianDistortionAnchorDispRatio, groups);

flagTbl = table(sessionNames, pctBad, pctInterp, pctJump, pctGeom, pctArena, pctLowConf, distort);

flagTbl.flag_high_bad = flagTbl.pctBad > 25;
flagTbl.flag_high_interp = flagTbl.pctInterp > 20;
flagTbl.flag_high_jump = flagTbl.pctJump > 20;
flagTbl.flag_high_geom = flagTbl.pctGeom > 15;
flagTbl.flag_high_arena = flagTbl.pctArena > 10;
flagTbl.flag_high_lowconf = flagTbl.pctLowConf > 15;
flagTbl.flag_high_distortion = flagTbl.distort > 0.10;

flagTbl.nFlags = ...
    double(flagTbl.flag_high_bad) + ...
    double(flagTbl.flag_high_interp) + ...
    double(flagTbl.flag_high_jump) + ...
    double(flagTbl.flag_high_geom) + ...
    double(flagTbl.flag_high_arena) + ...
    double(flagTbl.flag_high_lowconf) + ...
    double(flagTbl.flag_high_distortion);

flagTbl = sortrows(flagTbl, {'nFlags','pctBad'}, {'descend','descend'});

writetable(flagTbl, fullfile(outDir, 'qc_flag_report.csv'));

fig = figure('Visible','off','Position',[100 100 1200 500]);
bar(flagTbl.nFlags);
xlabel('Ranked session');
ylabel('# triggered flags');
title('Manual review priority by QC flags');
box off;

exportgraphics(fig, fullfile(outDir, 'qc_flag_priority.png'), 'Resolution', 220);
exportgraphics(fig, fullfile(outDir, 'qc_flag_priority.pdf'));
close(fig);
end