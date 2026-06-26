%RUN_04_PREPROCESS_QC_REVIEW Summarize preprocessing QC and paper cohorts.
%
% Plot colors
% and thresholds are loaded from visible files in config/.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
cd(repoRoot);
addpath(genpath(repoRoot));

paths = paper_paths();
cfg = load_preprocessing_pipeline_config(paths.preprocessingPipelineConfigPath);
manifestPath = paths.manifestPath;
preprocRoot = local_repo_path(repoRoot, cfg.paths.preprocessed_output_dir);
summaryPath = fullfile(preprocRoot, char(cfg.paths.preprocess_summary_file));
outDir = local_repo_path(repoRoot, cfg.paths.preprocess_qc_output_dir);
figDir = fullfile(outDir, char(cfg.paths.preprocess_qc_figure_subdir));
if ~exist(outDir, 'dir'); mkdir(outDir); end
if ~exist(figDir, 'dir'); mkdir(figDir); end

if ~isfile(manifestPath)
    run(fullfile(repoRoot, 'paper', 'run_01_manifest.m'));
end
assert(isfile(summaryPath), 'run_04_preprocess_qc_review:MissingSummary', ...
    'Missing preprocessing summary: %s. Run paper/run_03_preprocess_all.m first.', summaryPath);

M = readtable(manifestPath, 'TextType', 'string');
S = readtable(summaryPath, 'TextType', 'string');
stylePath = fullfile(paths.configDir, char(cfg.qc.condition_style_file));
styles = load_experimental_group_styles(stylePath, M);

T = local_merge_manifest_preprocess_summary(M, S);
T.analysis_class = local_analysis_class(T);

motifWarnBadframeFraction = cfg.qc.motif_warn_badframe_fraction;
motifMaxBadframeFraction = cfg.qc.motif_max_badframe_fraction;
motifMask = select_analysis_cohort(T, "motif_discovery");
T.qc_pass_motif_discovery = motifMask & T.preprocess_success & ...
    T.badframe_fraction <= motifMaxBadframeFraction;
T.qc_warn_motif_discovery = motifMask & T.preprocess_success & ...
    T.badframe_fraction > motifWarnBadframeFraction & ...
    T.badframe_fraction <= motifMaxBadframeFraction;
T.qc_requires_manual_review = motifMask & T.preprocess_success & ...
    T.badframe_fraction > motifWarnBadframeFraction;
T.qc_notes = local_qc_notes(T, motifMask, motifWarnBadframeFraction, motifMaxBadframeFraction);

[animalTable, pointTable, detailAudit] = local_build_detail_tables(T, repoRoot, preprocRoot, cfg);

sessionTablePath = fullfile(outDir, 'preprocess_qc_session_table.csv');
animalTablePath = fullfile(outDir, 'preprocess_qc_animal_table.csv');
pointTablePath = fullfile(outDir, 'preprocess_qc_point_table.csv');
motifTablePath = fullfile(outDir, 'motif_discovery_qc_sessions.csv');
anesthTablePath = fullfile(outDir, 'anesthetized_context_sessions.csv');
conditionSummaryPath = fullfile(outDir, 'preprocess_qc_condition_summary.csv');
analysisSummaryPath = fullfile(outDir, 'preprocess_qc_analysis_class_summary.csv');
animalSummaryPath = fullfile(outDir, 'preprocess_qc_animal_group_summary.csv');
pointSummaryPath = fullfile(outDir, 'preprocess_qc_point_summary.csv');
figureManifestPath = fullfile(outDir, 'preprocess_qc_figure_manifest.csv');

writetable(T, sessionTablePath);
writetable(animalTable, animalTablePath);
writetable(pointTable, pointTablePath);
writetable(T(motifMask,:), motifTablePath);
writetable(T(select_analysis_cohort(T, "anesthetized_context"),:), anesthTablePath);

summaryMetrics = local_summary_metric_vars(T);
conditionSummary = groupsummary(T, {'analysis_class','arena_condition','condition'}, ...
    {'mean','median','max'}, summaryMetrics);
analysisSummary = groupsummary(T, 'analysis_class', ...
    {'mean','median','max'}, summaryMetrics);
writetable(conditionSummary, conditionSummaryPath);
writetable(analysisSummary, analysisSummaryPath);

animalSummary = local_animal_group_summary(animalTable);
pointSummary = local_point_summary(pointTable);
writetable(animalSummary, animalSummaryPath);
writetable(pointSummary, pointSummaryPath);

figureManifest = local_make_figures(T, animalTable, pointTable, pointSummary, styles, cfg, figDir);
writetable(figureManifest, figureManifestPath);

qc = struct();
qc.sessionTable = T;
qc.animalTable = animalTable;
qc.pointTable = pointTable;
qc.conditionSummary = conditionSummary;
qc.analysisSummary = analysisSummary;
qc.animalSummary = animalSummary;
qc.pointSummary = pointSummary;
qc.detailAudit = detailAudit;
qc.params = cfg;
qc.styles = styles.table;
save(fullfile(outDir, 'preprocess_qc_review.mat'), 'qc', '-v7.3');

fprintf('run_04_preprocess_qc_review\n');
fprintf('Rows: %d\n', height(T));
fprintf('Motif discovery candidates: %d\n', nnz(motifMask));
fprintf('Motif QC pass: %d\n', nnz(T.qc_pass_motif_discovery));
fprintf('Motif manual-review warnings: %d\n', nnz(T.qc_requires_manual_review));
fprintf('Animal QC rows: %d\n', height(animalTable));
fprintf('Point QC rows: %d\n', height(pointTable));
fprintf('Anesthetized-center Block 2 sessions: %d\n', ...
    nnz(select_analysis_cohort(T, "anesthetized_context")));
fprintf('Wrote session QC table: %s\n', sessionTablePath);
fprintf('Wrote figure manifest: %s\n', figureManifestPath);
disp(analysisSummary);

function T = local_merge_manifest_preprocess_summary(M, S)
T = M;
T.preprocess_status = repmat("missing", height(T), 1);
T.preprocess_runtime_sec = nan(height(T), 1);
T.badframe_fraction = nan(height(T), 1);
auditVars = local_preprocess_audit_vars();
for i = 1:numel(auditVars)
    T.(char(auditVars(i))) = nan(height(T), 1);
end
T.preprocess_output_file = strings(height(T), 1);
T.preprocess_error_message = strings(height(T), 1);

[tf, loc] = ismember(T.raw_index, S.raw_index);
T.preprocess_status(tf) = S.status(loc(tf));
T.preprocess_runtime_sec(tf) = S.runtime_sec(loc(tf));
T.badframe_fraction(tf) = S.badframe_fraction(loc(tf));
for i = 1:numel(auditVars)
    v = char(auditVars(i));
    if ismember(v, S.Properties.VariableNames)
        T.(v)(tf) = S.(v)(loc(tf));
    end
end
T.preprocess_output_file(tf) = S.output_file(loc(tf));
if ismember('error_message', S.Properties.VariableNames)
    T.preprocess_error_message(tf) = S.error_message(loc(tf));
end

T.preprocess_success = T.preprocess_status == "success" | ...
    T.preprocess_status == "skipped_existing";
end

function vars = local_preprocess_audit_vars()
vars = [ ...
    "n_prediction_issue_animal_frames"; ...
    "n_interpolated_prediction_issue_animal_frames"; ...
    "n_usable_prediction_issue_animal_frames"; ...
    "n_repaired_prediction_issue_animal_frames"; ...
    "n_unresolved_prediction_issue_animal_frames"; ...
    "prediction_issue_fraction"; ...
    "interpolated_prediction_issue_fraction"; ...
    "usable_prediction_issue_fraction"; ...
    "repaired_prediction_issue_fraction"; ...
    "unresolved_prediction_issue_fraction"; ...
    "prediction_issue_usable_rate"; ...
    "prediction_issue_repair_rate"; ...
    "prediction_issue_unresolved_rate"; ...
    "n_prediction_issue_frames_any_animal"; ...
    "n_repaired_prediction_issue_frames_any_animal"; ...
    "n_unresolved_prediction_issue_frames_any_animal"; ...
    "prediction_issue_any_animal_fraction"; ...
    "repaired_prediction_issue_any_animal_fraction"; ...
    "unresolved_prediction_issue_any_animal_fraction"];
end

function vars = local_summary_metric_vars(T)
vars = { ...
    'badframe_fraction', ...
    'prediction_issue_fraction', ...
    'interpolated_prediction_issue_fraction', ...
    'repaired_prediction_issue_fraction', ...
    'unresolved_prediction_issue_fraction', ...
    'prediction_issue_repair_rate', ...
    'prediction_issue_unresolved_rate'};
vars = vars(ismember(vars, T.Properties.VariableNames));
end

function analysisClass = local_analysis_class(T)
analysisClass = repmat("excluded_or_unsupported", height(T), 1);
analysisClass(T.include_block2_egocentric == 1) = "block2_egocentric_context";
analysisClass(T.effective_n_animals == 1 & T.include_block2_egocentric == 1) = ...
    "block2_single_mouse";
analysisClass(select_analysis_cohort(T, "anesthetized_context")) = ...
    "block2_anesthetized_center";
analysisClass(select_analysis_cohort(T, "motif_discovery")) = ...
    "block1_motif_discovery";
end

function notes = local_qc_notes(T, motifMask, warnThresh, maxThresh)
notes = strings(height(T), 1);
for i = 1:height(T)
    parts = strings(0,1);
    if ~T.preprocess_success(i)
        parts(end+1,1) = "preprocessing_missing_or_failed"; 
    end
    if motifMask(i) && T.badframe_fraction(i) > maxThresh
        parts(end+1,1) = "exclude_motif_high_badframe_fraction";
    elseif motifMask(i) && T.badframe_fraction(i) > warnThresh
        parts(end+1,1) = "review_motif_badframe_fraction";
    end
    if string(T.motif_exclusion_reason(i)) == "anesthetized_center_special_egocentric_context"
        parts(end+1,1) = "reserved_for_block2_anesthetized_context"; 
    end
    notes(i) = strjoin(parts, ';');
end
end

function [animalTable, pointTable, detailAudit] = local_build_detail_tables(T, repoRoot, preprocRoot, cfg)
animalRows = struct([]);
pointRows = struct([]);
loaded = false(height(T), 1);
loadMessage = strings(height(T), 1);

for i = 1:height(T)
    if ~T.preprocess_success(i)
        loadMessage(i) = "preprocessing_not_successful";
        continue
    end

    outFile = local_resolve_output_file(T, i, repoRoot, preprocRoot, cfg);
    if ~isfile(outFile)
        loadMessage(i) = "missing_preprocessed_mat";
        continue
    end

    try
        S = load(outFile, 'sessionPreproc');
        assert(isfield(S, 'sessionPreproc'), 'Missing sessionPreproc');
        P = S.sessionPreproc;
        if ~isfield(P.qc, 'sessionStats') || ...
                ~isfield(P.qc.sessionStats, 'nPredictionIssueAnimalFrames')
            P.qc.sessionStats = summarize_preprocessing_qc(P);
        end
        [animalRows, pointRows] = local_append_session_detail_rows( ...
            animalRows, pointRows, T, i, P, cfg);
        loaded(i) = true;
        loadMessage(i) = "ok";
    catch ME
        loadMessage(i) = string(ME.message);
    end
end

animalTable = local_struct_to_table(animalRows);
pointTable = local_struct_to_table(pointRows);
detailAudit = table(T.raw_index, T.session_id, loaded, loadMessage, ...
    'VariableNames', {'raw_index','session_id','loaded_preprocessed_mat','load_message'});
end

function [animalRows, pointRows] = local_append_session_detail_rows(animalRows, pointRows, T, rowIdx, P, cfg)
nFrames = size(P.clean.tracks, 1);
nNodes = size(P.clean.tracks, 2);
nAnimals = size(P.clean.tracks, 4);
badframes = P.qc.badframes;
if isvector(badframes)
    badframes = badframes(:);
end
stats = P.qc.sessionStats;

for animal = 1:nAnimals
    base = local_base_row(T, rowIdx);
    base.animal = animal;
    base.n_frames = nFrames;
    base.n_nodes = nNodes;

    if isfield(stats, 'animal') && numel(stats.animal) >= animal
        st = stats.animal(animal);
    else
        st = struct();
    end

    animalRow = base;
    animalRow.pctBadframes = local_get_stat(st, 'pctBadframes');
    animalRow.pctInterpFrames = local_get_stat(st, 'pctInterpFrames');
    animalRow.pctJumpFrames = local_get_stat(st, 'pctJumpFrames');
    animalRow.pctGeomFrames = local_get_stat(st, 'pctGeomFrames');
    animalRow.pctArenaFrames = local_get_stat(st, 'pctArenaFrames');
    animalRow.pctLowConfFrames = local_get_stat(st, 'pctLowConfFrames');
    animalRow.pctJumpSamples = local_get_stat(st, 'pctJumpSamples');
    animalRow.pctInterpSamples = local_get_stat(st, 'pctInterpSamples');
    animalRow.nPredictionIssueFrames = local_get_stat(st, 'nPredictionIssueFrames');
    animalRow.nInterpolatedPredictionIssueFrames = local_get_stat(st, 'nInterpolatedPredictionIssueFrames');
    animalRow.nUsablePredictionIssueFrames = local_get_stat(st, 'nUsablePredictionIssueFrames');
    animalRow.nRepairedPredictionIssueFrames = local_get_stat(st, 'nRepairedPredictionIssueFrames');
    animalRow.nUnresolvedPredictionIssueFrames = local_get_stat(st, 'nUnresolvedPredictionIssueFrames');
    animalRow.pctPredictionIssueFrames = local_get_stat(st, 'pctPredictionIssueFrames');
    animalRow.pctPredictionIssueUsable = local_get_stat(st, 'pctPredictionIssueUsable');
    animalRow.pctPredictionIssueRepaired = local_get_stat(st, 'pctPredictionIssueRepaired');
    animalRow.pctPredictionIssueUnresolved = local_get_stat(st, 'pctPredictionIssueUnresolved');
    animalRows = [animalRows; animalRow]; 

    qa = P.qc.animals(animal);
    masks = local_full_qc_masks(qa, nFrames, nNodes);
    animalBad = badframes(:, min(animal, size(badframes, 2)));
    for node = 1:nNodes
        issue = masks.lowConf(:,node) | masks.jump(:,node) | masks.geom(:,node);
        repaired = issue & masks.interp(:,node) & ~animalBad & ~masks.finalNan(:,node);
        unresolved = issue & (animalBad | masks.finalNan(:,node));

        pointRow = base;
        pointRow.node = node;
        pointRow.node_label = local_node_label(cfg, node);
        pointRow.n_prediction_issue_samples = nnz(issue);
        pointRow.n_interpolated_prediction_issue_samples = nnz(issue & masks.interp(:,node));
        pointRow.n_repaired_prediction_issue_samples = nnz(repaired);
        pointRow.n_unresolved_prediction_issue_samples = nnz(unresolved);
        pointRow.fraction_low_confidence_samples = mean(masks.lowConf(:,node));
        pointRow.fraction_jump_samples = mean(masks.jump(:,node));
        pointRow.fraction_geometry_samples = mean(masks.geom(:,node));
        pointRow.fraction_arena_samples = mean(masks.arena(:,node));
        pointRow.fraction_interpolated_samples = mean(masks.interp(:,node));
        pointRow.fraction_final_nan_samples = mean(masks.finalNan(:,node));
        pointRow.fraction_prediction_issue_samples = mean(issue);
        pointRow.fraction_repaired_prediction_issue_samples = mean(repaired);
        pointRow.fraction_unresolved_prediction_issue_samples = mean(unresolved);
        pointRow.prediction_issue_repair_rate = local_fraction(nnz(repaired), nnz(issue));
        pointRow.prediction_issue_unresolved_rate = local_fraction(nnz(unresolved), nnz(issue));
        pointRows = [pointRows; pointRow]; 
    end
end
end

function base = local_base_row(T, rowIdx)
base = struct();
base.raw_index = T.raw_index(rowIdx);
base.session_id = string(T.session_id(rowIdx));
base.condition = string(T.condition(rowIdx));
base.condition_id = local_table_string(T, 'condition_id', rowIdx, string(T.condition(rowIdx)));
base.condition_label = local_table_string(T, 'condition_label', rowIdx, string(T.condition(rowIdx)));
base.arena = string(T.arena(rowIdx));
base.arena_condition = local_table_string(T, 'arena_condition', rowIdx, "");
base.arena_condition_label = local_table_string(T, 'arena_condition_label', rowIdx, "");
base.analysis_group_id = local_table_string(T, 'analysis_group_id', rowIdx, "");
base.analysis_group_label = local_table_string(T, 'analysis_group_label', rowIdx, "");
base.plot_group_label = local_table_string(T, 'plot_group_label', rowIdx, "");
base.analysis_class = string(T.analysis_class(rowIdx));
end

function masks = local_full_qc_masks(qa, nFrames, nNodes)
masks.lowConf = local_mask_field(qa, 'lowConfMask', nFrames, nNodes);
masks.jump = local_mask_field(qa, 'jumpMask', nFrames, nNodes);
masks.interp = local_mask_field(qa, 'interpMask', nFrames, nNodes);
masks.geom = local_mask_field(qa, 'geomMask', nFrames, nNodes);
masks.arena = local_mask_field(qa, 'arenaMask', nFrames, nNodes);
masks.finalNan = local_mask_field(qa, 'finalNanMask', nFrames, nNodes);
end

function mask = local_mask_field(qa, fieldName, nFrames, nNodes)
if isfield(qa, fieldName) && ~isempty(qa.(fieldName))
    mask = logical(full(qa.(fieldName)));
else
    mask = false(nFrames, nNodes);
end
if ~isequal(size(mask), [nFrames nNodes])
    mask = false(nFrames, nNodes);
end
end

function label = local_node_label(cfg, node)
labels = cfg.bodypoints.bodypoint_labels;
if node <= numel(labels)
    label = string(labels(node));
else
    label = "Point " + string(sprintf('%02d', node));
end
end

function summary = local_animal_group_summary(animalTable)
if isempty(animalTable) || height(animalTable) == 0
    summary = table();
    return
end
metrics = {'pctBadframes','pctPredictionIssueFrames','pctPredictionIssueRepaired', ...
    'pctPredictionIssueUnresolved','pctInterpFrames','pctLowConfFrames', ...
    'pctJumpFrames','pctGeomFrames','pctArenaFrames'};
metrics = metrics(ismember(metrics, animalTable.Properties.VariableNames));
summary = groupsummary(animalTable, {'analysis_class','arena_condition','condition'}, ...
    {'mean','median','max'}, metrics);
end

function summary = local_point_summary(pointTable)
if isempty(pointTable) || height(pointTable) == 0
    summary = table();
    return
end
metrics = {'fraction_low_confidence_samples','fraction_jump_samples', ...
    'fraction_geometry_samples','fraction_interpolated_samples', ...
    'fraction_final_nan_samples','fraction_prediction_issue_samples', ...
    'fraction_repaired_prediction_issue_samples', ...
    'fraction_unresolved_prediction_issue_samples'};
summary = groupsummary(pointTable, {'node','node_label'}, {'mean','median','max'}, metrics);
summary = sortrows(summary, 'node');
end

function figureManifest = local_make_figures(T, animalTable, ~, pointSummary, styles, cfg, figDir)
rows = struct([]);
rows = local_add_figure_row(rows, local_plot_story_overview(T, styles, cfg, figDir), ...
    'Overall preprocessing story: prediction flags, repair, unresolved burden, and final bad frames.');
rows = local_add_figure_row(rows, local_plot_point_burden(pointSummary, cfg, figDir), ...
    'Point-level QC burden across SLEAP body points.');
rows = local_add_figure_row(rows, local_plot_animal_repair(animalTable, styles, cfg, figDir), ...
    'Animal-level prediction-issue repair rates by non-pooled analysis group.');
rows = local_add_figure_row(rows, local_plot_session_rank(T, styles, cfg, figDir), ...
    'Dyadic/session-level QC ranking with stable condition colors.');
rows = local_add_figure_row(rows, local_plot_condition_arena_summary(T, styles, cfg, figDir), ...
    'Condition-by-arena QC summary without pooling arenas.');
rows = local_add_figure_row(rows, local_plot_color_key(styles, cfg, figDir), ...
    'Stable experimental-group color key used by paper pipeline figures.');
figureManifest = struct2table(rows);
end

function rows = local_add_figure_row(rows, files, description)
if isempty(files)
    return
end
for i = 1:numel(files)
    row = struct();
    row.figure_file = string(files(i));
    row.description = string(description);
    rows = [rows; row]; %#ok<AGROW>
end
end

function files = local_plot_story_overview(T, styles, cfg, figDir)
ok = isfinite(T.badframe_fraction);
if ~any(ok)
    files = strings(0,1);
    return
end

fig = figure('Visible','off', 'Color','w', 'Position',[80 80 1500 950]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');

nexttile(tl);
metricNames = ["prediction_issue_fraction","repaired_prediction_issue_fraction", ...
    "unresolved_prediction_issue_fraction","badframe_fraction"];
metricLabels = ["Prediction issue","Repaired","Unresolved","Final bad"];
stageColors = [ ...
    local_cfg_rgb(cfg, 'story_color_prediction_issue'); ...
    local_cfg_rgb(cfg, 'story_color_repaired'); ...
    local_cfg_rgb(cfg, 'story_color_unresolved'); ...
    local_cfg_rgb(cfg, 'story_color_badframe')];
values = nan(1, numel(metricNames));
lo = values;
hi = values;
for i = 1:numel(metricNames)
    x = T.(metricNames(i));
    x = x(isfinite(x)) .* 100;
    values(i) = median(x);
    lo(i) = prctile(x, 25);
    hi(i) = prctile(x, 75);
end
b = bar(values, 'FaceColor','flat');
b.CData = stageColors;
hold on;
errorbar(1:numel(values), values, values-lo, hi-values, 'k.', 'LineWidth', 0.9);
set(gca, 'XTick', 1:numel(metricLabels), 'XTickLabel', metricLabels);
ylabel('% animal-frames');
title('A. What preprocessing encountered and retained');
local_style_axes(gca, cfg);

nexttile(tl);
hold on;
conds = unique(string(T.condition_id(ok)), 'stable');
for i = 1:numel(conds)
    idx = ok & string(T.condition_id) == conds(i);
    scatter(T.prediction_issue_fraction(idx).*100, T.badframe_fraction(idx).*100, ...
        cfg.figures.marker_size, local_condition_rgb(styles, conds(i)), 'filled');
end
plot([0 100], [0 100], '--', 'Color', local_cfg_rgb(cfg, 'grid_color'), 'LineWidth', 1);
xlabel('% animal-frames with prediction issue');
ylabel('% final bad animal-frames');
title('B. Flagged prediction issues are audited separately from final exclusions');
xlim([0 max(65, max(T.prediction_issue_fraction(ok).*100) * 1.05)]);
ylim([0 max(60, max(T.badframe_fraction(ok).*100) * 1.05)]);
local_style_axes(gca, cfg);

nexttile(tl);
repair = T.prediction_issue_repair_rate(isfinite(T.prediction_issue_repair_rate));
usable = T.prediction_issue_usable_rate(isfinite(T.prediction_issue_usable_rate));
unresolved = T.prediction_issue_unresolved_rate(isfinite(T.prediction_issue_unresolved_rate));
usableNoInterp = max(median(usable) - median(repair), 0);
stackVals = 100 .* [median(repair), usableNoInterp, median(unresolved)];
b = bar(1, stackVals, 'stacked');
b(1).FaceColor = local_cfg_rgb(cfg, 'story_color_repaired');
b(2).FaceColor = local_cfg_rgb(cfg, 'story_color_interpolated');
b(3).FaceColor = local_cfg_rgb(cfg, 'story_color_unresolved');
set(gca, 'XTick', 1, 'XTickLabel', "Prediction-issue frames");
ylabel('% of prediction-issue animal-frames');
legend({'repaired and usable','usable without interpolation','unresolved'}, ...
    'Location','southoutside', 'Orientation','horizontal', 'Box','off');
title('C. Fate of prediction-issue frames');
ylim([0 100]);
local_style_axes(gca, cfg);

nexttile(tl);
motifMask = select_analysis_cohort(T, "motif_discovery");
reviewMask = T.qc_warn_motif_discovery;
passCleanMask = T.qc_pass_motif_discovery & ~reviewMask;
failMask = motifMask & T.preprocess_success & ~T.qc_pass_motif_discovery;
counts = [nnz(passCleanMask), nnz(reviewMask), nnz(failMask)];
b = bar(counts, 'FaceColor','flat');
b.CData = [local_cfg_rgb(cfg, 'story_color_repaired'); ...
    local_cfg_rgb(cfg, 'story_color_interpolated'); ...
    local_cfg_rgb(cfg, 'story_color_badframe')];
set(gca, 'XTick', 1:3, 'XTickLabel', ["pass","review","fail"]);
text(1:3, counts + max(counts) * 0.03, string(counts), ...
    'HorizontalAlignment','center', 'FontName', char(cfg.figures.font_name), ...
    'FontSize', cfg.figures.font_size);
ylabel('Motif-discovery sessions');
title('D. Motif-discovery QC decision');
ylim([0 max(counts) * 1.15]);
local_style_axes(gca, cfg);

files = local_export_figure(fig, figDir, 'preprocess_qc_story_overview', cfg);
close(fig);
end

function files = local_plot_point_burden(pointSummary, cfg, figDir)
if isempty(pointSummary) || height(pointSummary) == 0
    files = strings(0,1);
    return
end

metricVars = {'median_fraction_low_confidence_samples','median_fraction_jump_samples', ...
    'median_fraction_geometry_samples','median_fraction_interpolated_samples', ...
    'median_fraction_final_nan_samples','median_fraction_repaired_prediction_issue_samples', ...
    'median_fraction_unresolved_prediction_issue_samples'};
metricLabels = {'Low confidence','Jump','Geometry','Interpolated','Final NaN', ...
    'Repaired issue','Unresolved issue'};
X = zeros(height(pointSummary), numel(metricVars));
for j = 1:numel(metricVars)
    if ismember(metricVars{j}, pointSummary.Properties.VariableNames)
        X(:,j) = pointSummary.(metricVars{j}) .* 100;
    else
        X(:,j) = NaN;
    end
end

fig = figure('Visible','off', 'Color','w', 'Position',[80 80 1100 850]);
imagesc(X);
colormap(parula);
cb = colorbar;
cb.Label.String = 'Median % point samples';
set(gca, 'XTick', 1:numel(metricLabels), 'XTickLabel', metricLabels, ...
    'YTick', 1:height(pointSummary), 'YTickLabel', pointSummary.node_label);
xtickangle(35);
ylabel('SLEAP point');
title('Point-level preprocessing burden');
local_style_axes(gca, cfg);

files = local_export_figure(fig, figDir, 'preprocess_qc_point_level_bodypoints', cfg);
close(fig);
end

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

function G = local_sort_group_summary(G, styles)
conditionOrder = nan(height(G), 1);
for i = 1:height(G)
    conditionOrder(i) = find(styles.conditionIds == string(G.condition_id(i)), 1, 'first');
end
conditionOrder(~isfinite(conditionOrder)) = height(styles.table) + 1;
G.condition_order = conditionOrder;
sortVars = {'condition_order'};
if ismember('arena_condition_label', G.Properties.VariableNames)
    sortVars = [sortVars, {'arena_condition_label'}];
elseif ismember('plot_group_label', G.Properties.VariableNames)
    sortVars = [sortVars, {'plot_group_label'}];
end
G = sortrows(G, sortVars);
end

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

function local_style_axes(ax, cfg)
set(ax, 'FontName', char(cfg.figures.font_name), 'FontSize', cfg.figures.font_size, ...
    'LineWidth', 0.8, 'Box','off', 'TickDir','out');
grid(ax, 'on');
ax.GridColor = local_cfg_rgb(cfg, 'grid_color');
ax.GridAlpha = 0.35;
title(ax, ax.Title.String, 'FontSize', cfg.figures.title_font_size, ...
    'FontWeight','bold');
end

function rgb = local_condition_rgb(styles, conditionId)
idx = find(styles.conditionIds == string(conditionId), 1, 'first');
if isempty(idx)
    rgb = [0.35 0.35 0.35];
else
    rgb = styles.colors(idx,:);
end
end

function rgb = local_cfg_rgb(cfg, key)
rgb = local_hex_to_rgb(cfg.figures.(key));
end

function rgb = local_hex_to_rgb(hex)
hex = char(strtrim(string(hex)));
rgb = sscanf(hex(2:end), '%2x%2x%2x', [1 3]) ./ 255;
end

function y = local_jitter(n, cfg)
if n == 0
    y = [];
else
    y = linspace(-cfg.figures.jitter_width, cfg.figures.jitter_width, n)';
end
end

function pathOut = local_resolve_output_file(T, rowIdx, repoRoot, preprocRoot, cfg)
candidate = string(T.preprocess_output_file(rowIdx));
if strlength(candidate) > 0 && isfile(candidate)
    pathOut = char(candidate);
    return
end

sessionDir = fullfile(preprocRoot, char(cfg.paths.preprocess_session_subdir));
pathOut = fullfile(sessionDir, sprintf('session_%04d_preprocessed.mat', T.raw_index(rowIdx)));
if isfile(pathOut)
    return
end

if strlength(candidate) > 0
    pathOut = local_repo_path(repoRoot, candidate);
end
end

function pathOut = local_repo_path(repoRoot, relativePath)
relativePath = char(relativePath);
if startsWith(relativePath, filesep) || ~isempty(regexp(relativePath, '^[A-Za-z]:[\\/]', 'once'))
    pathOut = relativePath;
else
    pathOut = fullfile(repoRoot, relativePath);
end
end

function v = local_get_stat(stats, fieldName)
if isstruct(stats) && isfield(stats, fieldName)
    v = stats.(fieldName);
else
    v = NaN;
end
end

function v = local_table_string(T, varName, rowIdx, defaultValue)
if ismember(varName, T.Properties.VariableNames)
    v = string(T.(varName)(rowIdx));
else
    v = string(defaultValue);
end
end

function f = local_fraction(num, den)
if den > 0
    f = num ./ den;
else
    f = NaN;
end
end

function T = local_struct_to_table(rows)
if isempty(rows)
    T = table();
else
    T = struct2table(rows);
end
end
