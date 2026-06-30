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
preprocRoot = resolve_repo_path(repoRoot, cfg.paths.preprocessed_output_dir);
summaryPath = fullfile(preprocRoot, char(cfg.paths.preprocess_summary_file));
outDir = resolve_repo_path(repoRoot, cfg.paths.preprocess_qc_output_dir);
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
