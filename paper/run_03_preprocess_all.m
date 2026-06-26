%RUN_03_PREPROCESS_ALL Preprocess all manifest-supported sessions.
%
% Outputs one file per session so the run can resume after interruption.
% Existing successful outputs are skipped by default.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
cd(repoRoot);
addpath(genpath(repoRoot));

paths = paper_paths('requireRawData', true);
cfg = load_preprocessing_pipeline_config(paths.preprocessingPipelineConfigPath);
if ~isfile(paths.manifestPath)
    run(fullfile(repoRoot, 'paper', 'run_01_manifest.m'));
end

M = readtable(paths.manifestPath, 'TextType', 'string');

outRoot = local_repo_path(repoRoot, cfg.paths.preprocessed_output_dir);
sessionOutDir = fullfile(outRoot, char(cfg.paths.preprocess_session_subdir));
logDir = fullfile(outRoot, 'logs');
if ~exist(sessionOutDir, 'dir'); mkdir(sessionOutDir); end
if ~exist(logDir, 'dir'); mkdir(logDir); end

diary(fullfile(logDir, 'run_03_preprocess_all_latest.log'));
cleanup = onCleanup(@() diary('off'));

skipExisting = cfg.run_03.skip_existing;
targetMask = M.paper_include == 1;
idx = find(targetMask);

P = default_preprocessing_params();
P.output.return_raw = cfg.run_03.return_raw;
P.output.make_plots = cfg.run_03.make_plots;
P.output.store_full_masks = cfg.run_03.store_full_masks;
P.debug.verbose = cfg.run_03.debug_verbose;

fprintf('run_03_preprocess_all\n');
fprintf('Repo root: %s\n', repoRoot);
fprintf('Output root: %s\n', outRoot);
fprintf('Manifest rows: %d\n', height(M));
fprintf('Target sessions: %d\n', numel(idx));
fprintf('Skip existing successful outputs: %d\n', skipExisting);

summary = initialize_summary(M, idx);
summaryPath = fullfile(outRoot, char(cfg.paths.preprocess_summary_file));

for k = 1:numel(idx)
    rowIdx = idx(k);
    row = M(rowIdx,:);
    outFile = fullfile(sessionOutDir, sprintf('session_%04d_preprocessed.mat', row.raw_index));

    fprintf('\n[%d/%d] raw_index=%d file=%s condition=%s animals raw/effective=%d/%d\n', ...
        k, numel(idx), row.raw_index, row.session_file, row.condition, ...
        row.n_animals, row.effective_n_animals);

    if skipExisting && isfile(outFile)
        try
            S = load(outFile, 'status', 'runtime_sec', 'badframe_fraction', ...
                'preprocess_qc_audit', 'sessionPreproc');
            if isfield(S, 'status') && string(S.status) == "success"
                fprintf('  skipping existing successful output\n');
                audit = local_audit_from_saved_output(S);
                summary.status(k) = "skipped_existing";
                if isfield(S, 'runtime_sec')
                    summary.runtime_sec(k) = S.runtime_sec;
                end
                summary.output_file(k) = local_repo_relative_path(repoRoot, outFile);
                summary = local_assign_audit_summary(summary, k, audit);
                writetable(summary, summaryPath);
                continue
            end
        catch
            fprintf('  existing output could not be read; recomputing\n');
        end
    end

    tStart = tic;
    try
        [sessionRaw, manifestRow, animalQC] = load_session_raw(M, rowIdx, ...
            'applyAnimalQC', true, ...
            'sessionRoot', paths.files2runDir); 
        sessionPreproc = preprocess_session(sessionRaw, P); 
        runtime_sec = toc(tStart); 
        status = "success"; 

        preprocess_qc_audit = local_preprocessing_audit(sessionPreproc);
        badframe_fraction = preprocess_qc_audit.badframe_fraction;

        save(outFile, 'sessionPreproc', 'manifestRow', 'animalQC', ...
            'P', 'status', 'runtime_sec', 'badframe_fraction', ...
            'preprocess_qc_audit', '-v7.3');

        summary.status(k) = "success";
        summary.runtime_sec(k) = runtime_sec;
        summary = local_assign_audit_summary(summary, k, preprocess_qc_audit);
        summary.output_file(k) = local_repo_relative_path(repoRoot, outFile);
        fprintf(['  success runtime=%.1fs badframe_fraction=%.4f ' ...
            'prediction_issue_repair_rate=%.4f\n'], ...
            runtime_sec, badframe_fraction, ...
            preprocess_qc_audit.prediction_issue_repair_rate);
    catch ME
        runtime_sec = toc(tStart); 
        status = "failed"; 
        error_message = string(ME.message); 
        save(outFile, 'row', 'P', 'status', 'runtime_sec', 'error_message', '-v7.3');

        summary.status(k) = "failed";
        summary.runtime_sec(k) = runtime_sec;
        summary.error_message(k) = string(ME.message);
        summary.output_file(k) = local_repo_relative_path(repoRoot, outFile);
        fprintf('  FAILED runtime=%.1fs error=%s\n', runtime_sec, ME.message);
    end

    writetable(summary, summaryPath);
end

writetable(summary, summaryPath);
disp(groupsummary(summary, 'status'));
fprintf('Wrote summary: %s\n', summaryPath);

nFailed = nnz(summary.status == "failed");
assert(nFailed == 0, 'run_03_preprocess_all:Failures', ...
    'Preprocessing finished with %d failed sessions. See %s.', nFailed, summaryPath);

function summary = initialize_summary(M, idx)
summary = table();
summary.raw_index = M.raw_index(idx);
summary.session_file = M.session_file(idx);
summary.session_id = M.session_id(idx);
summary.condition = M.condition(idx);
summary.arena = M.arena(idx);
summary.n_animals = M.n_animals(idx);
summary.effective_n_animals = M.effective_n_animals(idx);
summary.animal_qc_status = M.animal_qc_status(idx);
summary.status = repmat("pending", numel(idx), 1);
summary.runtime_sec = nan(numel(idx),1);
summary.badframe_fraction = nan(numel(idx),1);
summary.n_prediction_issue_animal_frames = nan(numel(idx),1);
summary.n_interpolated_prediction_issue_animal_frames = nan(numel(idx),1);
summary.n_usable_prediction_issue_animal_frames = nan(numel(idx),1);
summary.n_repaired_prediction_issue_animal_frames = nan(numel(idx),1);
summary.n_unresolved_prediction_issue_animal_frames = nan(numel(idx),1);
summary.prediction_issue_fraction = nan(numel(idx),1);
summary.interpolated_prediction_issue_fraction = nan(numel(idx),1);
summary.usable_prediction_issue_fraction = nan(numel(idx),1);
summary.repaired_prediction_issue_fraction = nan(numel(idx),1);
summary.unresolved_prediction_issue_fraction = nan(numel(idx),1);
summary.prediction_issue_usable_rate = nan(numel(idx),1);
summary.prediction_issue_repair_rate = nan(numel(idx),1);
summary.prediction_issue_unresolved_rate = nan(numel(idx),1);
summary.n_prediction_issue_frames_any_animal = nan(numel(idx),1);
summary.n_repaired_prediction_issue_frames_any_animal = nan(numel(idx),1);
summary.n_unresolved_prediction_issue_frames_any_animal = nan(numel(idx),1);
summary.prediction_issue_any_animal_fraction = nan(numel(idx),1);
summary.repaired_prediction_issue_any_animal_fraction = nan(numel(idx),1);
summary.unresolved_prediction_issue_any_animal_fraction = nan(numel(idx),1);
summary.output_file = strings(numel(idx),1);
summary.error_message = strings(numel(idx),1);
end

function pathOut = local_repo_path(repoRoot, relativePath)
relativePath = char(relativePath);
if isfolder(relativePath) || startsWith(relativePath, filesep) || ...
        ~isempty(regexp(relativePath, '^[A-Za-z]:[\\/]', 'once'))
    pathOut = relativePath;
else
    pathOut = fullfile(repoRoot, relativePath);
end
end

function relPath = local_repo_relative_path(repoRoot, absPath)
repoRoot = char(repoRoot);
absPath = char(absPath);
if startsWith(absPath, repoRoot)
    relPath = string(extractAfter(string(absPath), strlength(string(repoRoot)) + 1));
else
    relPath = string(absPath);
end
relPath = replace(relPath, '\', '/');
end

function summary = local_assign_audit_summary(summary, rowIdx, audit)
fields = fieldnames(audit);
for i = 1:numel(fields)
    if ismember(fields{i}, summary.Properties.VariableNames)
        summary.(fields{i})(rowIdx) = audit.(fields{i});
    end
end
end

function audit = local_audit_from_saved_output(S)
if isfield(S, 'preprocess_qc_audit')
    audit = S.preprocess_qc_audit;
    if ~isfield(audit, 'badframe_fraction') && isfield(S, 'badframe_fraction')
        audit.badframe_fraction = S.badframe_fraction;
    end
elseif isfield(S, 'sessionPreproc')
    audit = local_preprocessing_audit(S.sessionPreproc);
elseif isfield(S, 'badframe_fraction')
    audit = local_empty_audit();
    audit.badframe_fraction = S.badframe_fraction;
else
    audit = local_empty_audit();
end
end

function audit = local_preprocessing_audit(sessionPreproc)
if ~isfield(sessionPreproc.qc, 'sessionStats') || ...
        ~isfield(sessionPreproc.qc.sessionStats, 'nPredictionIssueAnimalFrames')
    sessionPreproc.qc.sessionStats = summarize_preprocessing_qc(sessionPreproc);
end

stats = sessionPreproc.qc.sessionStats;
audit = local_empty_audit();
audit.badframe_fraction = local_badframe_fraction(sessionPreproc, stats);
audit.n_prediction_issue_animal_frames = local_get_stat(stats, 'nPredictionIssueAnimalFrames');
audit.n_interpolated_prediction_issue_animal_frames = local_get_stat(stats, 'nInterpolatedPredictionIssueAnimalFrames');
audit.n_usable_prediction_issue_animal_frames = local_get_stat(stats, 'nUsablePredictionIssueAnimalFrames');
audit.n_repaired_prediction_issue_animal_frames = local_get_stat(stats, 'nRepairedPredictionIssueAnimalFrames');
audit.n_unresolved_prediction_issue_animal_frames = local_get_stat(stats, 'nUnresolvedPredictionIssueAnimalFrames');
audit.prediction_issue_fraction = local_get_stat(stats, 'predictionIssueFraction');
audit.interpolated_prediction_issue_fraction = local_get_stat(stats, 'interpolatedPredictionIssueFraction');
audit.usable_prediction_issue_fraction = local_get_stat(stats, 'usablePredictionIssueFraction');
audit.repaired_prediction_issue_fraction = local_get_stat(stats, 'repairedPredictionIssueFraction');
audit.unresolved_prediction_issue_fraction = local_get_stat(stats, 'unresolvedPredictionIssueFraction');
audit.prediction_issue_usable_rate = local_get_stat(stats, 'predictionIssueUsableRate');
audit.prediction_issue_repair_rate = local_get_stat(stats, 'predictionIssueRepairRate');
audit.prediction_issue_unresolved_rate = local_get_stat(stats, 'predictionIssueUnresolvedRate');
audit.n_prediction_issue_frames_any_animal = local_get_stat(stats, 'nPredictionIssueFramesAnyAnimal');
audit.n_repaired_prediction_issue_frames_any_animal = local_get_stat(stats, 'nRepairedPredictionIssueFramesAnyAnimal');
audit.n_unresolved_prediction_issue_frames_any_animal = local_get_stat(stats, 'nUnresolvedPredictionIssueFramesAnyAnimal');
audit.prediction_issue_any_animal_fraction = local_get_stat(stats, 'predictionIssueAnyAnimalFraction');
audit.repaired_prediction_issue_any_animal_fraction = local_get_stat(stats, 'repairedPredictionIssueAnyAnimalFraction');
audit.unresolved_prediction_issue_any_animal_fraction = local_get_stat(stats, 'unresolvedPredictionIssueAnyAnimalFraction');
end

function audit = local_empty_audit()
audit = struct( ...
    'badframe_fraction', NaN, ...
    'n_prediction_issue_animal_frames', NaN, ...
    'n_interpolated_prediction_issue_animal_frames', NaN, ...
    'n_usable_prediction_issue_animal_frames', NaN, ...
    'n_repaired_prediction_issue_animal_frames', NaN, ...
    'n_unresolved_prediction_issue_animal_frames', NaN, ...
    'prediction_issue_fraction', NaN, ...
    'interpolated_prediction_issue_fraction', NaN, ...
    'usable_prediction_issue_fraction', NaN, ...
    'repaired_prediction_issue_fraction', NaN, ...
    'unresolved_prediction_issue_fraction', NaN, ...
    'prediction_issue_usable_rate', NaN, ...
    'prediction_issue_repair_rate', NaN, ...
    'prediction_issue_unresolved_rate', NaN, ...
    'n_prediction_issue_frames_any_animal', NaN, ...
    'n_repaired_prediction_issue_frames_any_animal', NaN, ...
    'n_unresolved_prediction_issue_frames_any_animal', NaN, ...
    'prediction_issue_any_animal_fraction', NaN, ...
    'repaired_prediction_issue_any_animal_fraction', NaN, ...
    'unresolved_prediction_issue_any_animal_fraction', NaN);
end

function v = local_get_stat(stats, fieldName)
if isstruct(stats) && isfield(stats, fieldName)
    v = stats.(fieldName);
else
    v = NaN;
end
end

function badframe_fraction = local_badframe_fraction(sessionPreproc, stats)
if isstruct(stats) && isfield(stats, 'badframeFraction')
    badframe_fraction = stats.badframeFraction;
else
    badframe_fraction = mean(sessionPreproc.qc.badframes(:));
end
end
