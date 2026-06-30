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

outRoot = resolve_repo_path(repoRoot, cfg.paths.preprocessed_output_dir);
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
