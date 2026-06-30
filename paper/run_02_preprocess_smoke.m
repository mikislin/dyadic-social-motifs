%RUN_02_PREPROCESS_SMOKE Preprocess a small representative real-data set.
%
% This is the first real-data gate after manifest generation. It intentionally
% includes one reduced-from-three-tracks session to verify animal QC.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
cd(repoRoot);
addpath(genpath(repoRoot));

paths = paper_paths('requireRawData', true);
if ~isfile(paths.manifestPath)
    run(fullfile(repoRoot, 'paper', 'run_01_manifest.m'));
end

M = readtable(paths.manifestPath, 'TextType', 'string');
outDir = fullfile(repoRoot, 'derived', 'preprocess_smoke');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

labels = ["single-wt", "wt-wt", "m-wt", "m-m", "reduced-three-wt-wt"];
rowIdx = nan(numel(labels),1);
rowIdx(1) = first_manifest_row(M, M.effective_n_animals == 1 & M.condition == "WT", ...
    'run_02_preprocess_smoke:MissingSmokeCase', labels(1));
rowIdx(2) = first_manifest_row(M, M.effective_n_animals == 2 & M.condition == "WT_WT" & M.animal_qc_status == "ok", ...
    'run_02_preprocess_smoke:MissingSmokeCase', labels(2));
rowIdx(3) = first_manifest_row(M, M.effective_n_animals == 2 & M.condition == "M_WT", ...
    'run_02_preprocess_smoke:MissingSmokeCase', labels(3));
rowIdx(4) = first_manifest_row(M, M.effective_n_animals == 2 & M.condition == "M_M", ...
    'run_02_preprocess_smoke:MissingSmokeCase', labels(4));
rowIdx(5) = first_manifest_row(M, M.effective_n_animals == 2 & M.animal_qc_status == "reduce_to_dyad", ...
    'run_02_preprocess_smoke:MissingSmokeCase', labels(5));

P = default_preprocessing_params();
P.output.return_raw = false;
P.output.make_plots = false;
P.debug.verbose = false;

summary = table();
for i = 1:numel(rowIdx)
    r = rowIdx(i);
    [sessionRaw, row, animalQC] = load_session_raw(M, r, ...
        'applyAnimalQC', true, ...
        'sessionRoot', paths.files2runDir);
    fprintf('Smoke preprocessing %s: %s (%s), raw animals=%d, effective=%d, keep=%s\n', ...
        labels(i), row.session_file, row.condition, animalQC.n_animals_raw, ...
        animalQC.n_animals_effective, strjoin(string(animalQC.keep_indices), ';'));

    sessionPreproc = preprocess_session(sessionRaw, P);
    outPath = fullfile(outDir, labels(i) + "_" + row.session_file);
    save(outPath, 'sessionPreproc', 'row', 'animalQC', '-v7.3');

    stats = sessionPreproc.qc.sessionStats;
    one = table();
    one.label = labels(i);
    one.raw_index = row.raw_index;
    one.session_file = row.session_file;
    one.session_id = row.session_id;
    one.condition = row.condition;
    one.arena = row.arena;
    one.n_animals_raw = animalQC.n_animals_raw;
    one.effective_n_animals = animalQC.n_animals_effective;
    one.animal_keep_indices = string(strjoin(string(animalQC.keep_indices), ';'));
    one.animal_qc_status = string(animalQC.status);
    one.n_frames = row.n_frames;
    if isfield(stats, 'badframeFraction')
        one.badframe_fraction = stats.badframeFraction;
    else
        one.badframe_fraction = mean(sessionPreproc.qc.badframes(:));
    end
    summary = [summary; one]; %#ok<AGROW>
end

summaryPath = fullfile(outDir, 'preprocess_smoke_summary.csv');
writetable(summary, summaryPath);
disp(summary);
fprintf('Wrote smoke preprocessing summary: %s\n', summaryPath);
