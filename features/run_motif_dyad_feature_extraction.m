function outputs = run_motif_dyad_feature_extraction(repoRoot, opts)
%RUN_MOTIF_DYAD_FEATURE_EXTRACTION Extract dyadic motif-discovery features.
%
% Condition-blind contract
%   - Session selection uses motif-discovery QC pass, preprocessing success,
%     effective dyad status, and arena-balanced smoke coverage only.
%   - Condition/cohort labels are carried only in audit tables for later
%     statistical joins. They are not used in feature computation, missingness
%     summaries, transforms, normalization, chunking, embedding, clustering,
%     topology, or motif-family definition.
%   - Small and large arenas are processed by the same feature code and remain
%     separable only in provenance and QC/distribution summaries.
%
% Parameters are loaded from config/motif_dyad_feature_extraction_config.csv.

if nargin < 1 || isempty(repoRoot)
    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
if nargin < 2 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'configPath') || strlength(string(opts.configPath)) == 0
    opts.configPath = fullfile(repoRoot, 'config', 'motif_dyad_feature_extraction_config.csv');
end

cd(repoRoot);
addpath(genpath(repoRoot));

paths = paper_paths();
cfg = load_preprocessing_pipeline_config(paths.preprocessingPipelineConfigPath);

qcDir = resolve_repo_path(repoRoot, cfg.paths.preprocess_qc_output_dir);
qcPath = fullfile(qcDir, 'preprocess_qc_session_table.csv');
if ~isfile(qcPath)
    run(fullfile(repoRoot, 'paper', 'run_04_preprocess_qc_review.m'));
end
assert(isfile(qcPath), 'run_05_extract_motif_dyad_features:MissingQC', ...
    'Missing preprocessing QC table: %s', qcPath);

outRoot = fullfile(paths.derivedDir, 'features_motif_discovery');
sessionOutDir = fullfile(outRoot, 'sessions');
figDir = fullfile(outRoot, 'figures');
logDir = fullfile(outRoot, 'logs');
if ~exist(outRoot, 'dir'); mkdir(outRoot); end
if ~exist(sessionOutDir, 'dir'); mkdir(sessionOutDir); end
if ~exist(figDir, 'dir'); mkdir(figDir); end
if ~exist(logDir, 'dir'); mkdir(logDir); end

diary(fullfile(logDir, 'run_05_extract_motif_dyad_features_latest.log'));
cleanup = onCleanup(@() diary('off'));

params = load_motif_dyad_feature_extraction_config(opts.configPath);
[nodeMap, partNames, nodeMetadata] = default_sleap_node_map(paths.preprocessingPipelineConfigPath);
[featureNames, featureMeta] = default_dyad_feature_metadata();

parameterAuditPath = fullfile(outRoot, 'feature_extraction_parameter_audit.csv');
featureDictionaryPath = fullfile(outRoot, 'feature_dictionary.csv');
sessionAuditPath = fullfile(outRoot, 'feature_extraction_session_audit.csv');
missingnessPath = fullfile(outRoot, 'feature_missingness_validity_summary.csv');
distributionPath = fullfile(outRoot, 'feature_distribution_summary.csv');
distributionDensityPath = fullfile(outRoot, 'feature_distribution_density_summary.csv');
sessionDistributionPath = fullfile(outRoot, 'feature_session_distribution_summary.csv');
correlationPath = fullfile(outRoot, 'feature_correlation_summary.csv');
arenaShiftPath = fullfile(outRoot, 'feature_arena_shift_summary.csv');
figureManifestPath = fullfile(outRoot, 'feature_qc_figure_manifest.csv');

writetable(local_parameter_audit_table(params), parameterAuditPath);
writetable(featureMeta, featureDictionaryPath);

Q = readtable(qcPath, 'TextType', 'string');
motifMask = select_analysis_cohort(Q, "motif_discovery");
targetMask = motifMask & Q.preprocess_success == 1 & ...
    Q.qc_pass_motif_discovery == 1 & Q.effective_n_animals == 2;
targetIdx = find(targetMask);
assert(~isempty(targetIdx), 'run_05_extract_motif_dyad_features:NoSessions', ...
    'No QC-pass motif-discovery dyad sessions are available.');

selectedIdx = local_select_sessions(Q, targetIdx, params);
assert(~isempty(selectedIdx), 'run_05_extract_motif_dyad_features:NoSelectedSessions', ...
    'Feature extraction selected zero sessions.');

fprintf('run_05_extract_motif_dyad_features\n');
fprintf('Repo root: %s\n', repoRoot);
fprintf('Output root: %s\n', outRoot);
fprintf('QC table: %s\n', qcPath);
fprintf('Parameter audit: %s\n', parameterAuditPath);
fprintf('Feature dictionary: %s\n', featureDictionaryPath);
fprintf('Run mode: %s\n', params.run_mode);
fprintf('QC-pass motif dyad sessions available: %d\n', numel(targetIdx));
fprintf('Sessions selected for this run: %d\n', numel(selectedIdx));
fprintf('Skip existing: %d\n', params.skip_existing);

sessionAudit = local_initialize_session_audit(Q, selectedIdx, params);
arenaLevels = local_arena_levels(Q, selectedIdx);
featureValues = cell(numel(arenaLevels), numel(featureNames));
featureTotalFrames = zeros(numel(featureNames), 1);
featureValidFrames = zeros(numel(featureNames), 1);
featureMissingAll = zeros(numel(featureNames), 1);
featureMissingValid = zeros(numel(featureNames), 1);
featureMatrixSamples = cell(numel(arenaLevels), 1);
sessionDistributionSummary = table();

for k = 1:numel(selectedIdx)
    rowIdx = selectedIdx(k);
    row = Q(rowIdx,:);
    outFile = fullfile(sessionOutDir, sprintf('session_%04d_dyad_features.mat', row.raw_index));

    fprintf('\n[%d/%d] raw_index=%d session=%s arena=%s qc_set=%s\n', ...
        k, numel(selectedIdx), row.raw_index, row.session_id, ...
        local_arena_label(row), sessionAudit.qc_set(k));

    tStart = tic;

    try
        if params.skip_existing && isfile(outFile)
            S = load(outFile, 'dyad', 'status');
            assert(isfield(S, 'dyad') && string(S.status) == "success", ...
                'Existing feature file is not a successful dyad output.');
            dyad = S.dyad;
            local_assert_unit_params(dyad.params.fps, dyad.params.pixelSizeMM, params, ...
                sprintf('existing feature file %s', outFile));
            runtime_sec = 0;
            fprintf('  loaded existing successful output\n');
        else
            preprocFile = local_resolve_preprocess_output_file(Q, rowIdx, repoRoot, cfg);
            S = load(preprocFile, 'sessionPreproc', 'manifestRow', 'animalQC');
            assert(isfield(S, 'sessionPreproc'), 'Missing sessionPreproc in %s', preprocFile);
            sessionPreproc = S.sessionPreproc;
            tracks = sessionPreproc.clean.tracks;
            assert(ndims(tracks) == 4 && size(tracks,4) == 2, ...
                'Expected cleaned dyad tracks [T x N x xy x 2].');

            [fps, pixelSizeMM] = local_unit_params_from_preprocessing(sessionPreproc, params, preprocFile);
            badframes = sessionPreproc.qc.badframes;

            dyad = compute_dyad_features(tracks, fps, nodeMap, ...
                'badframes', badframes, ...
                'pixelSizeMM', pixelSizeMM, ...
                'contactThresholdMM', params.contact_threshold_mm, ...
                'closeThresholdMM', params.close_threshold_mm, ...
                'smoothSpanFrames', params.smooth_span_frames);
            runtime_sec = toc(tStart);

            manifestRow = row;
            if isfield(S, 'animalQC')
                animalQC = S.animalQC;
            else
                animalQC = struct();
            end
            status = "success";
            save(outFile, 'dyad', 'manifestRow', 'animalQC', 'nodeMap', ...
                'partNames', 'nodeMetadata', 'featureMeta', 'params', ...
                'status', 'runtime_sec', '-v7.3');
        end

        sessionAudit.status(k) = "success";
        sessionAudit.runtime_sec(k) = runtime_sec;
        sessionAudit.feature_output_file(k) = local_repo_relative_path(repoRoot, outFile);
        sessionAudit.n_frames(k) = dyad.maskAudit.nFrames;
        sessionAudit.fps(k) = dyad.params.fps;
        sessionAudit.pixel_size_mm(k) = dyad.params.pixelSizeMM;
        sessionAudit.n_features(k) = numel(dyad.featureNames);
        sessionAudit.valid_frame_fraction(k) = dyad.maskAudit.validFrameFraction;
        sessionAudit.badframe_mask_fraction(k) = dyad.maskAudit.badframeFraction;
        sessionAudit.core_nan_fraction(k) = dyad.maskAudit.coreNanFraction;
        sessionAudit.core_nan_only_fraction(k) = dyad.maskAudit.coreNanOnlyFraction;
        sessionAudit.feature_matrix_nan_fraction(k) = dyad.maskAudit.featureMatrixNanFraction;

        [featureTotalFrames, featureValidFrames, featureMissingAll, ...
                featureMissingValid, featureValues, featureMatrixSamples] = ...
            local_accumulate_feature_qc(dyad, row, arenaLevels, featureNames, ...
            featureTotalFrames, featureValidFrames, featureMissingAll, ...
            featureMissingValid, featureValues, featureMatrixSamples, params);
        sessionDistributionSummary = [sessionDistributionSummary; ...
            local_session_distribution_summary(dyad, row, sessionAudit(k,:), featureMeta)]; %#ok<AGROW>

        fprintf('  success runtime=%.1fs valid_frame_fraction=%.4f feature_nan_fraction=%.4f\n', ...
            runtime_sec, dyad.maskAudit.validFrameFraction, dyad.maskAudit.featureMatrixNanFraction);
    catch ME
        runtime_sec = toc(tStart);
        status = "failed";
        error_message = string(ME.message);
        if params.skip_existing && isfile(outFile)
            failureFile = fullfile(logDir, sprintf('session_%04d_feature_failure.mat', row.raw_index));
            save(failureFile, 'row', 'nodeMap', 'nodeMetadata', 'featureMeta', ...
                'params', 'status', 'runtime_sec', 'error_message', '-v7.3');
        else
            save(outFile, 'row', 'nodeMap', 'nodeMetadata', 'featureMeta', ...
                'params', 'status', 'runtime_sec', 'error_message', '-v7.3');
        end
        sessionAudit.status(k) = "failed";
        sessionAudit.runtime_sec(k) = runtime_sec;
        sessionAudit.error_message(k) = string(ME.message);
        sessionAudit.feature_output_file(k) = local_repo_relative_path(repoRoot, outFile);
        fprintf('  FAILED runtime=%.1fs error=%s\n', runtime_sec, ME.message);
    end

    writetable(sessionAudit, sessionAuditPath);
end

writetable(sessionAudit, sessionAuditPath);
missingnessSummary = local_missingness_summary(featureMeta, sessionAudit, ...
    featureTotalFrames, featureValidFrames, featureMissingAll, featureMissingValid);
distributionSummary = local_distribution_summary(featureMeta, arenaLevels, featureValues);
distributionDensitySummary = local_distribution_density_summary(featureMeta, arenaLevels, ...
    featureValues, params);
correlationSummary = local_correlation_summary(featureMeta, arenaLevels, ...
    featureMatrixSamples, params);
arenaShiftSummary = local_arena_shift_summary(featureMeta, arenaLevels, featureValues);
writetable(missingnessSummary, missingnessPath);
writetable(distributionSummary, distributionPath);
writetable(distributionDensitySummary, distributionDensityPath);
writetable(sessionDistributionSummary, sessionDistributionPath);
writetable(correlationSummary, correlationPath);
writetable(arenaShiftSummary, arenaShiftPath);

figureManifest = local_make_figures(missingnessPath, distributionDensityPath, ...
    sessionDistributionPath, sessionAuditPath, correlationPath, arenaShiftPath, ...
    figDir, cfg, params, repoRoot);
writetable(figureManifest, figureManifestPath);

fprintf('\nWrote feature dictionary: %s\n', featureDictionaryPath);
fprintf('Wrote parameter audit: %s\n', parameterAuditPath);
fprintf('Wrote session audit: %s\n', sessionAuditPath);
fprintf('Wrote missingness/validity summary: %s\n', missingnessPath);
fprintf('Wrote distribution summary: %s\n', distributionPath);
fprintf('Wrote distribution density summary: %s\n', distributionDensityPath);
fprintf('Wrote session distribution summary: %s\n', sessionDistributionPath);
fprintf('Wrote correlation summary: %s\n', correlationPath);
fprintf('Wrote arena shift summary: %s\n', arenaShiftPath);
fprintf('Wrote figure manifest: %s\n', figureManifestPath);
disp(groupsummary(sessionAudit, 'status'));

nFailed = nnz(sessionAudit.status == "failed");
assert(nFailed == 0, 'run_05_extract_motif_dyad_features:Failures', ...
    'Feature extraction finished with %d failed sessions. See %s.', nFailed, sessionAuditPath);

outputs = struct();
outputs.output_root = outRoot;
outputs.parameter_audit_path = parameterAuditPath;
outputs.feature_dictionary_path = featureDictionaryPath;
outputs.session_audit_path = sessionAuditPath;
outputs.missingness_path = missingnessPath;
outputs.distribution_path = distributionPath;
outputs.distribution_density_path = distributionDensityPath;
outputs.session_distribution_path = sessionDistributionPath;
outputs.correlation_path = correlationPath;
outputs.arena_shift_path = arenaShiftPath;
outputs.figure_manifest_path = figureManifestPath;
outputs.n_sessions_selected = numel(selectedIdx);
outputs.n_sessions_failed = nFailed;
end

function T = local_parameter_audit_table(params)
T = params.config_table;
effectiveValue = strings(height(T), 1);
for i = 1:height(T)
    fieldName = matlab.lang.makeValidName(T.parameter(i));
    if isfield(params, fieldName)
        effectiveValue(i) = local_value_to_string(params.(fieldName));
    else
        effectiveValue(i) = "";
    end
end
T.effective_value = effectiveValue;
T.config_path = repmat(string(params.config_path), height(T), 1);
end

function value = local_value_to_string(x)
if islogical(x)
    value = string(x);
elseif isnumeric(x)
    if isscalar(x) && isnan(x)
        value = "NaN";
    elseif isscalar(x)
        value = string(sprintf('%.15g', x));
    else
        value = strjoin(string(compose('%.15g', x(:)')), ";");
    end
elseif isstring(x) || ischar(x)
    value = strjoin(string(x), ";");
else
    value = string(evalc('disp(x)'));
    value = strtrim(value);
end
end

function [fps, pixelSizeMM] = local_unit_params_from_preprocessing(sessionPreproc, params, sourceFile)
assert(isfield(sessionPreproc, 'params') && isfield(sessionPreproc.params, 'data'), ...
    'run_motif_dyad_feature_extraction:MissingUnitParams', ...
    'Missing preprocessing unit parameters in %s.', sourceFile);
assert(isfield(sessionPreproc.params.data, 'fps') && ...
    isfield(sessionPreproc.params.data, 'pixel_size_mm'), ...
    'run_motif_dyad_feature_extraction:MissingUnitParams', ...
    'Missing fps or pixel_size_mm in %s.', sourceFile);

local_assert_unit_params(sessionPreproc.params.data.fps, ...
    sessionPreproc.params.data.pixel_size_mm, params, sourceFile);

% Use config values after the assertion so the feature-stage unit contract
% is controlled from config/motif_dyad_feature_extraction_config.csv.
fps = params.fps;
pixelSizeMM = params.pixel_size_mm;
end

function local_assert_unit_params(fps, pixelSizeMM, params, sourceLabel)
tol = params.unit_parameter_tolerance;
assert(abs(double(fps) - params.fps) <= tol, ...
    'run_motif_dyad_feature_extraction:FpsMismatch', ...
    'fps in %s is %.15g but config expects %.15g.', ...
    sourceLabel, double(fps), params.fps);
assert(abs(double(pixelSizeMM) - params.pixel_size_mm) <= tol, ...
    'run_motif_dyad_feature_extraction:PixelSizeMismatch', ...
    'pixel_size_mm in %s is %.15g but config expects %.15g.', ...
    sourceLabel, double(pixelSizeMM), params.pixel_size_mm);
end

function selectedIdx = local_select_sessions(Q, targetIdx, params)
targetIdx = targetIdx(:);
if params.run_mode == "full"
    selectedIdx = targetIdx;
else
    arenas = local_arena_levels(Q, targetIdx);
    selectedIdx = [];
    for a = 1:numel(arenas)
        inArena = targetIdx(local_row_arena_values(Q, targetIdx) == arenas(a));
        clean = inArena(Q.qc_warn_motif_discovery(inArena) == 0);
        review = inArena(Q.qc_warn_motif_discovery(inArena) == 1);
        take = unique([clean(:); review(:)], 'stable');
        take = sortrows(Q(take, {'raw_index'}), 'raw_index');
        [~, loc] = ismember(take.raw_index, Q.raw_index);
        nTake = min(numel(loc), params.smoke_sessions_per_arena);
        selectedIdx = [selectedIdx; loc(1:nTake)]; %#ok<AGROW>
    end
    selectedIdx = unique(selectedIdx, 'stable');
end

if isfinite(params.max_sessions)
    selectedIdx = selectedIdx(1:min(numel(selectedIdx), params.max_sessions));
end
end

function sessionAudit = local_initialize_session_audit(Q, selectedIdx, params)
n = numel(selectedIdx);
sessionAudit = table();
sessionAudit.selection_rank = (1:n)';
sessionAudit.run_mode = repmat(params.run_mode, n, 1);
sessionAudit.selection_rule = repmat(params.selection_rule, n, 1);
sessionAudit.raw_index = Q.raw_index(selectedIdx);
sessionAudit.session_file = Q.session_file(selectedIdx);
sessionAudit.session_id = Q.session_id(selectedIdx);
sessionAudit.arena = local_table_string(Q, 'arena', selectedIdx, "");
sessionAudit.arena_label = local_table_string(Q, 'arena_label', selectedIdx, "");
sessionAudit.cohort_id = local_table_string(Q, 'cohort_id', selectedIdx, "");
sessionAudit.cohort_label = local_table_string(Q, 'cohort_label', selectedIdx, "");
sessionAudit.condition_id = local_table_string(Q, 'condition_id', selectedIdx, "");
sessionAudit.condition_label = local_table_string(Q, 'condition_label', selectedIdx, "");
sessionAudit.condition_role = repmat("provenance_only_not_used_for_features", n, 1);
sessionAudit.include_motif_discovery = Q.include_motif_discovery(selectedIdx);
sessionAudit.qc_pass_motif_discovery = Q.qc_pass_motif_discovery(selectedIdx);
sessionAudit.qc_requires_manual_review = Q.qc_requires_manual_review(selectedIdx);
sessionAudit.qc_warn_motif_discovery = Q.qc_warn_motif_discovery(selectedIdx);
sessionAudit.qc_set = repmat("primary_qc_pass_clean", n, 1);
sessionAudit.qc_set(sessionAudit.qc_warn_motif_discovery == 1) = "primary_qc_pass_review_band";
sessionAudit.badframe_fraction = Q.badframe_fraction(selectedIdx);
sessionAudit.preprocess_output_file = Q.preprocess_output_file(selectedIdx);
sessionAudit.feature_output_file = strings(n, 1);
sessionAudit.contact_threshold_mm = repmat(params.contact_threshold_mm, n, 1);
sessionAudit.close_threshold_mm = repmat(params.close_threshold_mm, n, 1);
sessionAudit.smooth_span_frames = repmat(params.smooth_span_frames, n, 1);
sessionAudit.status = repmat("pending", n, 1);
sessionAudit.runtime_sec = nan(n, 1);
sessionAudit.n_frames = nan(n, 1);
sessionAudit.fps = nan(n, 1);
sessionAudit.pixel_size_mm = nan(n, 1);
sessionAudit.n_features = nan(n, 1);
sessionAudit.valid_frame_fraction = nan(n, 1);
sessionAudit.badframe_mask_fraction = nan(n, 1);
sessionAudit.core_nan_fraction = nan(n, 1);
sessionAudit.core_nan_only_fraction = nan(n, 1);
sessionAudit.feature_matrix_nan_fraction = nan(n, 1);
sessionAudit.error_message = strings(n, 1);
end

function values = local_table_string(T, varName, rows, defaultValue)
if ismember(varName, T.Properties.VariableNames)
    values = string(T.(varName)(rows));
else
    values = repmat(string(defaultValue), numel(rows), 1);
end
end

function arenaLevels = local_arena_levels(Q, rows)
arenaValues = local_row_arena_values(Q, rows);
arenaLevels = unique(arenaValues, 'stable');
arenaLevels = arenaLevels(arenaLevels ~= "");
end

function arenaValues = local_row_arena_values(Q, rows)
if ismember('arena_label', Q.Properties.VariableNames)
    arenaValues = string(Q.arena_label(rows));
else
    arenaValues = string(Q.arena(rows));
end
end

function label = local_arena_label(row)
if ismember('arena_label', row.Properties.VariableNames)
    label = string(row.arena_label);
else
    label = string(row.arena);
end
end

function pathOut = local_resolve_preprocess_output_file(Q, rowIdx, repoRoot, cfg)
candidate = string(Q.preprocess_output_file(rowIdx));
if strlength(candidate) > 0
    candidateAbs = resolve_repo_path(repoRoot, candidate);
    if isfile(candidateAbs)
        pathOut = candidateAbs;
        return
    end
end
preprocRoot = resolve_repo_path(repoRoot, cfg.paths.preprocessed_output_dir);
sessionDir = fullfile(preprocRoot, char(cfg.paths.preprocess_session_subdir));
pathOut = fullfile(sessionDir, sprintf('session_%04d_preprocessed.mat', Q.raw_index(rowIdx)));
assert(isfile(pathOut), 'run_05_extract_motif_dyad_features:MissingPreprocessedMat', ...
    'Missing preprocessed session MAT: %s', pathOut);
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

function [totalFrames, validFrames, missingAll, missingValid, featureValues, featureMatrixSamples] = ...
    local_accumulate_feature_qc(dyad, row, arenaLevels, featureNames, ...
    totalFrames, validFrames, missingAll, missingValid, featureValues, featureMatrixSamples, params)
arena = local_arena_label(row);
arenaIdx = find(arenaLevels == arena, 1);
assert(~isempty(arenaIdx), 'Unknown arena label: %s', arena);

valid = logical(dyad.frameMask(:));
for f = 1:numel(featureNames)
    x = dyad.X(:, f);
    totalFrames(f) = totalFrames(f) + numel(x);
    validFrames(f) = validFrames(f) + nnz(valid);
    missingAll(f) = missingAll(f) + nnz(~isfinite(x));
    missingValid(f) = missingValid(f) + nnz(valid & ~isfinite(x));
    xv = x(isfinite(x));
    featureValues{arenaIdx, f} = [featureValues{arenaIdx, f}; xv];
end

Xvalid = dyad.X(valid, :);
Xvalid = Xvalid(any(isfinite(Xvalid), 2), :);
if ~isempty(Xvalid)
    featureMatrixSamples{arenaIdx} = [featureMatrixSamples{arenaIdx}; Xvalid];
    featureMatrixSamples{arenaIdx} = local_downsample_rows( ...
        featureMatrixSamples{arenaIdx}, params.max_correlation_frames_per_arena);
end
end

function S = local_missingness_summary(featureMeta, sessionAudit, totalFrames, validFrames, missingAll, missingValid)
S = featureMeta(:, {'FeatureIndex','Name','Family','FamilyLabel', ...
    'FeatureFamilyRole','FamilyDefinition','Unit','UnitLabel', ...
    'IsDirected','Directionality','IsCircular','IsBoolean','TransformHint', ...
    'BiologicalInterpretation','ClusteringCandidate','FeatureLayerRole', ...
    'WindowSummaryStats'});
nSuccess = nnz(sessionAudit.status == "success");
S.n_sessions = repmat(nSuccess, height(S), 1);
S.n_frames_total = totalFrames;
S.n_valid_frames = validFrames;
S.n_nonfinite_all_frames = missingAll;
S.missing_fraction_all_frames = missingAll ./ max(totalFrames, 1);
S.n_nonfinite_valid_frames = missingValid;
S.missing_fraction_valid_frames = missingValid ./ max(validFrames, 1);
S.valid_frame_fraction_overall = validFrames ./ max(totalFrames, 1);
end

function D = local_distribution_summary(featureMeta, arenaLevels, featureValues)
rows = struct([]);
for a = 1:numel(arenaLevels)
    for f = 1:height(featureMeta)
        x = featureValues{a, f};
        x = x(isfinite(x));
        row = struct();
        row.arena_label = arenaLevels(a);
        row.feature_index = featureMeta.FeatureIndex(f);
        row.feature_name = string(featureMeta.Name(f));
        row.family = string(featureMeta.Family(f));
        row.family_label = string(featureMeta.FamilyLabel(f));
        row.feature_family_role = string(featureMeta.FeatureFamilyRole(f));
        row.family_definition = string(featureMeta.FamilyDefinition(f));
        row.unit = string(featureMeta.Unit(f));
        row.feature_type = local_feature_type_from_meta(featureMeta, f);
        row.transform_hint = string(featureMeta.TransformHint(f));
        row.is_circular = logical(featureMeta.IsCircular(f));
        row.is_boolean = logical(featureMeta.IsBoolean(f));
        row.n_values = numel(x);
        if isempty(x)
            row.mean = NaN;
            row.sd = NaN;
            row.min = NaN;
            row.q05 = NaN;
            row.q25 = NaN;
            row.median = NaN;
            row.q75 = NaN;
            row.q95 = NaN;
            row.max = NaN;
            row.circular_mean_deg = NaN;
            row.circular_resultant = NaN;
        else
            row.mean = mean(x);
            row.sd = std(x, 0);
            row.min = min(x);
            q = prctile(x, [5 25 50 75 95]);
            row.q05 = q(1);
            row.q25 = q(2);
            row.median = q(3);
            row.q75 = q(4);
            row.q95 = q(5);
            row.max = max(x);
            if row.is_circular
                r = deg2rad(x);
                c = mean(cos(r));
                s = mean(sin(r));
                row.circular_mean_deg = rad2deg(atan2(s, c));
                row.circular_resultant = hypot(c, s);
            else
                row.circular_mean_deg = NaN;
                row.circular_resultant = NaN;
            end
        end
        rows = [rows; row]; %#ok<AGROW>
    end
end
D = struct2table(rows);
end

function S = local_session_distribution_summary(dyad, row, sessionAuditRow, featureMeta)
rows = struct([]);
featureNames = string(dyad.featureNames(:));
validFrameMask = dyad.frameMask(:);
for f = 1:numel(featureNames)
    metaIdx = find(featureMeta.Name == featureNames(f), 1);
    assert(~isempty(metaIdx), 'run_motif_dyad_feature_extraction:UnknownFeature', ...
        'Feature "%s" is missing from canonical metadata.', featureNames(f));
    x = dyad.X(:, f);
    x = x(validFrameMask & isfinite(x));
    stats = local_distribution_row_stats(x, featureMeta.IsCircular(metaIdx) == 1);

    rowOut = struct();
    rowOut.selection_rank = sessionAuditRow.selection_rank;
    rowOut.raw_index = row.raw_index;
    rowOut.session_id = string(row.session_id);
    rowOut.arena = string(row.arena);
    rowOut.arena_label = local_arena_label(row);
    rowOut.qc_set = string(sessionAuditRow.qc_set);
    rowOut.feature_index = featureMeta.FeatureIndex(metaIdx);
    rowOut.feature_name = featureMeta.Name(metaIdx);
    rowOut.family = featureMeta.Family(metaIdx);
    rowOut.family_label = featureMeta.FamilyLabel(metaIdx);
    rowOut.feature_family_role = featureMeta.FeatureFamilyRole(metaIdx);
    rowOut.family_definition = featureMeta.FamilyDefinition(metaIdx);
    rowOut.unit = featureMeta.Unit(metaIdx);
    rowOut.feature_type = local_feature_type_from_meta(featureMeta, metaIdx);
    rowOut.transform_hint = featureMeta.TransformHint(metaIdx);
    rowOut.is_circular = featureMeta.IsCircular(metaIdx);
    rowOut.is_boolean = featureMeta.IsBoolean(metaIdx);
    rowOut.n_valid_frames = nnz(validFrameMask);
    rowOut.n_valid_values = numel(x);
    if isempty(x)
        rowOut.mean = NaN;
        rowOut.sd = NaN;
        rowOut.min = NaN;
        rowOut.max = NaN;
    else
        rowOut.mean = mean(x, 'omitnan');
        rowOut.sd = std(x, 'omitnan');
        rowOut.min = min(x, [], 'omitnan');
        rowOut.max = max(x, [], 'omitnan');
    end
    rowOut.q05 = stats.q05;
    rowOut.q25 = stats.q25;
    rowOut.median = stats.median;
    rowOut.q75 = stats.q75;
    rowOut.q95 = stats.q95;
    rowOut.prevalence_true = NaN;
    if featureMeta.IsBoolean(metaIdx) == 1 && ~isempty(x)
        rowOut.prevalence_true = mean(x ~= 0, 'omitnan');
    end
    rowOut.circular_mean_deg = stats.circular_mean_deg;
    rowOut.circular_resultant = stats.circular_resultant;
    rows = [rows; rowOut]; %#ok<AGROW>
end
S = struct2table(rows);
end

function B = local_distribution_density_summary(featureMeta, arenaLevels, featureValues, params)
rows = struct([]);
for f = 1:height(featureMeta)
    allX = [];
    for a = 1:numel(arenaLevels)
        x = featureValues{a, f};
        allX = [allX; x(isfinite(x))]; %#ok<AGROW>
    end
    if isempty(allX)
        continue
    end

    isCircular = logical(featureMeta.IsCircular(f));
    isBoolean = logical(featureMeta.IsBoolean(f));
    if isBoolean
        edges = [-0.5 0.5 1.5];
        centers = [0 1];
        displayMin = -0.5;
        displayMax = 1.5;
        binMode = "boolean_probability";
    elseif isCircular
        displayMin = -180;
        displayMax = 180;
        edges = linspace(displayMin, displayMax, params.circular_density_bins + 1);
        centers = edges(1:end-1) + diff(edges)./2;
        binMode = "circular_degree_density";
    else
        q = prctile(allX, [0.5 99.5]);
        displayMin = q(1);
        displayMax = q(2);
        if ~isfinite(displayMin) || ~isfinite(displayMax) || displayMin >= displayMax
            displayMin = min(allX);
            displayMax = max(allX);
        end
        if displayMin == displayMax
            displayMin = displayMin - 0.5;
            displayMax = displayMax + 0.5;
        end
        edges = linspace(displayMin, displayMax, params.distribution_density_bins + 1);
        centers = edges(1:end-1) + diff(edges)./2;
        binMode = "winsorized_common_axis_density";
    end

    for a = 1:numel(arenaLevels)
        x = featureValues{a, f};
        x = x(isfinite(x));
        n = numel(x);
        nBelowDisplay = 0;
        nAboveDisplay = 0;
        if isBoolean
            counts = [nnz(x == 0), nnz(x == 1)];
            nDensity = sum(counts);
        elseif isCircular
            xPlot = mod(x + 180, 360) - 180;
            xPlot(xPlot == 180) = -180;
            counts = histcounts(xPlot, edges);
            nDensity = sum(counts);
        else
            nBelowDisplay = nnz(x < displayMin);
            nAboveDisplay = nnz(x > displayMax);
            xPlot = x(x >= displayMin & x <= displayMax);
            counts = histcounts(xPlot, edges);
            nDensity = sum(counts);
        end

        stats = local_distribution_row_stats(x, isCircular);
        for b = 1:numel(centers)
            binWidth = edges(b+1) - edges(b);
            row = struct();
            row.feature_index = featureMeta.FeatureIndex(f);
            row.feature_name = string(featureMeta.Name(f));
            row.family = string(featureMeta.Family(f));
            row.family_label = string(featureMeta.FamilyLabel(f));
            row.feature_family_role = string(featureMeta.FeatureFamilyRole(f));
            row.family_definition = string(featureMeta.FamilyDefinition(f));
            row.unit = string(featureMeta.Unit(f));
            row.feature_type = local_feature_type_from_meta(featureMeta, f);
            row.transform_hint = string(featureMeta.TransformHint(f));
            row.is_circular = isCircular;
            row.is_boolean = isBoolean;
            row.arena_label = arenaLevels(a);
            row.bin_mode = binMode;
            row.bin_index = b;
            row.bin_left = edges(b);
            row.bin_right = edges(b+1);
            row.bin_center = centers(b);
            row.bin_width = binWidth;
            row.count = counts(b);
            row.n_values = n;
            row.n_density_values = nDensity;
            row.n_below_display = nBelowDisplay;
            row.n_above_display = nAboveDisplay;
            row.probability = counts(b) ./ max(n, 1);
            row.density = counts(b) ./ max(nDensity .* binWidth, eps);
            row.display_min = displayMin;
            row.display_max = displayMax;
            row.q05 = stats.q05;
            row.q25 = stats.q25;
            row.median = stats.median;
            row.q75 = stats.q75;
            row.q95 = stats.q95;
            row.circular_mean_deg = stats.circular_mean_deg;
            row.circular_resultant = stats.circular_resultant;
            rows = [rows; row]; %#ok<AGROW>
        end
    end
end
B = struct2table(rows);
end

function stats = local_distribution_row_stats(x, isCircular)
stats = struct('q05', NaN, 'q25', NaN, 'median', NaN, 'q75', NaN, ...
    'q95', NaN, 'circular_mean_deg', NaN, 'circular_resultant', NaN);
if isempty(x)
    return
end
q = prctile(x, [5 25 50 75 95]);
stats.q05 = q(1);
stats.q25 = q(2);
stats.median = q(3);
stats.q75 = q(4);
stats.q95 = q(5);
if isCircular
    r = deg2rad(x);
    c = mean(cos(r));
    s = mean(sin(r));
    stats.circular_mean_deg = rad2deg(atan2(s, c));
    stats.circular_resultant = hypot(c, s);
end
end

function C = local_correlation_summary(featureMeta, arenaLevels, featureMatrixSamples, params)
rows = struct([]);
scopeLabels = ["all"; arenaLevels(:)];
for s = 1:numel(scopeLabels)
    if scopeLabels(s) == "all"
        nonEmpty = ~cellfun(@isempty, featureMatrixSamples);
        if any(nonEmpty)
            X = vertcat(featureMatrixSamples{nonEmpty});
        else
            X = zeros(0, height(featureMeta));
        end
    else
        arenaIdx = find(arenaLevels == scopeLabels(s), 1);
        X = featureMatrixSamples{arenaIdx};
    end
    X = local_downsample_rows(X, params.max_correlation_frames_per_arena);
    [Y, channelMeta] = local_correlation_channels(X, featureMeta);
    for i = 1:height(channelMeta)
        for j = i:height(channelMeta)
            ok = isfinite(Y(:, i)) & isfinite(Y(:, j));
            nPair = nnz(ok);
            rho = NaN;
            if nPair >= 3 && numel(unique(Y(ok, i))) > 1 && numel(unique(Y(ok, j))) > 1
                rho = corr(Y(ok, i), Y(ok, j), 'Type','Spearman', 'Rows','complete');
            end
            row = struct();
            row.arena_label = scopeLabels(s);
            row.diagnostic_sample_rows = size(Y, 1);
            row.correlation_method = "spearman_pairwise_complete_on_diagnostic_channels";
            row.channel_i_index = channelMeta.ChannelIndex(i);
            row.channel_i = channelMeta.ChannelName(i);
            row.source_feature_i_index = channelMeta.FeatureIndex(i);
            row.source_feature_i = channelMeta.FeatureName(i);
            row.family_i = channelMeta.Family(i);
            row.unit_i = channelMeta.Unit(i);
            row.transform_i = channelMeta.Transform(i);
            row.channel_j_index = channelMeta.ChannelIndex(j);
            row.channel_j = channelMeta.ChannelName(j);
            row.source_feature_j_index = channelMeta.FeatureIndex(j);
            row.source_feature_j = channelMeta.FeatureName(j);
            row.family_j = channelMeta.Family(j);
            row.unit_j = channelMeta.Unit(j);
            row.transform_j = channelMeta.Transform(j);
            row.n_pairwise = nPair;
            row.spearman_r = rho;
            rows = [rows; row]; %#ok<AGROW>
        end
    end
end
C = struct2table(rows);
end

function [Y, channelMeta] = local_correlation_channels(X, featureMeta)
Y = zeros(size(X, 1), 0);
rows = struct([]);
for f = 1:height(featureMeta)
    x = X(:, f);
    if featureMeta.IsCircular(f)
        Y(:, end+1) = sind(x); %#ok<AGROW>
        rows = local_add_channel_row(rows, size(Y,2), featureMeta, f, ...
            string(featureMeta.Name(f)) + "__sin", "sind_circular_projection");
        Y(:, end+1) = cosd(x); %#ok<AGROW>
        rows = local_add_channel_row(rows, size(Y,2), featureMeta, f, ...
            string(featureMeta.Name(f)) + "__cos", "cosd_circular_projection");
    else
        transformHint = string(featureMeta.TransformHint(f));
        if transformHint == "log1p"
            y = log1p(max(x, 0));
            transform = "log1p_nonnegative";
        elseif featureMeta.IsBoolean(f)
            y = x;
            transform = "identity_binary";
        else
            y = x;
            transform = "identity_rank_diagnostic";
        end
        Y(:, end+1) = y; %#ok<AGROW>
        rows = local_add_channel_row(rows, size(Y,2), featureMeta, f, ...
            string(featureMeta.Name(f)), transform);
    end
end
channelMeta = struct2table(rows);
end

function rows = local_add_channel_row(rows, channelIndex, featureMeta, featureIndex, channelName, transform)
row = struct();
row.ChannelIndex = channelIndex;
row.ChannelName = string(channelName);
row.FeatureIndex = featureMeta.FeatureIndex(featureIndex);
row.FeatureName = string(featureMeta.Name(featureIndex));
row.Family = string(featureMeta.Family(featureIndex));
row.FamilyLabel = string(featureMeta.FamilyLabel(featureIndex));
row.FeatureFamilyRole = string(featureMeta.FeatureFamilyRole(featureIndex));
row.FamilyDefinition = string(featureMeta.FamilyDefinition(featureIndex));
row.Unit = string(featureMeta.Unit(featureIndex));
row.FeatureType = local_feature_type_from_meta(featureMeta, featureIndex);
row.Transform = string(transform);
rows = [rows; row];
end

function A = local_arena_shift_summary(featureMeta, arenaLevels, featureValues)
rows = struct([]);
if numel(arenaLevels) < 2
    A = table();
    return
end

for a = 1:numel(arenaLevels)-1
    for b = a+1:numel(arenaLevels)
        for f = 1:height(featureMeta)
            xA = featureValues{a, f};
            xB = featureValues{b, f};
            xA = xA(isfinite(xA));
            xB = xB(isfinite(xB));

            row = struct();
            row.comparison = string(arenaLevels(b)) + "_minus_" + string(arenaLevels(a));
            row.arena_a = arenaLevels(a);
            row.arena_b = arenaLevels(b);
            row.feature_index = featureMeta.FeatureIndex(f);
            row.feature_name = string(featureMeta.Name(f));
            row.family = string(featureMeta.Family(f));
            row.family_label = string(featureMeta.FamilyLabel(f));
            row.feature_family_role = string(featureMeta.FeatureFamilyRole(f));
            row.family_definition = string(featureMeta.FamilyDefinition(f));
            row.unit = string(featureMeta.Unit(f));
            row.feature_type = local_feature_type_from_meta(featureMeta, f);
            row.transform_hint = string(featureMeta.TransformHint(f));
            row.is_circular = logical(featureMeta.IsCircular(f));
            row.is_boolean = logical(featureMeta.IsBoolean(f));
            row.n_a = numel(xA);
            row.n_b = numel(xB);
            [statsA, circA] = local_distribution_stats(xA, row.is_circular);
            [statsB, circB] = local_distribution_stats(xB, row.is_circular);
            row.mean_a = statsA.mean;
            row.mean_b = statsB.mean;
            row.median_a = statsA.median;
            row.median_b = statsB.median;
            row.q25_a = statsA.q25;
            row.q75_a = statsA.q75;
            row.q25_b = statsB.q25;
            row.q75_b = statsB.q75;
            row.raw_median_difference = statsB.median - statsA.median;
            row.mean_difference = statsB.mean - statsA.mean;
            row.prevalence_difference = NaN;
            row.circular_mean_a_deg = circA.mean_deg;
            row.circular_mean_b_deg = circB.mean_deg;
            row.circular_mean_difference_deg = NaN;
            row.robust_median_difference_iqr = NaN;

            pooledIqr = mean([statsA.q75 - statsA.q25, statsB.q75 - statsB.q25], 'omitnan');
            if row.is_circular
                row.circular_mean_difference_deg = local_angle_diff_deg( ...
                    circA.mean_deg, circB.mean_deg);
                row.metric_name = "circular_mean_difference_deg";
                row.metric_value = row.circular_mean_difference_deg;
            elseif row.is_boolean
                row.prevalence_difference = row.mean_difference;
                row.metric_name = "prevalence_difference";
                row.metric_value = row.prevalence_difference;
            else
                if isfinite(pooledIqr) && pooledIqr > 0
                    row.robust_median_difference_iqr = row.raw_median_difference ./ pooledIqr;
                end
                row.metric_name = "robust_median_difference_iqr";
                row.metric_value = row.robust_median_difference_iqr;
            end
            rows = [rows; row]; %#ok<AGROW>
        end
    end
end
A = struct2table(rows);
end

function [stats, circ] = local_distribution_stats(x, isCircular)
stats = struct('mean', NaN, 'median', NaN, 'q25', NaN, 'q75', NaN);
circ = struct('mean_deg', NaN, 'resultant', NaN);
if isempty(x)
    return
end
stats.mean = mean(x);
q = prctile(x, [25 50 75]);
stats.q25 = q(1);
stats.median = q(2);
stats.q75 = q(3);
if isCircular
    r = deg2rad(x);
    c = mean(cos(r));
    s = mean(sin(r));
    circ.mean_deg = rad2deg(atan2(s, c));
    circ.resultant = hypot(c, s);
end
end

function d = local_angle_diff_deg(a, b)
if ~isfinite(a) || ~isfinite(b)
    d = NaN;
else
    d = mod((b - a) + 180, 360) - 180;
end
end

function X = local_downsample_rows(X, maxRows)
if isempty(X) || ~isfinite(maxRows) || size(X, 1) <= maxRows
    return
end
idx = unique(round(linspace(1, size(X, 1), maxRows)));
X = X(idx, :);
end

function figureManifest = local_make_figures(missingnessPath, distributionDensityPath, sessionDistributionPath, sessionAuditPath, correlationPath, arenaShiftPath, figDir, cfg, params, repoRoot)
rows = struct([]);
rows = local_add_figure_rows(rows, local_plot_missingness_by_family(missingnessPath, figDir, cfg), ...
    missingnessPath, 'Feature-level missingness for canonical chunking/clustering candidate features.');
rows = local_add_figure_rows(rows, local_plot_distribution_examples(distributionDensityPath, figDir, cfg, params), ...
    distributionDensityPath, 'All canonical feature distributions stratified by arena with binned densities and quantile markers.');
rows = local_add_figure_rows(rows, local_plot_binary_prevalence_by_arena(sessionDistributionPath, figDir, cfg), ...
    sessionDistributionPath, 'Session-level prevalence of boolean clustering candidate features by arena.');
rows = local_add_figure_rows(rows, local_plot_frame_validity(sessionAuditPath, figDir, cfg), ...
    sessionAuditPath, 'Session-level frame validity and bad-frame propagation summary.');
rows = local_add_figure_rows(rows, local_plot_correlation_heatmaps(correlationPath, figDir, cfg), ...
    correlationPath, 'Row-aligned feature redundancy diagnostics using reviewer-safe transformed channels.');
rows = local_add_figure_rows(rows, local_plot_arena_shift(arenaShiftPath, figDir, cfg), ...
    arenaShiftPath, 'Arena-separable feature shift diagnostics without condition labels.');
figureManifest = struct2table(rows);
for i = 1:height(figureManifest)
    figureManifest.figure_file(i) = local_repo_relative_path(repoRoot, figureManifest.figure_file(i));
    figureManifest.source_csv(i) = local_repo_relative_path(repoRoot, figureManifest.source_csv(i));
end
end

function rows = local_add_figure_rows(rows, files, sourceCsv, description)
for i = 1:numel(files)
    row = struct();
    row.figure_file = string(files(i));
    row.source_csv = string(sourceCsv);
    row.description = string(description);
    rows = [rows; row]; %#ok<AGROW>
end
end

function files = local_plot_missingness_by_family(missingnessPath, figDir, cfg)
T = readtable(missingnessPath, 'TextType', 'string');
if any(strcmp(T.Properties.VariableNames, 'ClusteringCandidate'))
    T = T(T.ClusteringCandidate == 1, :);
end
T.feature_type = repmat("continuous", height(T), 1);
T.feature_type(T.IsCircular == 1) = "circular";
T.feature_type(T.IsBoolean == 1) = "boolean";
T.feature_type(T.Unit == "signed_index") = "signed_index";

commonValid = mean(T.valid_frame_fraction_overall, 'omitnan') .* 100;
commonInvalid = 100 - commonValid;
T = sortrows(T, {'missing_fraction_valid_frames','feature_type','FeatureIndex'}, ...
    {'ascend','ascend','ascend'});

fig = figure('Visible','off', 'Color','w', 'Position',[80 80 1800 880]);
tl = tiledlayout(fig, 1, 3, 'TileSpacing','compact', 'Padding','compact');
title(tl, 'Frame-mask propagation and residual missingness across clustering-candidate features', ...
    'FontName', char(cfg.figures.font_name), ...
    'FontSize', cfg.figures.title_font_size, 'FontWeight','bold', ...
    'Interpreter','none');

nexttile(tl, 1);
hold on;
rectangle('Position', [0 0.75 commonValid 0.5], ...
    'FaceColor', [0.00 0.62 0.45], 'EdgeColor','none');
rectangle('Position', [commonValid 0.75 commonInvalid 0.5], ...
    'FaceColor', [0.25 0.25 0.25], 'EdgeColor','none');
xlim([0 100]);
ylim([0.4 1.6]);
set(gca, 'YTick', 1, 'YTickLabel', "all feature frames", ...
    'TickLabelInterpreter','none');
xlabel('% of frame-feature positions', 'Interpreter','none');
title('Common frame mask before feature-specific checks', 'Interpreter','none');
text(commonValid/2, 1, sprintf('valid %.1f%%', commonValid), ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
    'Color','w', 'FontName', char(cfg.figures.font_name), ...
    'FontSize', cfg.figures.font_size, 'Interpreter','none');
text(commonValid + commonInvalid/2, 1, sprintf('masked %.1f%%', commonInvalid), ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
    'Color','w', 'FontName', char(cfg.figures.font_name), ...
    'FontSize', cfg.figures.font_size, 'Interpreter','none');
hFrameLegend = gobjects(2, 1);
hFrameLegend(1) = plot(nan, nan, 's', 'MarkerFaceColor', [0.00 0.62 0.45], ...
    'MarkerEdgeColor','none', 'MarkerSize', 8);
hFrameLegend(2) = plot(nan, nan, 's', 'MarkerFaceColor', [0.25 0.25 0.25], ...
    'MarkerEdgeColor','none', 'MarkerSize', 8);
legend(hFrameLegend, {'valid frames','preprocessing/core-coordinate mask'}, ...
    'Location','southoutside', 'Orientation','horizontal', ...
    'Box','off', 'Interpreter','none');
local_style_axes(gca, cfg);

nexttile(tl, 2, [1 2]);
hold on;
y = T.missing_fraction_valid_frames .* 100;
x = 1:height(T);
colors = local_feature_type_colors(T.feature_type);
xMax = max([0.35; y .* 1.25]);
b = barh(x, y, 'FaceColor','flat', 'EdgeColor','none');
b.CData = colors;
scatter(repmat(xMax * 0.006, height(T), 1), x(:), 18, colors, ...
    'filled', 'MarkerEdgeColor','w', 'LineWidth', 0.25);
set(gca, 'YTick', x, 'YTickLabel', local_display_feature_label(T.Name), ...
    'TickLabelInterpreter','none');
xlim([0 xMax]);
ylim([0.25 height(T)+0.75]);
xlabel('Residual nonfinite values (% of valid frames)', 'Interpreter','none');
title(sprintf('Feature-specific missingness after frame masking (n = %d features)', height(T)), ...
    'Interpreter','none');
for i = 1:height(T)
    if y(i) > 0
        text(y(i) + xMax * 0.01, x(i), sprintf('%.3f', y(i)), ...
            'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
            'FontName', char(cfg.figures.font_name), ...
            'FontSize', max(6, cfg.figures.font_size - 2), ...
            'Interpreter','none');
    end
end
legendTypes = ["continuous","circular","boolean","signed_index"];
legendHandles = gobjects(numel(legendTypes), 1);
for iType = 1:numel(legendTypes)
    cType = local_feature_type_colors(legendTypes(iType));
    legendHandles(iType) = plot(nan, nan, 's', 'MarkerFaceColor', cType, ...
        'MarkerEdgeColor','none', 'MarkerSize', 8);
end
legendLabels = cellstr(replace(legendTypes, "_", " "));
legend(legendHandles, legendLabels, 'Location','southoutside', ...
    'Orientation','horizontal', 'Box','off', 'Interpreter','none');
local_style_axes(gca, cfg);
files = local_export_figure(fig, figDir, 'feature_missingness_by_family', cfg);
close(fig);
end

function files = local_plot_distribution_examples(distributionPath, figDir, cfg, params)
D = readtable(distributionPath, 'TextType', 'string');
if isstring(params.distribution_example_features) && ...
        isscalar(params.distribution_example_features) && ...
        lower(params.distribution_example_features) == "all"
    [~, firstIdx] = unique(D.feature_index, 'stable');
    featureOrder = sortrows(D(firstIdx, {'feature_index','feature_name'}), 'feature_index');
    exampleFeatures = string(featureOrder.feature_name);
else
    exampleFeatures = params.distribution_example_features(:);
    exampleFeatures = exampleFeatures(ismember(exampleFeatures, string(D.feature_name)));
end
if isempty(exampleFeatures)
    files = strings(0,1);
    return
end

nFeat = numel(exampleFeatures);
if nFeat > 18
    nCol = 3;
else
    nCol = min(4, nFeat);
end
nRow = ceil(nFeat / nCol);
fig = figure('Visible','off', 'Color','w', 'Position',[80 80 2100 315*nRow]);
tl = tiledlayout(fig, nRow, nCol, 'TileSpacing','compact', 'Padding','compact');
title(tl, 'All canonical dyadic feature distributions by arena', ...
    'FontName', char(cfg.figures.font_name), ...
    'FontSize', cfg.figures.title_font_size, 'FontWeight','bold', ...
    'Interpreter','none');
arenaLevels = unique(D.arena_label, 'stable');
arenaColors = [0.00 0.45 0.70; 0.80 0.47 0.10; 0.35 0.35 0.35];
legendHandles = gobjects(numel(arenaLevels), 1);

for i = 1:nFeat
    nexttile(tl);
    hold on;
    if i == 1
        for a = 1:numel(arenaLevels)
            c = arenaColors(1 + mod(a-1, size(arenaColors,1)), :);
            legendHandles(a) = plot(nan, nan, '-', 'Color', c, 'LineWidth', 2);
        end
    end
    F = D(D.feature_name == exampleFeatures(i), :);
    if isempty(F)
        continue
    end
    isBoolean = logical(F.is_boolean(1));
    maxY = 0;
    for a = 1:numel(arenaLevels)
        A = F(F.arena_label == arenaLevels(a), :);
        if isempty(A)
            continue
        end
        c = arenaColors(1 + mod(a-1, size(arenaColors,1)), :);
        A = sortrows(A, 'bin_center');
        if isBoolean
            oneBin = A(A.bin_center == 1, :);
            prevalence = 0;
            if ~isempty(oneBin)
                prevalence = oneBin.probability(1);
            end
            bar(a, prevalence, 0.62, 'FaceColor', c, 'FaceAlpha', 0.82, ...
                'EdgeColor', 'none');
            plot(a, prevalence, 'o', 'Color', c .* 0.65, ...
                'MarkerFaceColor', c, 'MarkerSize', 4);
            maxY = max(maxY, prevalence);
        else
            x = A.bin_center;
            y = A.density;
            fill([x; flipud(x)], [zeros(size(y)); flipud(y)], c, ...
                'FaceAlpha', 0.18, 'EdgeColor','none');
            plot(x, y, '-', 'Color', c, 'LineWidth', 1.7);
            maxY = max(maxY, max(y));
            q25 = A.q25(1);
            q75 = A.q75(1);
            med = A.median(1);
            if isfinite(q25) && isfinite(q75)
                plot([q25 q75], [0 0] + max(y)*0.08, '-', 'Color', c, ...
                    'LineWidth', 5);
            end
            if isfinite(med)
                plot([med med], [0 max(y)*0.20], '-', 'Color', c, ...
                    'LineWidth', 1.5);
            end
        end
    end

    if isBoolean
        set(gca, 'XTick', 1:numel(arenaLevels), 'XTickLabel', arenaLevels, ...
            'TickLabelInterpreter','none');
        xlim([0.4 numel(arenaLevels)+0.6]);
        ylim([0 max(0.01, min(1, maxY*1.30))]);
        ylabel('P(feature = 1)', 'Interpreter','none');
        xlabel('arena', 'Interpreter','none');
    else
        xlim([F.display_min(1), F.display_max(1)]);
        ylim([0 max(eps, maxY*1.16)]);
        ylabel('density', 'Interpreter','none');
        xlabel(local_display_unit(F.unit(1)), 'Interpreter','none');
    end
    title(local_display_feature_label(exampleFeatures(i)), 'Interpreter','none');
    local_style_axes(gca, cfg);
end
lgd = legend(legendHandles, cellstr(arenaLevels), 'Location','southoutside', ...
    'Orientation','horizontal', 'Box','off');
lgd.Layout.Tile = 'south';
files = local_export_figure(fig, figDir, 'feature_distribution_examples_by_arena', cfg);
close(fig);
end

function files = local_plot_binary_prevalence_by_arena(sessionDistributionPath, figDir, cfg)
T = readtable(sessionDistributionPath, 'TextType', 'string');
T = T(T.is_boolean == 1, :);
if isempty(T)
    files = strings(0,1);
    return
end

[~, firstIdx] = unique(T.feature_index, 'stable');
featureOrder = sortrows(T(firstIdx, {'feature_index','feature_name'}), 'feature_index');
featureNames = string(featureOrder.feature_name);
arenaLevels = unique(T.arena_label, 'stable');
arenaColors = [0.00 0.45 0.70; 0.80 0.47 0.10; 0.35 0.35 0.35];
offsets = linspace(-0.22, 0.22, max(numel(arenaLevels), 1));

fig = figure('Visible','off', 'Color','w', 'Position',[80 80 1250 620]);
hold on;
legendHandles = gobjects(numel(arenaLevels), 1);
for a = 1:numel(arenaLevels)
    c = arenaColors(1 + mod(a-1, size(arenaColors,1)), :);
    legendHandles(a) = plot(nan, nan, 'o-', 'Color', c, ...
        'MarkerFaceColor', c, 'MarkerEdgeColor','w', 'LineWidth', 2);
end

for f = 1:numel(featureNames)
    for a = 1:numel(arenaLevels)
        S = T(T.feature_name == featureNames(f) & T.arena_label == arenaLevels(a), :);
        if isempty(S)
            continue
        end
        c = arenaColors(1 + mod(a-1, size(arenaColors,1)), :);
        y = S.prevalence_true .* 100;
        xCenter = f + offsets(a);
        jitter = linspace(-0.035, 0.035, max(numel(y), 1))';
        scatter(repmat(xCenter, numel(y), 1) + jitter, y, 42, ...
            'MarkerFaceColor', c, 'MarkerEdgeColor','w', ...
            'MarkerFaceAlpha', 0.82, 'LineWidth', 0.6);
        med = median(y, 'omitnan');
        if isfinite(med)
            plot([xCenter-0.11 xCenter+0.11], [med med], '-', ...
                'Color', c, 'LineWidth', 2.4);
        end
    end
end

set(gca, 'XTick', 1:numel(featureNames), ...
    'XTickLabel', local_display_feature_label(featureNames), ...
    'TickLabelInterpreter','none');
xtickangle(35);
xlim([0.45 numel(featureNames)+0.55]);
ylim([0 100]);
ylabel('Session prevalence (% valid frames)', 'Interpreter','none');
xlabel('Boolean frame-level features used as clustering candidates', 'Interpreter','none');
title('Boolean feature prevalence by arena');
legend(legendHandles, cellstr(arenaLevels), 'Location','northoutside', ...
    'Orientation','horizontal', 'Box','off');
local_style_axes(gca, cfg);
files = local_export_figure(fig, figDir, 'feature_binary_prevalence_by_arena', cfg);
close(fig);
end

function files = local_plot_frame_validity(sessionAuditPath, figDir, cfg)
T = readtable(sessionAuditPath, 'TextType', 'string');
ok = T.status == "success";
T = T(ok,:);
if isempty(T)
    files = strings(0,1);
    return
end
T = sortrows(T, {'arena_label','raw_index'});
valid = T.valid_frame_fraction .* 100;
bad = T.badframe_mask_fraction .* 100;
coreOnly = T.core_nan_only_fraction .* 100;
otherInvalid = max(0, 100 - valid - bad - coreOnly);
Y = [valid, bad, coreOnly, otherInvalid];

fig = figure('Visible','off', 'Color','w', 'Position',[80 80 1400 760]);
b = bar(Y, 'stacked', 'EdgeColor','none');
b(1).FaceColor = [0.00 0.62 0.45];
b(2).FaceColor = [0.20 0.20 0.20];
b(3).FaceColor = [0.82 0.37 0.00];
b(4).FaceColor = [0.65 0.65 0.65];
set(gca, 'XTick', 1:height(T), 'XTickLabel', string(T.raw_index));
xtickangle(45);
ylabel('% frames', 'Interpreter','none');
xlabel('Raw session index', 'Interpreter','none');
ylim([0 100]);
legend({'valid','preprocessing bad frame','core coordinate missing','other invalid'}, ...
    'Location','southoutside', 'Orientation','horizontal', 'Box','off');
title('Frame validity and bad-frame propagation');
local_style_axes(gca, cfg);
files = local_export_figure(fig, figDir, 'feature_frame_validity_summary', cfg);
close(fig);
end

function files = local_plot_correlation_heatmaps(correlationPath, figDir, cfg)
C = readtable(correlationPath, 'TextType', 'string');
if isempty(C)
    files = strings(0,1);
    return
end

files = strings(0,1);
Call = C(C.arena_label == "all", :);
if ~isempty(Call)
    [M, labels] = local_correlation_matrix_from_table(Call);
    fig = figure('Visible','off', 'Color','w', 'Position',[80 80 1300 1100]);
    imagesc(M, [-1 1]);
    axis square;
    colormap(local_correlation_colormap());
    colorbar;
    set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels, ...
        'YTick', 1:numel(labels), 'YTickLabel', labels, 'FontSize', 6, ...
        'TickLabelInterpreter','none');
    xtickangle(90);
    title('Feature redundancy diagnostic channels: all arenas');
    local_style_axes(gca, cfg);
    files = [files; local_export_figure(fig, figDir, 'feature_correlation_heatmap_all_arenas', cfg)];
    close(fig);
end

arenaLevels = unique(C.arena_label(C.arena_label ~= "all"), 'stable');
if ~isempty(arenaLevels)
    fig = figure('Visible','off', 'Color','w', ...
        'Position',[80 80 max(900, 760*numel(arenaLevels)) 900]);
    tl = tiledlayout(fig, 1, numel(arenaLevels), 'TileSpacing','compact', 'Padding','compact');
    for a = 1:numel(arenaLevels)
        nexttile(tl);
        Ca = C(C.arena_label == arenaLevels(a), :);
        [M, labels] = local_correlation_matrix_from_table(Ca);
        imagesc(M, [-1 1]);
        axis square;
        colormap(local_correlation_colormap());
        colorbar;
        set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels, ...
            'YTick', 1:numel(labels), 'YTickLabel', labels, 'FontSize', 5, ...
            'TickLabelInterpreter','none');
        xtickangle(90);
        title(sprintf('%s arena', arenaLevels(a)));
        local_style_axes(gca, cfg);
    end
    files = [files; local_export_figure(fig, figDir, 'feature_correlation_heatmap_by_arena', cfg)];
    close(fig);
end
end

function [M, labels] = local_correlation_matrix_from_table(T)
n = max([T.channel_i_index; T.channel_j_index]);
M = nan(n, n);
labels = strings(n, 1);
for r = 1:height(T)
    i = T.channel_i_index(r);
    j = T.channel_j_index(r);
    M(i, j) = T.spearman_r(r);
    M(j, i) = T.spearman_r(r);
    labels(i) = T.channel_i(r);
    labels(j) = T.channel_j(r);
end
labels = replace(labels, "_", " ");
end

function files = local_plot_arena_shift(arenaShiftPath, figDir, cfg)
T = readtable(arenaShiftPath, 'TextType', 'string');
if isempty(T)
    files = strings(0,1);
    return
end

metrics = ["robust_median_difference_iqr", "prevalence_difference", ...
    "circular_mean_difference_deg"];
metricLabels = ["Continuous features: median shift / pooled IQR", ...
    "Boolean features: prevalence shift", ...
    "Circular features: mean direction shift (degrees)"];
metricAxisLabels = ["median shift / pooled IQR", ...
    "prevalence shift", ...
    "circular mean shift (degrees)"];

fig = figure('Visible','off', 'Color','w', 'Position',[80 80 1500 1050]);
tl = tiledlayout(fig, 1, 3, 'TileSpacing','compact', 'Padding','compact');
for m = 1:numel(metrics)
    nexttile(tl);
    S = T(T.metric_name == metrics(m) & isfinite(T.metric_value), :);
    if isempty(S)
        axis off;
        title(metricLabels(m), 'Interpreter','none');
        continue
    end
    [~, ord] = sort(abs(S.metric_value), 'ascend');
    S = S(ord, :);
    colors = local_family_colors(S.family);
    b = barh(S.metric_value, 'FaceColor','flat', 'EdgeColor','none');
    b.CData = colors;
    set(gca, 'YTick', 1:height(S), 'YTickLabel', replace(S.feature_name, "_", " "), ...
        'TickLabelInterpreter','none');
    hold on;
    yl = ylim;
    plot([0 0], yl, 'k-', 'LineWidth', 0.8);
    ylim(yl);
    xlabel(metricAxisLabels(m), 'Interpreter','none');
    title(metricLabels(m), 'Interpreter','none');
    local_style_axes(gca, cfg);
end
files = local_export_figure(fig, figDir, 'feature_arena_shift_summary', cfg);
close(fig);
end

function label = local_display_feature_label(name)
label = replace(string(name), "__", " ");
label = replace(label, "_", " ");
end

function label = local_display_unit(unitName)
unitName = string(unitName);
switch unitName
    case "mm_per_s"
        label = "mm/s";
    case "signed_index"
        label = "signed index";
    case "unitless"
        label = "unitless";
    otherwise
        label = replace(unitName, "_", " ");
end
end

function featureType = local_feature_type_from_meta(featureMeta, rowIdx)
if featureMeta.IsBoolean(rowIdx) == 1
    featureType = "boolean";
elseif featureMeta.IsCircular(rowIdx) == 1
    featureType = "circular";
elseif string(featureMeta.Unit(rowIdx)) == "signed_index"
    featureType = "signed_index";
else
    featureType = "continuous";
end
end

function colors = local_feature_type_colors(featureTypes)
featureTypes = string(featureTypes);
colors = zeros(numel(featureTypes), 3);
for i = 1:numel(featureTypes)
    switch featureTypes(i)
        case "boolean"
            colors(i,:) = [0.80 0.47 0.10];
        case "circular"
            colors(i,:) = [0.35 0.55 0.35];
        case "signed_index"
            colors(i,:) = [0.45 0.30 0.65];
        otherwise
            colors(i,:) = [0.00 0.45 0.70];
    end
end
end

function colors = local_family_colors(families)
[familyLevels, ~, groupIdx] = unique(families, 'stable');
palette = lines(max(numel(familyLevels), 1));
colors = palette(groupIdx, :);
end

function cmap = local_correlation_colormap()
n = 256;
neg = [linspace(0.10, 0.95, n/2)' linspace(0.32, 0.95, n/2)' linspace(0.65, 0.95, n/2)'];
pos = [linspace(0.95, 0.70, n/2)' linspace(0.95, 0.12, n/2)' linspace(0.95, 0.10, n/2)'];
cmap = [neg; pos];
end

function files = local_export_figure(fig, figDir, baseName, cfg)
files = strings(0,1);
if cfg.figures.export_png
    pngFile = fullfile(figDir, [baseName '.png']);
    exportgraphics(fig, pngFile, 'Resolution', cfg.figures.dpi);
    files(end+1,1) = string(pngFile);
end
if cfg.figures.export_pdf
    pdfFile = fullfile(figDir, [baseName '.pdf']);
    exportgraphics(fig, pdfFile, 'ContentType','vector');
    files(end+1,1) = string(pdfFile);
end
end

function local_style_axes(ax, cfg)
set(ax, 'FontName', char(cfg.figures.font_name), 'FontSize', cfg.figures.font_size, ...
    'LineWidth', 0.8, 'Box','off', 'TickDir','out');
grid(ax, 'on');
ax.GridColor = local_hex_to_rgb(cfg.figures.grid_color);
ax.GridAlpha = 0.35;
title(ax, ax.Title.String, 'FontSize', cfg.figures.title_font_size, ...
    'FontWeight','bold', 'Interpreter','none');
end

function rgb = local_hex_to_rgb(hex)
hex = char(strtrim(string(hex)));
rgb = sscanf(hex(2:end), '%2x%2x%2x', [1 3]) ./ 255;
end
