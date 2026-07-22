function outputs = run_condition_blind_embedding_build(repoRoot, opts)
%RUN_CONDITION_BLIND_EMBEDDING_BUILD Build run-07 embedding artifacts.
%
% Condition-blind contract
%   - Primary scale and anchor inputs come from run_06 primary scale-specific
%     outputs.
%   - Condition, cohort, group, drug, genotype, outcome, and arena labels do
%     not enter matrix construction, PCA, global PCA, scale weights, or
%     stability resampling.
%   - Arena labels are joined only after fitting for audit diagnostics.

if nargin < 1 || strlength(string(repoRoot)) == 0
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
end
if nargin < 2 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'configPath') || strlength(string(opts.configPath)) == 0
    opts.configPath = fullfile(repoRoot, 'config', 'multiscale_embedding_config.csv');
end

repoRoot = string(repoRoot);
cd(repoRoot);
addpath(genpath(repoRoot));

params = load_multiscale_embedding_config(opts.configPath);
outRoot = resolve_repo_path(repoRoot, params.output_dir);
figDir = fullfile(outRoot, 'figures');
logDir = fullfile(outRoot, 'logs');
local_ensure_dir(outRoot);
local_ensure_dir(figDir);
local_ensure_dir(logDir);

paths = local_output_paths(outRoot, params);
diary(fullfile(logDir, 'run_07_build_condition_blind_embedding_latest.log'));
cleanup = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('run_07_build_condition_blind_embedding\n');
fprintf('Repo root: %s\n', repoRoot);
fprintf('Output root: %s\n', outRoot);
fprintf('Run mode: %s\n', params.run_mode);
fprintf('Anchor manifest mode: %s\n', params.anchor_manifest_mode);

[inputRoot, usedFallback] = local_resolve_input_root(repoRoot, params);
fprintf('Run_06 input root: %s\n', inputRoot);
fprintf('Reviewed snapshot fallback used: %d\n', usedFallback);

writetable(params.config_table, paths.parameterAudit);

[primaryScales, primaryAnchorsFull, primaryBankSummary, eventSummary, ...
    dimensionAudit, pcaLoadingStability, arenaSensitivityRun06, transformAudit, ...
    featureDict, featureSessionAudit] = local_load_inputs(repoRoot, inputRoot, params);
local_assert_label_free_inputs(primaryScales, primaryAnchorsFull, eventSummary, transformAudit);

[primaryAnchors, anchorAudit] = local_select_run07_anchors(primaryAnchorsFull, params);
inputScaleAudit = local_input_scale_audit(primaryScales, dimensionAudit, pcaLoadingStability, ...
    arenaSensitivityRun06, inputRoot, usedFallback);
inputAnchorAudit = local_input_anchor_audit(primaryAnchorsFull, primaryAnchors, primaryBankSummary, inputRoot, params);
anchorWeightAudit = local_anchor_weight_audit(primaryAnchors, params);
rareCoverageAudit = local_rare_coverage_audit(primaryAnchors, params);
featureDictionaryAudit = local_feature_dictionary_header(featureDict, featureSessionAudit);

writetable(inputScaleAudit, paths.inputScaleAudit);
writetable(inputAnchorAudit, paths.inputAnchorAudit);
writetable(anchorWeightAudit, paths.anchorWeightAudit);
writetable(rareCoverageAudit, paths.rareCoverageAudit);

fprintf('Primary scales selected for embedding: %d\n', height(primaryScales));
fprintf('Primary anchor rows available: %d\n', height(primaryAnchorsFull));
fprintf('Run_07 anchor rows used: %d\n', height(primaryAnchors));

Input = build_primary_scale_embedding_inputs(repoRoot, primaryScales, primaryAnchors, ...
    eventSummary, transformAudit, featureDict, dimensionAudit, params);
[Embedding, Audit] = fit_condition_blind_embedding_pca(Input, dimensionAudit, params);

featureDictionary = Input.featureDictionary;
featureDictionary.run07_feature_dictionary_role = repmat("embedding_matrix_column_dictionary", height(featureDictionary), 1);
featureDictionary.labels_used_for_feature_construction = repmat("none", height(featureDictionary), 1);
featureDictionary.arena_used_for_feature_construction = false(height(featureDictionary), 1);
featureDictionary.condition_used_for_feature_construction = false(height(featureDictionary), 1);
featureDictionaryAudit = [featureDictionaryAudit; local_feature_dictionary_summary_rows(featureDictionary)]; %#ok<AGROW>

writetable(featureDictionary, paths.featureDictionary);
writetable(Audit.pcaByScale, paths.pcaByScale);
writetable(Audit.scaleWeights, paths.scaleWeights);
writetable(Audit.stability, paths.stabilityAudit);
writetable(Audit.arenaSensitivity, paths.arenaSensitivity);
writetable(Audit.preprocessDimensions, paths.preprocessDimensionAudit);
writetable(Audit.anchorStageSensitivity, paths.anchorStageSensitivityAudit);
writetable(Embedding.rowManifest, paths.rowManifest);
writetable(Embedding.scoreTable, paths.globalScores);

save(paths.summaryMatrixMat, 'Input', 'params', 'inputScaleAudit', 'inputAnchorAudit', '-v7.3');
save(paths.embeddingModelMat, 'Embedding', 'Audit', 'params', 'inputScaleAudit', ...
    'inputAnchorAudit', 'featureDictionary', 'featureDictionaryAudit', '-v7.3');

matrixManifest = local_matrix_manifest(paths, Input, Embedding, Audit, inputRoot, params);
writetable(matrixManifest, paths.matrixManifest);

figureManifest = make_run07_embedding_qc_figures(outRoot, params);
writetable(figureManifest, paths.figureManifest);

fprintf('Wrote parameter audit: %s\n', paths.parameterAudit);
fprintf('Wrote input scale audit: %s\n', paths.inputScaleAudit);
fprintf('Wrote input anchor audit: %s\n', paths.inputAnchorAudit);
fprintf('Wrote anchor weight audit: %s\n', paths.anchorWeightAudit);
fprintf('Wrote rare-strata coverage audit: %s\n', paths.rareCoverageAudit);
fprintf('Wrote feature dictionary: %s\n', paths.featureDictionary);
fprintf('Wrote matrix manifest: %s\n', paths.matrixManifest);
fprintf('Wrote PCA by scale: %s\n', paths.pcaByScale);
fprintf('Wrote scale weights: %s\n', paths.scaleWeights);
fprintf('Wrote stability audit: %s\n', paths.stabilityAudit);
fprintf('Wrote arena sensitivity audit: %s\n', paths.arenaSensitivity);
fprintf('Wrote preprocessing dimension audit: %s\n', paths.preprocessDimensionAudit);
fprintf('Wrote anchor-stage PCA sensitivity audit: %s\n', paths.anchorStageSensitivityAudit);
fprintf('Wrote figure manifest: %s\n', paths.figureManifest);
fprintf('Wrote ignored MAT artifacts: %s | %s\n', paths.summaryMatrixMat, paths.embeddingModelMat);

outputs = struct();
outputs.output_root = string(outRoot);
outputs.input_root = string(inputRoot);
outputs.used_reviewed_snapshot_fallback = usedFallback;
outputs.parameter_audit_path = string(paths.parameterAudit);
outputs.input_scale_audit_path = string(paths.inputScaleAudit);
outputs.input_anchor_audit_path = string(paths.inputAnchorAudit);
outputs.anchor_weight_audit_path = string(paths.anchorWeightAudit);
outputs.rare_coverage_audit_path = string(paths.rareCoverageAudit);
outputs.feature_dictionary_path = string(paths.featureDictionary);
outputs.matrix_manifest_path = string(paths.matrixManifest);
outputs.pca_by_scale_path = string(paths.pcaByScale);
outputs.scale_weight_audit_path = string(paths.scaleWeights);
outputs.stability_audit_path = string(paths.stabilityAudit);
outputs.arena_sensitivity_audit_path = string(paths.arenaSensitivity);
outputs.preprocess_dimension_audit_path = string(paths.preprocessDimensionAudit);
outputs.anchor_stage_sensitivity_audit_path = string(paths.anchorStageSensitivityAudit);
outputs.figure_manifest_path = string(paths.figureManifest);
outputs.global_scores_path = string(paths.globalScores);
outputs.embedding_model_mat_path = string(paths.embeddingModelMat);
outputs.summary_matrix_mat_path = string(paths.summaryMatrixMat);
outputs.n_primary_scales = height(primaryScales);
outputs.n_input_anchors = height(primaryAnchors);
outputs.n_global_rows = size(Embedding.global_matrix, 1);
outputs.n_global_dims = size(Embedding.global_matrix, 2);
outputs.n_global_pcs = size(Embedding.global_score, 2);
end

function paths = local_output_paths(outRoot, params)
paths = struct();
paths.parameterAudit = fullfile(outRoot, 'embedding_parameter_audit.csv');
paths.inputScaleAudit = fullfile(outRoot, 'embedding_input_scale_audit.csv');
paths.inputAnchorAudit = fullfile(outRoot, 'embedding_input_anchor_audit.csv');
paths.anchorWeightAudit = fullfile(outRoot, 'embedding_anchor_weight_audit.csv');
paths.rareCoverageAudit = fullfile(outRoot, 'embedding_rare_strata_coverage_audit.csv');
paths.featureDictionary = fullfile(outRoot, 'embedding_feature_dictionary.csv');
paths.matrixManifest = fullfile(outRoot, 'embedding_matrix_manifest.csv');
paths.pcaByScale = fullfile(outRoot, 'embedding_pca_by_scale.csv');
paths.scaleWeights = fullfile(outRoot, 'embedding_scale_weight_audit.csv');
paths.stabilityAudit = fullfile(outRoot, 'embedding_stability_audit.csv');
paths.arenaSensitivity = fullfile(outRoot, 'embedding_arena_sensitivity_audit.csv');
paths.preprocessDimensionAudit = fullfile(outRoot, 'embedding_preprocess_dimension_audit.csv');
paths.anchorStageSensitivityAudit = fullfile(outRoot, 'embedding_anchor_stage_pca_sensitivity_audit.csv');
paths.figureManifest = fullfile(outRoot, 'embedding_qc_figure_manifest.csv');
paths.rowManifest = fullfile(outRoot, 'embedding_row_manifest.csv');
paths.globalScores = fullfile(outRoot, 'embedding_global_scores.csv');
paths.embeddingModelMat = fullfile(outRoot, char(params.embedding_model_mat_file));
paths.summaryMatrixMat = fullfile(outRoot, char(params.summary_matrix_mat_file));
end

function [inputRoot, usedFallback] = local_resolve_input_root(repoRoot, params)
inputRoot = resolve_repo_path(repoRoot, params.chunk_input_dir);
usedFallback = false;
if local_has_required_run06_files(inputRoot, params)
    return
end
fallbackRoot = resolve_repo_path(repoRoot, params.fallback_chunk_input_dir);
if logical(params.allow_reviewed_snapshot_fallback) && local_has_required_run06_files(fallbackRoot, params)
    inputRoot = fallbackRoot;
    usedFallback = true;
    return
end
error('run_condition_blind_embedding_build:MissingRun06PrimaryOutputs', ...
    ['Preferred run_06 input root lacks required primary outputs: %s\n' ...
    'Fallback root checked: %s'], inputRoot, fallbackRoot);
end

function tf = local_has_required_run06_files(rootDir, params)
files = local_required_run06_files(params);
tf = true;
for i = 1:numel(files)
    tf = tf && isfile(fullfile(rootDir, char(files(i))));
end
end

function files = local_required_run06_files(params)
F = local_active_run06_files(params);
files = [params.primary_scales_file, F.anchor_manifest, ...
    F.bank_summary, F.event_summary, ...
    params.embedding_dimension_audit_file, params.pca_loading_stability_file, ...
    params.arena_sensitivity_audit_file, params.chunk_feature_transform_audit_file];
if params.anchor_manifest_mode == "rare_enriched"
    files = [files, params.expanded_weight_audit_file, params.rare_strata_sampling_summary_file];
end
end

function F = local_active_run06_files(params)
F = struct();
if params.anchor_manifest_mode == "rare_enriched"
    F.anchor_manifest = string(params.expanded_anchor_manifest_file);
    F.bank_summary = string(params.expanded_chunk_bank_summary_file);
    F.event_summary = string(params.expanded_event_summary_file);
else
    F.anchor_manifest = string(params.primary_anchor_manifest_file);
    F.bank_summary = string(params.primary_chunk_bank_summary_file);
    F.event_summary = string(params.primary_event_summary_file);
end
end

function [primaryScales, primaryAnchors, primaryBankSummary, eventSummary, ...
    dimensionAudit, pcaLoadingStability, arenaSensitivityRun06, transformAudit, ...
    featureDict, featureSessionAudit] = local_load_inputs(repoRoot, inputRoot, params)
F = local_active_run06_files(params);
primaryScales = local_read_csv(fullfile(inputRoot, char(params.primary_scales_file)));
primaryAnchors = local_read_csv(fullfile(inputRoot, char(F.anchor_manifest)));
primaryBankSummary = local_read_csv(fullfile(inputRoot, char(F.bank_summary)));
eventSummary = local_read_csv(fullfile(inputRoot, char(F.event_summary)));
dimensionAudit = local_read_csv(fullfile(inputRoot, char(params.embedding_dimension_audit_file)));
pcaLoadingStability = local_read_csv(fullfile(inputRoot, char(params.pca_loading_stability_file)));
arenaSensitivityRun06 = local_read_csv(fullfile(inputRoot, char(params.arena_sensitivity_audit_file)));
transformAudit = local_read_csv(fullfile(inputRoot, char(params.chunk_feature_transform_audit_file)));

featureRoot = resolve_repo_path(repoRoot, params.feature_output_dir);
featureDict = local_read_csv(fullfile(featureRoot, char(params.feature_dictionary_file)));
featureSessionAudit = local_read_csv(fullfile(featureRoot, char(params.feature_session_audit_file)));
end

function local_assert_label_free_inputs(primaryScales, primaryAnchors, eventSummary, transformAudit)
local_assert_string_none(primaryScales, 'labels_used_for_primary_scale_rule', 'primary scale rule');
local_assert_false(primaryScales, 'arena_used_for_primary_scale_rule', 'primary scale rule arena');
local_assert_false(primaryScales, 'condition_used_for_primary_scale_rule', 'primary scale rule condition');
local_assert_string_none(primaryAnchors, 'labels_used_for_primary_anchor_selection', 'primary anchor selection');
local_assert_false(primaryAnchors, 'arena_used_for_primary_anchor_selection', 'primary anchor selection arena');
local_assert_false(primaryAnchors, 'condition_used_for_primary_anchor_selection', 'primary anchor selection condition');
local_assert_string_none(primaryAnchors, 'labels_used_for_anchor_selection', 'expanded anchor selection');
local_assert_false(primaryAnchors, 'arena_used_for_anchor_selection', 'expanded anchor selection arena');
local_assert_false(primaryAnchors, 'condition_used_for_anchor_selection', 'expanded anchor selection condition');
local_assert_string_none(eventSummary, 'labels_used_for_event_summary', 'event summary');
local_assert_false(eventSummary, 'arena_used_for_event_summary', 'event summary arena');
local_assert_false(eventSummary, 'condition_used_for_event_summary', 'event summary condition');
if ismember('condition_blind_transform', transformAudit.Properties.VariableNames)
    assert(all(logical(transformAudit.condition_blind_transform)), ...
        'run_condition_blind_embedding_build:TransformNotConditionBlind', ...
        'Run_06 transform audit is not marked condition-blind.');
end
if ismember('final_inclusion_probability', primaryAnchors.Properties.VariableNames)
    requiredWeights = ["base_inclusion_probability","rare_inclusion_probability", ...
        "final_inclusion_probability","pca_training_weight","graph_training_weight", ...
        "audit_inverse_probability_weight"];
    missingWeights = setdiff(requiredWeights, string(primaryAnchors.Properties.VariableNames));
    assert(isempty(missingWeights), 'run_condition_blind_embedding_build:MissingAnchorWeights', ...
        'Expanded anchors are missing required weight columns: %s', strjoin(missingWeights, ', '));
    for name = requiredWeights
        assert(all(isfinite(double(primaryAnchors.(name)))), ...
            'run_condition_blind_embedding_build:NonfiniteAnchorWeights', ...
            'Expanded anchor weight %s must be finite.', name);
    end
end
end

function local_assert_string_none(T, varName, context)
if ismember(varName, T.Properties.VariableNames)
    vals = string(T.(varName));
    assert(all(vals == "none"), 'run_condition_blind_embedding_build:LabelLeakInput', ...
        'Expected %s to report labels_used=none.', context);
end
end

function local_assert_false(T, varName, context)
if ismember(varName, T.Properties.VariableNames)
    assert(~any(logical(T.(varName))), ...
        'run_condition_blind_embedding_build:LabelLeakInput', ...
        'Expected %s to report false.', context);
end
end

function [selected, audit] = local_select_run07_anchors(anchorFull, params)
selected = table();
audit = table();
scales = unique(anchorFull.scale_index, 'stable')';
for s = scales
    A = anchorFull(anchorFull.scale_index == s, :);
    A.source_primary_anchor_manifest_row = find(anchorFull.scale_index == s);
    if params.run_mode == "smoke"
        cap = min(height(A), floor(params.smoke_max_anchors_per_scale));
    elseif isfinite(params.max_anchors_per_scale)
        cap = min(height(A), floor(params.max_anchors_per_scale));
    else
        cap = height(A);
    end
    if cap < height(A)
        pickInput = A(:, {'session_index', 'raw_index', 'anchor_time_s', 'source_primary_anchor_manifest_row'});
        picked = select_condition_blind_anchor_subset(pickInput, cap, ...
            'minAnchorsPerSession', params.min_smoke_anchors_per_session);
        [~, loc] = ismember(picked.source_primary_anchor_manifest_row, A.source_primary_anchor_manifest_row);
        S = A(loc, :);
        rule = "run07_deterministic_time_even_subset_of_selected_anchor_manifest";
    else
        S = A;
        rule = "run07_all_rows_of_selected_anchor_manifest";
    end
    S = sortrows(S, {'raw_index', 'anchor_time_s'});
    S.run07_anchor_selection_rule = repmat(rule, height(S), 1);
    S.run07_labels_used_for_anchor_subset = repmat("none", height(S), 1);
    S.run07_arena_used_for_anchor_subset = false(height(S), 1);
    S.run07_condition_used_for_anchor_subset = false(height(S), 1);
    selected = [selected; S]; %#ok<AGROW>
end
selected = sortrows(selected, {'chunk_sec', 'raw_index', 'anchor_time_s'});
selected.run07_embedding_anchor_id = (1:height(selected))';
end

function T = local_input_scale_audit(primaryScales, dimensionAudit, pcaLoadingStability, arenaSensitivityRun06, inputRoot, usedFallback)
T = sortrows(primaryScales, 'chunk_sec');
T.run06_input_root = repmat(string(inputRoot), height(T), 1);
T.used_reviewed_snapshot_fallback = repmat(logical(usedFallback), height(T), 1);
T.labels_used_for_run07_input_scale_selection = repmat("none", height(T), 1);
T.arena_used_for_run07_input_scale_selection = false(height(T), 1);
T.condition_used_for_run07_input_scale_selection = false(height(T), 1);
T.run07_input_scale_role = repmat("run06_primary_operational_scale", height(T), 1);

[tf, loc] = ismember(T.scale_index, dimensionAudit.scale_index);
T.run06_recommended_pcs_for_next_embedding = nan(height(T), 1);
T.run06_n_pcs_90pct = nan(height(T), 1);
T.run06_summary_profile = strings(height(T), 1);
T.run06_n_temporal_bins = nan(height(T), 1);
T.run06_n_dct_coeffs = nan(height(T), 1);
T.run06_summary_dims = nan(height(T), 1);
T.run06_recommended_pcs_for_next_embedding(tf) = dimensionAudit.recommended_pcs_for_next_embedding(loc(tf));
T.run06_n_pcs_90pct(tf) = dimensionAudit.n_pcs_90pct(loc(tf));
T.run06_summary_profile(tf) = string(dimensionAudit.summary_profile(loc(tf)));
T.run06_n_temporal_bins(tf) = dimensionAudit.n_temporal_bins(loc(tf));
T.run06_n_dct_coeffs(tf) = dimensionAudit.n_dct_coeffs(loc(tf));
T.run06_summary_dims(tf) = dimensionAudit.n_summary_dims(loc(tf));

[tf, loc] = ismember(T.scale_index, pcaLoadingStability.scale_index);
T.run06_pca_loading_median_subspace_similarity = nan(height(T), 1);
T.run06_pca_loading_median_subspace_similarity(tf) = pcaLoadingStability.median_subspace_similarity(loc(tf));

[tf, loc] = ismember(T.scale_index, arenaSensitivityRun06.scale_index);
T.run06_arena_shift_audit_iqr_units = nan(height(T), 1);
T.run06_arena_shift_audit_iqr_units(tf) = arenaSensitivityRun06.embedding_arena_shift_iqr_units(loc(tf));
end

function audit = local_input_anchor_audit(anchorFull, selected, primaryBankSummary, inputRoot, params)
audit = table();
scales = unique(anchorFull.scale_index, 'stable')';
for s = scales
    A = anchorFull(anchorFull.scale_index == s, :);
    S = selected(selected.scale_index == s, :);
    one = table();
    one.scale_index = double(s);
    one.chunk_sec = double(A.chunk_sec(1));
    one.n_primary_anchor_rows_available = height(A);
    one.n_run07_anchor_rows_used = height(S);
    one.n_sessions_available = numel(unique(A.session_index));
    one.n_sessions_used = numel(unique(S.session_index));
    one.mean_valid_frac_used = mean(S.valid_frac, 'omitnan');
    one.min_valid_frac_used = min(S.valid_frac, [], 'omitnan');
    one.mean_feature_finite_frac_used = mean(S.feature_finite_frac, 'omitnan');
    one.min_feature_finite_frac_used = min(S.feature_finite_frac, [], 'omitnan');
    one.run_mode = string(params.run_mode);
    one.anchor_manifest_mode = string(params.anchor_manifest_mode);
    one.run06_input_root = string(inputRoot);
    one.run07_anchor_subset_rule = string(S.run07_anchor_selection_rule(1));
    F = local_active_run06_files(params);
    one.source_anchor_manifest = string(fullfile(inputRoot, char(F.anchor_manifest)));
    one.labels_used_for_anchor_subset = "none";
    one.arena_used_for_anchor_subset = false;
    one.condition_used_for_anchor_subset = false;
    audit = [audit; one]; %#ok<AGROW>
end

if ~isempty(primaryBankSummary) && all(ismember({'n_scale_specific_candidate_anchors', ...
        'n_sessions_with_primary_anchors'}, primaryBankSummary.Properties.VariableNames))
    [tf, loc] = ismember(audit.scale_index, primaryBankSummary.scale_index);
    audit.n_scale_specific_candidate_anchors_run06 = nan(height(audit), 1);
    audit.n_sessions_with_primary_anchors_run06 = nan(height(audit), 1);
    audit.n_scale_specific_candidate_anchors_run06(tf) = primaryBankSummary.n_scale_specific_candidate_anchors(loc(tf));
    audit.n_sessions_with_primary_anchors_run06(tf) = primaryBankSummary.n_sessions_with_primary_anchors(loc(tf));
end
audit = sortrows(audit, 'chunk_sec');
end

function Audit = local_anchor_weight_audit(selected, params)
Audit = table();
scales = unique(double(selected.scale_index), 'stable')';
for s = scales
    idx = double(selected.scale_index) == s;
    if ismember('anchor_stage', selected.Properties.VariableNames)
        stage = string(selected.anchor_stage(idx));
    else
        stage = repmat("base_time_even", nnz(idx), 1);
    end
    if ismember('final_inclusion_probability', selected.Properties.VariableNames)
        baseP = double(selected.base_inclusion_probability(idx));
        rareP = double(selected.rare_inclusion_probability(idx));
        finalP = double(selected.final_inclusion_probability(idx));
        auditW = double(selected.audit_inverse_probability_weight(idx));
        pcaW = double(selected.pca_training_weight(idx));
        graphW = double(selected.graph_training_weight(idx));
    else
        baseP = ones(nnz(idx), 1);
        rareP = zeros(nnz(idx), 1);
        finalP = ones(nnz(idx), 1);
        auditW = ones(nnz(idx), 1);
        pcaW = ones(nnz(idx), 1);
        graphW = ones(nnz(idx), 1);
    end
    assert(all(isfinite([baseP; rareP; finalP; auditW; pcaW; graphW])) && ...
        all(finalP > 0 & finalP <= 1) && all(auditW > 0), ...
        'run_condition_blind_embedding_build:InvalidSelectedAnchorWeights', ...
        'Selected anchor inclusion probabilities and weights must be finite and valid.');
    one = table();
    one.scale_index = s;
    one.chunk_sec = double(selected.chunk_sec(find(idx, 1)));
    one.anchor_manifest_mode = string(params.anchor_manifest_mode);
    one.n_rows = nnz(idx);
    one.n_base_rows = nnz(stage == "base_time_even");
    one.n_rare_enriched_rows = nnz(stage == "rare_strata_enriched");
    one.mean_base_inclusion_probability = mean(baseP, 'omitnan');
    one.mean_rare_inclusion_probability = mean(rareP, 'omitnan');
    one.mean_final_inclusion_probability = mean(finalP, 'omitnan');
    one.min_final_inclusion_probability = min(finalP, [], 'omitnan');
    one.max_audit_inverse_probability_weight = max(auditW, [], 'omitnan');
    one.audit_inverse_probability_effective_sample_size = local_effective_sample_size(auditW);
    one.pca_training_weight_min = min(pcaW, [], 'omitnan');
    one.pca_training_weight_max = max(pcaW, [], 'omitnan');
    one.graph_training_weight_min = min(graphW, [], 'omitnan');
    one.graph_training_weight_max = max(graphW, [], 'omitnan');
    one.weights_used_for_pca = logical(params.use_anchor_weights_for_pca);
    one.weights_used_for_graph_edges = false;
    one.labels_used_for_weight_audit = "none";
    one.arena_used_for_weight_audit = false;
    one.condition_used_for_weight_audit = false;
    Audit = [Audit; one]; %#ok<AGROW>
end
end

function Audit = local_rare_coverage_audit(selected, params)
if ismember('rare_stratum_id', selected.Properties.VariableNames)
    stratum = string(selected.rare_stratum_id);
else
    stratum = repmat("none", height(selected), 1);
end
if ismember('anchor_stage', selected.Properties.VariableNames)
    stage = string(selected.anchor_stage);
else
    stage = repmat("base_time_even", height(selected), 1);
end
Audit = table();
scales = unique(double(selected.scale_index), 'stable')';
for s = scales
    scaleMask = double(selected.scale_index) == s;
    ids = unique(stratum(scaleMask), 'stable');
    for j = 1:numel(ids)
        idx = scaleMask & stratum == ids(j);
        one = table();
        one.scale_index = s;
        one.chunk_sec = double(selected.chunk_sec(find(scaleMask, 1)));
        one.anchor_manifest_mode = string(params.anchor_manifest_mode);
        one.rare_stratum_id = ids(j);
        one.n_selected_rows = nnz(idx);
        one.selected_fraction_within_scale = nnz(idx) ./ nnz(scaleMask);
        one.n_base_rows = nnz(idx & stage == "base_time_even");
        one.n_rare_enriched_rows = nnz(idx & stage == "rare_strata_enriched");
        one.labels_used_for_coverage_audit = "none";
        one.arena_used_for_coverage_audit = false;
        one.condition_used_for_coverage_audit = false;
        Audit = [Audit; one]; %#ok<AGROW>
    end
end
end

function ess = local_effective_sample_size(w)
w = double(w(:));
w = w(isfinite(w) & w > 0);
if isempty(w)
    ess = NaN;
else
    ess = sum(w).^2 ./ sum(w.^2);
end
end

function T = local_feature_dictionary_header(featureDict, featureSessionAudit)
T = table();
T.embedding_feature_index = NaN;
T.scale_index = NaN;
T.chunk_sec = NaN;
T.obs_name = "canonical_run05_feature_count";
T.base_feature = "all_34_canonical_dyadic_features";
T.channel_type = "feature_dictionary_header";
T.summary_kind = "source_feature_dictionary";
T.temporal_bin = NaN;
T.dct_coefficient = NaN;
T.source_tensor = "derived/features_motif_discovery/feature_dictionary.csv";
if istable(featureSessionAudit) && ismember('valid_frame_fraction', featureSessionAudit.Properties.VariableNames) && ~isempty(featureSessionAudit)
    T.mean_frame_valid_fraction_for_channel = mean(featureSessionAudit.valid_frame_fraction, 'omitnan');
else
    T.mean_frame_valid_fraction_for_channel = NaN;
end
T.uses_frame_mask = true;
T.labels_used_for_summary = "none";
T.feature_family = "all";
T.unit = "mixed";
T.run07_feature_dictionary_role = "source_feature_dictionary_header";
T.labels_used_for_feature_construction = "none";
T.arena_used_for_feature_construction = false;
T.condition_used_for_feature_construction = false;
T.n_canonical_features = height(featureDict);
if istable(featureSessionAudit)
    T.n_feature_sessions = height(featureSessionAudit);
else
    T.n_feature_sessions = NaN;
end
end

function T = local_feature_dictionary_summary_rows(featureDictionary)
T = table();
families = unique(string(featureDictionary.feature_family), 'stable');
for i = 1:numel(families)
    idx = string(featureDictionary.feature_family) == families(i);
    one = local_feature_dictionary_header(featureDictionary, table());
    one.obs_name = "run07_embedding_feature_family_count";
    one.base_feature = families(i);
    one.channel_type = "feature_dictionary_summary";
    one.summary_kind = "feature_family_count";
    one.source_tensor = "embedding_feature_dictionary";
    one.feature_family = families(i);
    one.n_canonical_features = nnz(idx);
    one.n_feature_sessions = NaN;
    T = [T; one]; %#ok<AGROW>
end
end

function manifest = local_matrix_manifest(paths, Input, Embedding, Audit, inputRoot, params)
manifest = table();
F = local_active_run06_files(params);
for s = 1:numel(Input.scale)
    S = Input.scale(s);
    one = table();
    one.matrix_name = "scale_" + string(S.scale_index) + "_summary_matrix";
    one.matrix_role = "per_scale_pca_input";
    one.artifact_path = string(paths.summaryMatrixMat) + "#scale_index=" + string(S.scale_index);
    one.summary_profile = string(S.summary_profile);
    one.n_temporal_bins = double(S.n_temporal_bins);
    one.n_dct_coeffs = double(S.n_dct_coeffs);
    one.run06_expected_multiresolution_dims = double(S.n_summary_dims_run06);
    one.n_event_summary_dims = nnz(string(S.dimMeta.channel_type) == "event_summary");
    one.n_rows = size(S.X, 1);
    one.n_columns = size(S.X, 2);
    one.finite_value_fraction = mean(isfinite(S.X), 'all');
    one.source_csv = string(fullfile(inputRoot, char(F.anchor_manifest)));
    one.labels_used_for_matrix = "none";
    one.arena_used_for_matrix = false;
    one.condition_used_for_matrix = false;
    manifest = [manifest; one]; %#ok<AGROW>
end

globalOne = table();
globalOne.matrix_name = "global_ordinal_pc_stack";
globalOne.matrix_role = "global_pca_input";
globalOne.artifact_path = string(paths.embeddingModelMat) + "#global_matrix";
globalOne.summary_profile = "not_applicable_global_pc_stack";
globalOne.n_temporal_bins = NaN;
globalOne.n_dct_coeffs = NaN;
globalOne.run06_expected_multiresolution_dims = NaN;
globalOne.n_event_summary_dims = NaN;
globalOne.n_rows = size(Embedding.global_matrix, 1);
globalOne.n_columns = size(Embedding.global_matrix, 2);
globalOne.finite_value_fraction = Audit.global.finite_value_fraction;
globalOne.source_csv = string(paths.rowManifest);
globalOne.labels_used_for_matrix = "none";
globalOne.arena_used_for_matrix = false;
globalOne.condition_used_for_matrix = false;
manifest = [manifest; globalOne];

scoreOne = table();
scoreOne.matrix_name = "global_embedding_scores";
scoreOne.matrix_role = "downstream_condition_blind_embedding";
scoreOne.artifact_path = string(paths.globalScores);
scoreOne.summary_profile = "not_applicable_global_embedding_scores";
scoreOne.n_temporal_bins = NaN;
scoreOne.n_dct_coeffs = NaN;
scoreOne.run06_expected_multiresolution_dims = NaN;
scoreOne.n_event_summary_dims = NaN;
scoreOne.n_rows = size(Embedding.global_score, 1);
scoreOne.n_columns = size(Embedding.global_score, 2);
scoreOne.finite_value_fraction = mean(isfinite(Embedding.global_score), 'all');
scoreOne.source_csv = string(paths.globalScores);
scoreOne.labels_used_for_matrix = "none";
scoreOne.arena_used_for_matrix = false;
scoreOne.condition_used_for_matrix = false;
manifest = [manifest; scoreOne];
end

function local_ensure_dir(pathText)
if ~exist(pathText, 'dir')
    mkdir(pathText);
end
end

function T = local_read_csv(pathText)
opts = detectImportOptions(pathText, 'FileType', 'text', ...
    'Delimiter', ',', 'TextType', 'string');
T = readtable(pathText, opts);
end
