function tests = test_15_run07_condition_blind_embedding
tests = functiontests(localfunctions);
end

function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));
end

function testRun07ConfigLoads(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
envNames = ["RUN07_EMBEDDING_RUN_MODE", "RUN07_EMBEDDING_OUTPUT_DIR", ...
    "RUN07_CHUNK_INPUT_DIR", "RUN07_FALLBACK_CHUNK_INPUT_DIR"];
oldValues = i_capture_env(envNames);
cleanup = onCleanup(@() i_restore_env(envNames, oldValues)); %#ok<NASGU>
i_clear_env(envNames);

params = load_multiscale_embedding_config(fullfile(repoRoot, 'config', 'multiscale_embedding_config.csv'));

verifyEqual(testCase, params.run_mode, "smoke");
verifyEqual(testCase, params.output_dir, "derived/embedding_motif_discovery_smoke");
verifyEqual(testCase, params.chunk_input_dir, "derived/chunks_motif_discovery");
verifyEqual(testCase, params.anchor_manifest_mode, "primary");
verifyFalse(testCase, params.use_anchor_weights_for_pca);
verifyTrue(testCase, params.allow_reviewed_snapshot_fallback);
verifyEqual(testCase, params.global_matrix_mode, "ordinal_pc_stack");
verifyEqual(testCase, params.scale_weight_mode, "equal_total_weight");
verifyTrue(testCase, params.preprocess_sparse_scale_safeguard_enabled);
verifyEqual(testCase, params.preprocess_min_iqr_to_std_ratio, 0.15, 'AbsTol', 1e-12);
verifyGreaterThan(testCase, params.preprocess_severe_tail_audit_abs_threshold, ...
    params.preprocess_tail_audit_abs_threshold);
verifyGreaterThan(testCase, params.micro_max_pcs, 5);
verifyGreaterThan(testCase, params.motif_max_pcs, params.micro_max_pcs);
verifyGreaterThan(testCase, params.context_max_pcs, params.micro_max_pcs);
verifyTrue(testCase, contains(params.embedding_dimension_rule, "embedding_dimension_audit"));
verifyTrue(testCase, contains(params.label_exclusion_rule, "no condition cohort"));
end

function testRun07FullModeAutoRoutesOutputDefault(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
envNames = ["RUN07_EMBEDDING_RUN_MODE", "RUN07_EMBEDDING_OUTPUT_DIR", ...
    "RUN07_CHUNK_INPUT_DIR", "RUN07_FALLBACK_CHUNK_INPUT_DIR"];
oldValues = i_capture_env(envNames);
cleanup = onCleanup(@() i_restore_env(envNames, oldValues)); %#ok<NASGU>

setenv('RUN07_EMBEDDING_RUN_MODE', 'full');
setenv('RUN07_EMBEDDING_OUTPUT_DIR', '');
setenv('RUN07_CHUNK_INPUT_DIR', '');
setenv('RUN07_FALLBACK_CHUNK_INPUT_DIR', '');

params = load_multiscale_embedding_config(fullfile(repoRoot, 'config', 'multiscale_embedding_config.csv'));

verifyEqual(testCase, params.run_mode, "full");
verifyEqual(testCase, params.output_dir, "derived/embedding_motif_discovery");
verifyEqual(testCase, params.chunk_input_dir, "derived/chunks_motif_discovery");
verifyFalse(testCase, contains(params.output_dir, "smoke"));
end

function testRun07RejectsExplicitFullModeWithSmokeRoots(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
envNames = ["RUN07_EMBEDDING_RUN_MODE", "RUN07_EMBEDDING_OUTPUT_DIR", ...
    "RUN07_CHUNK_INPUT_DIR", "RUN07_FALLBACK_CHUNK_INPUT_DIR"];
oldValues = i_capture_env(envNames);
cleanup = onCleanup(@() i_restore_env(envNames, oldValues)); %#ok<NASGU>

setenv('RUN07_EMBEDDING_RUN_MODE', 'full');
setenv('RUN07_EMBEDDING_OUTPUT_DIR', 'derived/embedding_motif_discovery_smoke');
setenv('RUN07_CHUNK_INPUT_DIR', 'derived/chunks_motif_discovery');
setenv('RUN07_FALLBACK_CHUNK_INPUT_DIR', 'derived/chunks_motif_discovery_scale_survey_multires_v1');

verifyError(testCase, ...
    @() load_multiscale_embedding_config(fullfile(repoRoot, 'config', 'multiscale_embedding_config.csv')), ...
    'load_multiscale_embedding_config:BadFullOutputDir');
end

function testRun07PaperScriptHasNoLocalFunctions(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
scriptPath = fullfile(repoRoot, 'paper', 'run_07_build_condition_blind_embedding.m');
txt = string(fileread(scriptPath));

verifyFalse(testCase, contains(txt, newline + "function "));
verifyTrue(testCase, contains(txt, "run_condition_blind_embedding_build"));
end

function testRun07DeclaresRequiredOutputs(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
runnerPath = fullfile(repoRoot, 'embedding', 'run_condition_blind_embedding_build.m');
figurePath = fullfile(repoRoot, 'embedding', 'make_run07_embedding_qc_figures.m');
runnerText = string(fileread(runnerPath));
figureText = string(fileread(figurePath));

required = ["embedding_parameter_audit.csv", "embedding_input_scale_audit.csv", ...
    "embedding_input_anchor_audit.csv", "embedding_feature_dictionary.csv", ...
    "embedding_matrix_manifest.csv", "embedding_pca_by_scale.csv", ...
    "embedding_scale_weight_audit.csv", "embedding_stability_audit.csv", ...
    "embedding_anchor_weight_audit.csv", "embedding_rare_strata_coverage_audit.csv", ...
    "embedding_preprocess_dimension_audit.csv", ...
    "embedding_anchor_stage_pca_sensitivity_audit.csv", ...
    "embedding_arena_sensitivity_audit.csv", "embedding_qc_figure_manifest.csv"];
for i = 1:numel(required)
    verifyTrue(testCase, contains(runnerText, required(i)), required(i));
end
verifyTrue(testCase, contains(figureText, "embedding_input_anchor_audit.csv"));
verifyTrue(testCase, contains(figureText, "embedding_pca_by_scale.csv"));
verifyTrue(testCase, contains(figureText, "embedding_global_scores.csv"));
end

function testRun07UsesRun06SummaryProfileAuditForMatrixConstruction(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
builderPath = fullfile(repoRoot, 'embedding', 'build_primary_scale_embedding_inputs.m');
runnerPath = fullfile(repoRoot, 'embedding', 'run_condition_blind_embedding_build.m');
builderText = string(fileread(builderPath));
runnerText = string(fileread(runnerPath));

verifyTrue(testCase, contains(runnerText, "featureDict, dimensionAudit, params"));
verifyTrue(testCase, contains(builderText, "local_summary_profile_for_scale"));
verifyTrue(testCase, contains(builderText, "run06_embedding_dimension_audit.csv"));
verifyTrue(testCase, contains(builderText, "BadMotifSummaryProfile"));
verifyTrue(testCase, contains(builderText, "SummaryDimensionMismatch"));
verifyTrue(testCase, contains(builderText, "profile.n_temporal_bins"));
verifyTrue(testCase, contains(builderText, "profile.n_dct_coeffs"));
verifyFalse(testCase, contains(builderText, "'nTemporalBins', params.summary_temporal_bins"));
verifyFalse(testCase, contains(builderText, "'nDctCoeffs', params.summary_dct_coeffs"));
end

function testRun07CurrentRun06AuditHasBandAdaptiveSummaryProfiles(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
auditPath = fullfile(repoRoot, 'derived', 'chunks_motif_discovery', 'embedding_dimension_audit.csv');
primaryPath = fullfile(repoRoot, 'derived', 'chunks_motif_discovery', 'primary_operational_scales.csv');
if ~(isfile(auditPath) && isfile(primaryPath))
    runnerText = string(fileread(fullfile(repoRoot, 'embedding', 'run_condition_blind_embedding_build.m')));
    verifyTrue(testCase, contains(runnerText, "derived/chunks_motif_discovery"));
    return
end

Audit = readtable(auditPath, 'TextType', 'string');
Primary = readtable(primaryPath, 'TextType', 'string');
Audit = Audit(ismember(Audit.scale_index, Primary.scale_index), :);
verifyFalse(testCase, isempty(Audit));

motif = Audit(string(Audit.initial_band) == "motif", :);
base = Audit(string(Audit.initial_band) == "micro" | string(Audit.initial_band) == "context", :);
verifyFalse(testCase, isempty(motif));
verifyFalse(testCase, isempty(base));

verifyTrue(testCase, all(motif.n_temporal_bins == 12));
verifyTrue(testCase, all(motif.n_dct_coeffs == 8));
verifyTrue(testCase, all(base.n_temporal_bins == 6));
verifyTrue(testCase, all(base.n_dct_coeffs == 4));
verifyTrue(testCase, all(contains(string(motif.summary_profile), "bins12_dct8")));
verifyTrue(testCase, all(contains(string(base.summary_profile), "bins6_dct4")));

staleMotifDims = i_expected_summary_dims(motif.n_input_channels(1), ...
    motif.n_boolean_transition_channels(1), 6, 4);
verifyTrue(testCase, all(motif.n_summary_dims ~= staleMotifDims), ...
    "Motif summaries must not collapse to stale global 6-bin/4-DCT dimensions.");
verifyGreaterThan(testCase, median(motif.n_summary_dims), median(base.n_summary_dims));
end

function testRun07EventSummaryChannelsRemainSeparateAndAuditable(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
builderText = string(fileread(fullfile(repoRoot, 'embedding', 'build_primary_scale_embedding_inputs.m')));
runnerText = string(fileread(fullfile(repoRoot, 'embedding', 'run_condition_blind_embedding_build.m')));

verifyTrue(testCase, contains(builderText, 'channel_type = repmat("event_summary"'));
verifyTrue(testCase, contains(builderText, 'summary_kind = repmat("condition_blind_event_summary"'));
verifyTrue(testCase, contains(builderText, 'event_summary_channels_separate'));
verifyTrue(testCase, contains(runnerText, "n_event_summary_dims"));
verifyTrue(testCase, contains(runnerText, "run06_expected_multiresolution_dims"));
end

function testPcaUsesScaleSpecificGuidance(testCase)
params = i_mock_params();
Input = i_mock_input();
dimensionAudit = table( ...
    [101; 202], [0.2; 8.0], ["micro"; "context"], [3; 7], [5; 9], ...
    'VariableNames', {'scale_index','chunk_sec','initial_band', ...
    'recommended_pcs_for_next_embedding','n_pcs_90pct'});

[Embedding, Audit] = fit_condition_blind_embedding_pca(Input, dimensionAudit, params);

verifyEqual(testCase, Audit.pcaByScale.n_pcs_selected(:), [3; 7]);
for s = 1:numel(Embedding.scale)
    nPC = Embedding.scale(s).n_pcs_selected;
    expected = sum(Embedding.scale(s).explained(1:nPC));
    verifyEqual(testCase, Audit.pcaByScale.cum_selected_explained(s), expected, 'AbsTol', 1e-10);
end
verifyNotEqual(testCase, Embedding.scale(1).n_pcs_selected, Embedding.scale(2).n_pcs_selected);
verifyTrue(testCase, all(Audit.pcaByScale.labels_used_for_pca == "none"));
verifyTrue(testCase, all(Audit.preprocessDimensions.labels_used_for_preprocessing == "none"));
verifyFalse(testCase, any(Audit.preprocessDimensions.condition_used_for_preprocessing));
verifyTrue(testCase, all(Audit.anchorStageSensitivity.audit_status == "complete"));
verifyFalse(testCase, any(Audit.anchorStageSensitivity.anchor_stage_used_for_primary_pca));
verifyFalse(testCase, any(Audit.scaleWeights.condition_used_for_scale_weight));
verifyFalse(testCase, ismember("condition_id", string(Embedding.scoreTable.Properties.VariableNames)));
verifyFalse(testCase, ismember("arena_label", string(Embedding.scoreTable.Properties.VariableNames)));
verifyTrue(testCase, ismember("arena_label", string(Embedding.rowManifest.Properties.VariableNames)));
end

function testSparseFeatureScaleSafeguardAndTailAudit(testCase)
params = i_mock_params();
Input = i_sparse_guard_input();
dimensionAudit = table(101, 1.25, "motif", 5, 8, ...
    'VariableNames', {'scale_index','chunk_sec','initial_band', ...
    'recommended_pcs_for_next_embedding','n_pcs_90pct'});

[~, Audit] = fit_condition_blind_embedding_pca(Input, dimensionAudit, params);
row = Audit.preprocessDimensions(Audit.preprocessDimensions.input_dimension_index == 1, :);
verifyEqual(testCase, height(row), 1);
verifyTrue(testCase, row.sparse_scale_safeguard_triggered);
verifyEqual(testCase, row.scale_method, "winsorized_std_sparse_guard");
verifyLessThan(testCase, row.iqr_to_std_ratio, params.preprocess_min_iqr_to_std_ratio);
verifyLessThan(testCase, abs(row.standardized_std - 1), 1e-10);
verifyEqual(testCase, row.n_abs_gt_severe_tail_threshold, 0);
verifyLessThan(testCase, Audit.pcaByScale.pc1_top_dimension_loading_fraction, 0.98);
end

function testConditionAndArenaMetadataCannotChangePcaFit(testCase)
params = i_mock_params();
InputA = i_mock_input();
InputB = InputA;
for s = 1:numel(InputB.scale)
    n = height(InputB.scale(s).rowMeta);
    InputB.scale(s).rowMeta.condition_id = "changed_condition_" + string((1:n)');
    InputB.scale(s).rowMeta.arena_label = repmat("changed_arena", n, 1);
end
dimensionAudit = table( ...
    [101; 202], [0.2; 8.0], ["micro"; "context"], [3; 7], [5; 9], ...
    'VariableNames', {'scale_index','chunk_sec','initial_band', ...
    'recommended_pcs_for_next_embedding','n_pcs_90pct'});

[EmbeddingA, ~] = fit_condition_blind_embedding_pca(InputA, dimensionAudit, params);
[EmbeddingB, ~] = fit_condition_blind_embedding_pca(InputB, dimensionAudit, params);
for s = 1:numel(EmbeddingA.scale)
    verifyEqual(testCase, abs(EmbeddingA.scale(s).coeff), ...
        abs(EmbeddingB.scale(s).coeff), 'AbsTol', 1e-10);
end
verifyEqual(testCase, abs(EmbeddingA.global_score), abs(EmbeddingB.global_score), 'AbsTol', 1e-10);
end

function testLegacyConditionAwareTokensNotPorted(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
files = [
    string(fullfile(repoRoot, 'paper', 'run_07_build_condition_blind_embedding.m'))
    string(fullfile(repoRoot, 'embedding', 'run_condition_blind_embedding_build.m'))
    string(fullfile(repoRoot, 'embedding', 'build_primary_scale_embedding_inputs.m'))
    string(fullfile(repoRoot, 'embedding', 'fit_condition_blind_embedding_pca.m'))
    ];
txt = "";
for i = 1:numel(files)
    txt = txt + newline + string(fileread(files(i)));
end

forbidden = ["conditionQuota", "minAnchorsPerCondition", "rare_aware", ...
    "rare_weighted", "is_training", "is_validation", "holdout_fold", ...
    "cluster(", "fitgmdist", "kmeans("];
for i = 1:numel(forbidden)
    verifyFalse(testCase, contains(txt, forbidden(i)), ...
        "Forbidden legacy token was ported into run_07 embedding: " + forbidden(i));
end
end

function params = i_mock_params()
params = struct();
params.rng_seed = 17;
params.preprocess_winsor_low = 0.01;
params.preprocess_winsor_high = 0.99;
params.preprocess_min_robust_scale = 1e-8;
params.preprocess_sparse_scale_safeguard_enabled = true;
params.preprocess_min_iqr_to_std_ratio = 0.15;
params.preprocess_tail_audit_abs_threshold = 10;
params.preprocess_severe_tail_audit_abs_threshold = 100;
params.min_pcs_per_scale = 2;
params.micro_max_pcs = 6;
params.motif_max_pcs = 8;
params.context_max_pcs = 10;
params.global_n_pcs = 5;
params.global_matrix_mode = "ordinal_pc_stack";
params.standardize_scale_scores_before_global = true;
params.scale_score_winsor_abs = 8;
params.scale_weight_mode = "equal_total_weight";
params.scale_weight_rule = "test_equal_weight";
params.run_mode = "smoke";
params.smoke_stability_splits = 2;
params.stability_splits = 2;
params.stability_n_pcs_compared = 3;
params.enrichment_sensitivity_n_pcs_compared = 3;
params.arena_sensitivity_n_pcs = 2;
end

function Input = i_mock_input()
rng(2);
Scale = repmat(struct('scale_index', [], 'chunk_sec', [], 'hierarchical_role', "", ...
    'rowMeta', table(), 'X', [], 'dimMeta', table(), 'summaryAudit', table()), 2, 1);
scaleIndex = [101, 202];
chunkSec = [0.2, 8.0];
role = ["micro", "context"];
for s = 1:2
    n = 40;
    d = 12;
    Scale(s).scale_index = scaleIndex(s);
    Scale(s).chunk_sec = chunkSec(s);
    Scale(s).hierarchical_role = role(s);
    Scale(s).X = randn(n, d) + 0.1 .* (1:n)' .* randn(1, d);
    Scale(s).rowMeta = table();
    Scale(s).rowMeta.embedding_row_id = ((s - 1) * n + (1:n))';
    Scale(s).rowMeta.scale_index = repmat(Scale(s).scale_index, n, 1);
    Scale(s).rowMeta.chunk_sec = repmat(Scale(s).chunk_sec, n, 1);
    Scale(s).rowMeta.arena_label = repmat(["big"; "small"], n / 2, 1);
    Scale(s).rowMeta.condition_id = repmat("provenance_only", n, 1);
    Scale(s).rowMeta.anchor_stage = [repmat("base_time_even", n / 2, 1); ...
        repmat("rare_strata_enriched", n / 2, 1)];
    Scale(s).dimMeta = table();
    Scale(s).dimMeta.embedding_feature_index = (1:d)';
    Scale(s).dimMeta.base_feature = "feature_" + string((1:d)');
    Scale(s).dimMeta.feature_family = repmat("mock", d, 1);
    Scale(s).dimMeta.summary_kind = repmat("mock_summary", d, 1);
    Scale(s).dimMeta.dct_coefficient = nan(d, 1);
end
Input = struct();
Input.scale = Scale;
Input.rowManifest = [Scale(1).rowMeta; Scale(2).rowMeta];
Input.featureDictionary = table();
Input.matrixAudit = table();
end

function Input = i_sparse_guard_input()
rng(31);
n = 400;
d = 16;
Scale = struct();
Scale.scale_index = 101;
Scale.chunk_sec = 1.25;
Scale.hierarchical_role = "motif";
Scale.X = randn(n, d);
Scale.X(:, 1) = [repmat(-0.0077, 100, 1); zeros(220, 1); ones(80, 1)];
Scale.rowMeta = table();
Scale.rowMeta.embedding_row_id = (1:n)';
Scale.rowMeta.scale_index = repmat(101, n, 1);
Scale.rowMeta.chunk_sec = repmat(1.25, n, 1);
Scale.rowMeta.arena_label = repmat(["big"; "small"], n / 2, 1);
Scale.rowMeta.condition_id = repmat("provenance_only", n, 1);
Scale.rowMeta.anchor_stage = [repmat("base_time_even", n / 2, 1); ...
    repmat("rare_strata_enriched", n / 2, 1)];
Scale.dimMeta = table();
Scale.dimMeta.embedding_feature_index = (1:d)';
Scale.dimMeta.base_feature = "feature_" + string((1:d)');
Scale.dimMeta.feature_family = repmat("mock", d, 1);
Scale.dimMeta.summary_kind = repmat("low_frequency_dct", d, 1);
Scale.dimMeta.dct_coefficient = (1:d)';
Input = struct();
Input.scale = Scale;
Input.rowManifest = Scale.rowMeta;
Input.featureDictionary = table();
Input.matrixAudit = table();
end

function i_restore_env(envNames, oldValues)
for i = 1:numel(envNames)
    setenv(char(envNames(i)), char(oldValues(i)));
end
end

function i_clear_env(envNames)
for i = 1:numel(envNames)
    setenv(char(envNames(i)), '');
end
end

function values = i_capture_env(envNames)
values = strings(size(envNames));
for i = 1:numel(envNames)
    values(i) = string(getenv(envNames(i)));
end
end

function nDims = i_expected_summary_dims(nChannels, nBool, nBins, nDct)
nBase = 7;
nDims = double(nChannels) .* (nBase + double(nBins) + double(nDct)) + double(nBool);
end
