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
verifyTrue(testCase, params.allow_reviewed_snapshot_fallback);
verifyEqual(testCase, params.global_matrix_mode, "ordinal_pc_stack");
verifyEqual(testCase, params.scale_weight_mode, "equal_total_weight");
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
    "embedding_arena_sensitivity_audit.csv", "embedding_qc_figure_manifest.csv"];
for i = 1:numel(required)
    verifyTrue(testCase, contains(runnerText, required(i)), required(i));
end
verifyTrue(testCase, contains(figureText, "embedding_input_anchor_audit.csv"));
verifyTrue(testCase, contains(figureText, "embedding_pca_by_scale.csv"));
verifyTrue(testCase, contains(figureText, "embedding_global_scores.csv"));
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
verifyNotEqual(testCase, Embedding.scale(1).n_pcs_selected, Embedding.scale(2).n_pcs_selected);
verifyTrue(testCase, all(Audit.pcaByScale.labels_used_for_pca == "none"));
verifyFalse(testCase, any(Audit.scaleWeights.condition_used_for_scale_weight));
verifyFalse(testCase, ismember("condition_id", string(Embedding.scoreTable.Properties.VariableNames)));
verifyFalse(testCase, ismember("arena_label", string(Embedding.scoreTable.Properties.VariableNames)));
verifyTrue(testCase, ismember("arena_label", string(Embedding.rowManifest.Properties.VariableNames)));
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
    Scale(s).dimMeta = table();
    Scale(s).dimMeta.base_feature = "feature_" + string((1:d)');
    Scale(s).dimMeta.feature_family = repmat("mock", d, 1);
end
Input = struct();
Input.scale = Scale;
Input.rowManifest = [Scale(1).rowMeta; Scale(2).rowMeta];
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
