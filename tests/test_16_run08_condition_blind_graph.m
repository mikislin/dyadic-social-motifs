function tests = test_16_run08_condition_blind_graph
tests = functiontests(localfunctions);
end

function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));
end

function testRun08ConfigLoads(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
envNames = ["RUN08_GRAPH_RUN_MODE", "RUN08_GRAPH_OUTPUT_DIR", "RUN08_EMBEDDING_INPUT_DIR"];
oldValues = i_capture_env(envNames);
cleanup = onCleanup(@() i_restore_env(envNames, oldValues)); %#ok<NASGU>
i_clear_env(envNames);

params = load_multiscale_graph_config(fullfile(repoRoot, 'config', 'multiscale_graph_config.csv'));

verifyEqual(testCase, params.run_mode, "smoke");
verifyEqual(testCase, params.output_dir, "derived/graph_motif_discovery_smoke");
verifyEqual(testCase, params.embedding_input_dir, "derived/embedding_motif_discovery_smoke");
verifyEqual(testCase, params.anchor_manifest_mode, "primary");
verifyGreaterThan(testCase, params.graph_n_pcs, 5);
verifyGreaterThan(testCase, params.k_neighbors, 5);
verifyTrue(testCase, any(params.k_sensitivity_values == params.k_neighbors));
verifyTrue(testCase, contains(params.graph_rule, "run_07 global PCs only"));
verifyTrue(testCase, contains(params.event_coverage_rule, "after graph construction"));
verifyTrue(testCase, contains(params.label_exclusion_rule, "no condition cohort"));
verifyGreaterThanOrEqual(testCase, params.audit_resample_replicates, 1);
verifyGreaterThan(testCase, params.audit_resample_panel_max_nodes, params.k_neighbors + 1);
verifyTrue(testCase, params.audit_session_excluded_enabled);
verifyTrue(testCase, params.umap_enabled);
verifyGreaterThanOrEqual(testCase, params.umap_num_epochs, 10);
end

function testRun08FullEnvRoutesToProductionRoots(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
envNames = ["RUN08_GRAPH_RUN_MODE", "RUN08_GRAPH_OUTPUT_DIR", "RUN08_EMBEDDING_INPUT_DIR"];
oldValues = i_capture_env(envNames);
cleanup = onCleanup(@() i_restore_env(envNames, oldValues)); %#ok<NASGU>

setenv('RUN08_GRAPH_RUN_MODE', 'full');
setenv('RUN08_GRAPH_OUTPUT_DIR', 'derived/graph_motif_discovery');
setenv('RUN08_EMBEDDING_INPUT_DIR', 'derived/embedding_motif_discovery');

params = load_multiscale_graph_config(fullfile(repoRoot, 'config', 'multiscale_graph_config.csv'));

verifyEqual(testCase, params.run_mode, "full");
verifyEqual(testCase, params.output_dir, "derived/graph_motif_discovery");
verifyEqual(testCase, params.embedding_input_dir, "derived/embedding_motif_discovery");
verifyFalse(testCase, contains(params.output_dir, "smoke"));
verifyFalse(testCase, contains(params.embedding_input_dir, "smoke"));
end

function testRun08FullModeAutoRoutesCanonicalDefaults(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
envNames = ["RUN08_GRAPH_RUN_MODE", "RUN08_GRAPH_OUTPUT_DIR", "RUN08_EMBEDDING_INPUT_DIR"];
oldValues = i_capture_env(envNames);
cleanup = onCleanup(@() i_restore_env(envNames, oldValues)); %#ok<NASGU>

setenv('RUN08_GRAPH_RUN_MODE', 'full');
setenv('RUN08_GRAPH_OUTPUT_DIR', '');
setenv('RUN08_EMBEDDING_INPUT_DIR', '');

params = load_multiscale_graph_config(fullfile(repoRoot, 'config', 'multiscale_graph_config.csv'));

verifyEqual(testCase, params.run_mode, "full");
verifyEqual(testCase, params.output_dir, "derived/graph_motif_discovery");
verifyEqual(testCase, params.embedding_input_dir, "derived/embedding_motif_discovery");
end

function testRun08RejectsExplicitFullModeWithSmokeRoots(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
envNames = ["RUN08_GRAPH_RUN_MODE", "RUN08_GRAPH_OUTPUT_DIR", "RUN08_EMBEDDING_INPUT_DIR"];
oldValues = i_capture_env(envNames);
cleanup = onCleanup(@() i_restore_env(envNames, oldValues)); %#ok<NASGU>

setenv('RUN08_GRAPH_RUN_MODE', 'full');
setenv('RUN08_GRAPH_OUTPUT_DIR', 'derived/graph_motif_discovery_smoke');
setenv('RUN08_EMBEDDING_INPUT_DIR', 'derived/embedding_motif_discovery');

verifyError(testCase, ...
    @() load_multiscale_graph_config(fullfile(repoRoot, 'config', 'multiscale_graph_config.csv')), ...
    'load_multiscale_graph_config:BadFullOutputDir');
end

function testRun08PaperScriptHasNoLocalFunctions(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
scriptPath = fullfile(repoRoot, 'paper', 'run_08_build_condition_blind_motif_graph.m');
txt = string(fileread(scriptPath));

verifyFalse(testCase, contains(txt, newline + "function "));
verifyTrue(testCase, contains(txt, "run_condition_blind_motif_graph_build"));
end

function testRun08DeclaresRequiredOutputs(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
runnerPath = fullfile(repoRoot, 'graph', 'run_condition_blind_motif_graph_build.m');
figurePath = fullfile(repoRoot, 'graph', 'make_run08_graph_qc_figures.m');
runnerText = string(fileread(runnerPath));
figureText = string(fileread(figurePath));
sensitivityText = string(fileread(fullfile(repoRoot, 'graph', ...
    'build_condition_blind_graph_sensitivity_audits.m')));
visualizationText = string(fileread(fullfile(repoRoot, 'graph', ...
    'build_run08_embedding_visualization_audits.m')));
rareStageText = string(fileread(fullfile(repoRoot, 'graph', ...
    'persist_run08_rare_definition_stages.m')));

required = ["graph_parameter_audit.csv", "graph_input_manifest.csv", ...
    "graph_score_preprocess_audit.csv", "graph_node_manifest.csv", ...
    "graph_edge_list.csv", "graph_topology_summary.csv", ...
    "graph_degree_audit.csv", "graph_component_audit.csv", ...
    "graph_neighbor_composition_audit.csv", "graph_k_sensitivity_audit.csv", ...
    "graph_rare_event_coverage_audit.csv", "graph_event_node_audit.csv", ...
    "graph_rare_strata_composition_audit.csv", ...
    "graph_rare_strata_neighbor_retention_audit.csv", ...
    "graph_baseline_vs_enriched_coverage_audit.csv", ...
    "graph_scale_mixing_matrix_audit.csv", ...
    "graph_neighbor_null_normalized_audit.csv", ...
    "graph_event_prevalence_fold_audit.csv", ...
    "graph_qc_figure_manifest.csv"];
for i = 1:numel(required)
    verifyTrue(testCase, contains(runnerText, required(i)), required(i));
end
verifyTrue(testCase, contains(figureText, "graph_node_manifest.csv"));
verifyTrue(testCase, contains(figureText, "graph_degree_audit.csv"));
verifyTrue(testCase, contains(figureText, "graph_rare_event_coverage_audit.csv"));
for required = ["graph_anchor_stage_sensitivity_audit.csv", ...
        "graph_primary_base_only_edge_list.csv", ...
        "graph_session_excluded_sensitivity_audit.csv", ...
        "graph_session_excluded_edge_list_audit.csv", ...
        "graph_neighborhood_resampling_audit.csv"]
    verifyTrue(testCase, contains(sensitivityText, required), required);
end
for required = ["graph_global_pca_cumulative_variance_audit.csv", ...
        "graph_umap_embedding_audit.csv", "visualization_only_not_graph_or_motif_input"]
    verifyTrue(testCase, contains(visualizationText, required), required);
end
for required = ["rare_strata_selection_definition_locked.csv", ...
        "rare_strata_postfit_definition.csv", ...
        "rare_strata_definition_provenance_audit.csv"]
    verifyTrue(testCase, contains(rareStageText, required), required);
end
end

function testConditionBlindKnnGraphMock(testCase)
rng(16);
n = 18;
X = randn(n, 4);
rowMeta = table();
rowMeta.embedding_row_id = (101:(100 + n))';
rowMeta.scale_index = repelem([1; 2; 3], 6, 1);
rowMeta.chunk_sec = repelem([0.2; 2.1; 13.2], 6, 1);

params = struct();
params.k_neighbors = 5;
params.knn_block_size = 16;
Graph = build_condition_blind_knn_graph(X, rowMeta, params);

verifyEqual(testCase, Graph.k, 5);
verifyEqual(testCase, height(Graph.Edges), n * params.k_neighbors);
verifyFalse(testCase, any(Graph.Edges.source_node_id == Graph.Edges.target_node_id));
verifyTrue(testCase, all(isfinite(Graph.Edges.neighbor_distance)));
verifyTrue(testCase, all(Graph.Edges.labels_used_for_graph == "none"));
verifyFalse(testCase, any(Graph.Edges.arena_used_for_graph));
verifyFalse(testCase, any(Graph.Edges.condition_used_for_graph));
verifyTrue(testCase, all(Graph.Edges.neighbor_rank >= 1 & ...
    Graph.Edges.neighbor_rank <= params.k_neighbors));
end

function testDemoDeclaresConditionBlindPrinciples(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
demoPath = fullfile(repoRoot, 'paper', 'demo_run08_condition_blind_graph_principles.m');
txt = string(fileread(demoPath));

verifyTrue(testCase, contains(txt, "build_condition_blind_knn_graph"));
verifyTrue(testCase, contains(txt, "event_flags_used_for_demo_graph = false"));
verifyTrue(testCase, contains(txt, 'labels_used_for_graph = "none"'));
verifyTrue(testCase, contains(txt, "demo_run08_event_coverage_by_scale.csv"));
end

function testScaleMixingAndNullAuditMock(testCase)
rng(116);
n = 24;
X = randn(n, 5);
Meta = table();
Meta.graph_node_id = (1:n)';
Meta.embedding_row_id = (101:(100 + n))';
Meta.scale_index = repelem((1:3)', 8);
Meta.chunk_sec = repelem([0.2; 1.2; 7.9], 8);
Meta.session_index = repmat(repelem((1:4)', 2), 3, 1);
params = struct('k_neighbors', 4, 'knn_block_size', 16);
G = build_condition_blind_knn_graph(X, Meta, params);
arena = repmat(["a";"b"], n / 2, 1);
[M, N] = compute_graph_scale_null_audits(G, Meta, arena);

verifyEqual(testCase, height(M), 9);
verifyEqual(testCase, numel(unique(M.source_scale_index)), 3);
verifyEqual(testCase, numel(unique(M.target_scale_index)), 3);
verifyTrue(testCase, all(isfinite(M.observed_over_random_ratio)));
verifyEqual(testCase, height(N), 12);
verifyTrue(testCase, all(isfinite(N.observed_over_random_ratio)));
verifyTrue(testCase, all(N.labels_used_for_graph == "none"));
verifyFalse(testCase, any(N.condition_used_for_audit));
end

function testSessionExcludedGraphIsAuditOnlyAndLabelInvariant(testCase)
rng(117);
n = 32;
X = randn(n, 6);
Meta = table();
Meta.embedding_row_id = (1:n)';
Meta.session_index = repelem((1:4)', 8);
Meta.condition_id = "condition_" + string(randperm(n))';
params = struct('k_neighbors', 5, 'knn_block_size', 16);
G1 = build_condition_blind_session_excluded_knn_audit(X, Meta, params);
Meta.condition_id = flipud(Meta.condition_id);
Meta.arena_label = "arena_" + string(randperm(n))';
G2 = build_condition_blind_session_excluded_knn_audit(X, Meta, params);

sourceSession = Meta.session_index(G1.Edges.source_node_id);
targetSession = Meta.session_index(G1.Edges.target_node_id);
verifyFalse(testCase, any(sourceSession == targetSession));
verifyEqual(testCase, G1.Edges.target_node_id, G2.Edges.target_node_id);
verifyEqual(testCase, G1.Edges.neighbor_distance, G2.Edges.neighbor_distance, 'AbsTol', 1e-12);
verifyTrue(testCase, all(G1.Edges.audit_only_not_primary_graph));
verifyTrue(testCase, all(G1.Edges.labels_used_for_graph_distance == "none"));
end

function testPrevalenceFoldAuditIsFinite(testCase)
C = table();
C.event_id = ["contact";"contact"];
C.event_rule = ["gt_zero";"gt_zero"];
C.scale_index = [1;2];
C.chunk_sec = [0.2;1.2];
C.n_graph_nodes = [100;100];
C.n_graph_event_nodes = [20;0];
C.baseline_n_graph_nodes = [25;25];
C.baseline_n_graph_event_nodes = [0;0];
A = compute_graph_event_prevalence_fold_audit(C);
verifyEqual(testCase, height(A), 3);
verifyTrue(testCase, all(isfinite(A.prevalence_fold_enriched_over_baseline)));
verifyTrue(testCase, all(isfinite(A.log2_prevalence_fold)));
verifyTrue(testCase, all(A.half_count_smoothing_used));
verifyFalse(testCase, any(A.condition_used_for_prevalence_audit));
end

function testNewAuditDistanceCodeDoesNotReadForbiddenLabels(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
files = ["build_condition_blind_graph_sensitivity_audits.m", ...
    "build_condition_blind_session_excluded_knn_audit.m", ...
    "build_run08_embedding_visualization_audits.m"];
txt = "";
for file = files
    txt = txt + newline + string(fileread(fullfile(repoRoot, 'graph', file)));
end
for token = [".condition_id", ".condition_label", ".cohort_id", ...
        ".arena_label", ".drug", ".genotype", ".outcome", ".rare_stratum_id"]
    verifyFalse(testCase, contains(lower(txt), token), ...
        "Distance/UMAP audit code reads a forbidden label-like value: " + token);
end
verifyTrue(testCase, contains(txt, "labels_used_for_distance"));
verifyTrue(testCase, contains(txt, "visualization_only_not_graph_or_motif_input"));
end

function testLegacyConditionAwareTokensNotPorted(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
files = [
    string(fullfile(repoRoot, 'paper', 'run_08_build_condition_blind_motif_graph.m'))
    string(fullfile(repoRoot, 'graph', 'run_condition_blind_motif_graph_build.m'))
    string(fullfile(repoRoot, 'graph', 'build_condition_blind_knn_graph.m'))
    ];
txt = "";
for i = 1:numel(files)
    txt = txt + newline + string(fileread(files(i)));
end

forbidden = ["conditionQuota", "minAnchorsPerCondition", "rare_aware", ...
    "rare_weighted", "is_training", "is_validation", "holdout_fold", ...
    "fitgmdist", "kmeans(", "motif_family_id", "condition_motif"];
for i = 1:numel(forbidden)
    verifyFalse(testCase, contains(txt, forbidden(i)), ...
        "Forbidden legacy token was ported into run_08 graph layer: " + forbidden(i));
end
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
