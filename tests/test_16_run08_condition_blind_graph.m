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
verifyGreaterThan(testCase, params.graph_n_pcs, 5);
verifyGreaterThan(testCase, params.k_neighbors, 5);
verifyTrue(testCase, any(params.k_sensitivity_values == params.k_neighbors));
verifyTrue(testCase, contains(params.graph_rule, "run_07 global PCs only"));
verifyTrue(testCase, contains(params.event_coverage_rule, "after graph construction"));
verifyTrue(testCase, contains(params.label_exclusion_rule, "no condition cohort"));
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

required = ["graph_parameter_audit.csv", "graph_input_manifest.csv", ...
    "graph_score_preprocess_audit.csv", "graph_node_manifest.csv", ...
    "graph_edge_list.csv", "graph_topology_summary.csv", ...
    "graph_degree_audit.csv", "graph_component_audit.csv", ...
    "graph_neighbor_composition_audit.csv", "graph_k_sensitivity_audit.csv", ...
    "graph_rare_event_coverage_audit.csv", "graph_event_node_audit.csv", ...
    "graph_qc_figure_manifest.csv"];
for i = 1:numel(required)
    verifyTrue(testCase, contains(runnerText, required(i)), required(i));
end
verifyTrue(testCase, contains(figureText, "graph_node_manifest.csv"));
verifyTrue(testCase, contains(figureText, "graph_degree_audit.csv"));
verifyTrue(testCase, contains(figureText, "graph_rare_event_coverage_audit.csv"));
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
