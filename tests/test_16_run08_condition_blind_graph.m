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
verifyTrue(testCase, params.consensus_enabled);
verifyGreaterThanOrEqual(testCase, params.consensus_replicates, 4);
verifyTrue(testCase, any(params.consensus_k_values == params.k_neighbors));
verifyEqual(testCase, params.consensus_support_thresholds, [0.3 0.5 0.7]);
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
consensusText = string(fileread(fullfile(repoRoot, 'graph', ...
    'build_condition_blind_consensus_neighborhood.m')));

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
for required = ["graph_dimension_resample_manifest.csv", ...
        "graph_stage_balanced_resample_manifest.csv", ...
        "graph_replicate_manifest.csv", ...
        "graph_resampled_edge_support_audit.csv", ...
        "graph_consensus_edge_list.csv", ...
        "graph_consensus_node_manifest.csv", ...
        "graph_consensus_node_stability_audit.csv", ...
        "graph_consensus_threshold_sensitivity.csv", ...
        "graph_consensus_topology_summary.csv", ...
        "graph_consensus_scale_mixing_audit.csv", ...
        "graph_consensus_session_sensitivity_audit.csv", ...
        "graph_consensus_rare_event_coverage_audit.csv", ...
        "graph_consensus_qc_figure_manifest.csv", ...
        "run08_to_run09_node_input.csv", ...
        "run08_to_run09_edge_list.csv", ...
        "run08_to_run09_handoff_manifest.csv"]
    verifyTrue(testCase, contains(consensusText, required), required);
end
end

function testConsensusNeighborhoodIsFiniteUniqueAndLabelInvariant(testCase)
rng(118);
n = 96;
d = 24;
X = randn(n, d);
Meta = table();
Meta.graph_node_id = (1:n)';
Meta.embedding_row_id = (1001:(1000 + n))';
Meta.scale_index = repelem((1:3)', n / 3);
Meta.chunk_sec = repelem([0.3;1.3;7.9], n / 3);
Meta.session_index = repmat(repelem((1:8)', 4), 3, 1);
Meta.anchor_stage = repmat("rare_strata_enriched", n, 1);
for s = 1:3
    idx = find(Meta.scale_index == s, 8, 'first');
    Meta.anchor_stage(idx) = "base_time_even";
end
Meta.condition_id = "condition_" + string(randperm(n))';
Meta.arena_label = "arena_" + string(mod((1:n)', 2));

baseParams = struct('k_neighbors', 5, 'knn_block_size', 32);
Primary = build_condition_blind_knn_graph(X, Meta, baseParams);
Event = table();
Event.graph_node_id = (1:n)';
Event.embedding_row_id = Meta.embedding_row_id;
Event.scale_index = Meta.scale_index;
Event.chunk_sec = Meta.chunk_sec;
Event.labels_used_for_event_node_audit = repmat("none", n, 1);
Event.arena_used_for_event_node_audit = false(n, 1);
Event.condition_used_for_event_node_audit = false(n, 1);
Event.contact_present = mod((1:n)', 7) == 0;

params = i_consensus_test_params();
root1 = string(tempname); mkdir(root1); mkdir(fullfile(root1, 'figures'));
root2 = string(tempname); mkdir(root2); mkdir(fullfile(root2, 'figures'));
cleanup = onCleanup(@() i_remove_dirs([root1 root2])); %#ok<NASGU>
A1 = build_condition_blind_consensus_neighborhood(X, Meta, Primary, Event, params, root1);
Meta.condition_id = flipud(Meta.condition_id);
Meta.arena_label = flipud(Meta.arena_label);
A2 = build_condition_blind_consensus_neighborhood(X, Meta, Primary, Event, params, root2);

E1 = readtable(fullfile(root1, 'graph_consensus_edge_list.csv'), 'TextType', 'string');
E2 = readtable(fullfile(root2, 'graph_consensus_edge_list.csv'), 'TextType', 'string');
verifyEqual(testCase, E1.source_node_id, E2.source_node_id);
verifyEqual(testCase, E1.target_node_id, E2.target_node_id);
verifyEqual(testCase, E1.consensus_edge_weight, E2.consensus_edge_weight, 'AbsTol', 1e-12);
verifyEqual(testCase, height(unique(E1(:, {'source_node_id','target_node_id'}))), height(E1));
verifyTrue(testCase, all(isfinite(E1.consensus_edge_weight)));
verifyTrue(testCase, all(E1.source_node_id < E1.target_node_id));

S = readtable(fullfile(root1, 'graph_resampled_edge_support_audit.csv'), 'TextType', 'string');
verifyTrue(testCase, all(isfinite(S.conditional_neighbor_support_k5)));
verifyTrue(testCase, all(S.co_inclusion_replicates >= 1));
verifyTrue(testCase, all(S.labels_used_for_edge_support == "none"));
TS = readtable(fullfile(root1, 'graph_consensus_threshold_sensitivity.csv'), 'TextType', 'string');
verifyFalse(testCase, any(logical(TS.session_or_scale_provenance_used_for_threshold_selection)));
H = readtable(fullfile(root1, 'run08_to_run09_handoff_manifest.csv'), 'TextType', 'string');
requiredHandoff = logical(H.required_for_run09);
verifyTrue(testCase, all(strlength(H.sha256(requiredHandoff)) == 64));
verifyTrue(testCase, all(H.artifact_status(requiredHandoff) == ...
    "present_hashed_byte_exact"));
verifyTrue(testCase, all(H.hash_algorithm == "SHA-256"));
verifyTrue(testCase, all(H.hash_scope == "exact_file_bytes"));
hashedRows = find(H.artifact_status == "present_hashed_byte_exact");
for i = hashedRows'
    verifyEqual(testCase, H.sha256(i), compute_file_sha256(H.artifact_path(i)));
end
verifyFalse(testCase, any(H.run09_may_modify_edges));
verifyFalse(testCase, any(H.motifs_defined_in_run08));
verifyEqual(testCase, A1.nConsensusEdges, A2.nConsensusEdges);
R9N = readtable(fullfile(root1, 'run08_to_run09_node_input.csv'), 'TextType', 'string');
R9E = readtable(fullfile(root1, 'run08_to_run09_edge_list.csv'), 'TextType', 'string');
verifyGreaterThanOrEqual(testCase, min(R9N.consensus_stable_induced_degree), ...
    params.consensus_stable_min_degree);
verifyTrue(testCase, all(ismember(R9E.source_node_id, R9N.graph_node_id)));
verifyTrue(testCase, all(ismember(R9E.target_node_id, R9N.graph_node_id)));
for forbiddenColumn = ["anchor_stage","session_index","rare_stratum_id", ...
        "condition_id","arena_label"]
    verifyFalse(testCase, ismember(forbiddenColumn, string(R9N.Properties.VariableNames)));
    verifyFalse(testCase, ismember(forbiddenColumn, string(R9E.Properties.VariableNames)));
end
Sel = readtable(fullfile(root1, 'graph_stage_balanced_resample_manifest.csv'), ...
    'TextType', 'string');
[g, ~, selectedStage] = findgroups(Sel.graph_node_id, Sel.anchor_stage);
exposure = splitapply(@numel, Sel.graph_node_id, g);
verifyTrue(testCase, all(exposure(selectedStage == "base_time_even") == 12));
verifyTrue(testCase, all(exposure(selectedStage == "rare_strata_enriched") == 4));
rareRows = Sel.anchor_stage == "rare_strata_enriched";
[gd, selectedNode, ~] = findgroups(Sel.graph_node_id(rareRows), Sel.dimension_design(rareRows));
designExposure = splitapply(@numel, Sel.graph_node_id(rareRows), gd);
verifyTrue(testCase, all(designExposure == 1));
[gn, ~] = findgroups(selectedNode);
designCount = splitapply(@numel, selectedNode, gn);
verifyTrue(testCase, all(designCount == 4));
figureParams = params;
figureParams.figure_export_png = true;
F = make_run08_consensus_qc_figures(root1, figureParams);
verifyEqual(testCase, height(F), 4);
verifyTrue(testCase, all(isfile(F.figure_file)));
end

function testSha256UsesExactFileBytes(testCase)
root = string(tempname);
mkdir(root);
cleanup = onCleanup(@() i_remove_dirs(root)); %#ok<NASGU>
pathA = fullfile(root, 'same_size_a.bin');
pathB = fullfile(root, 'same_size_b.bin');
i_write_test_bytes(pathA, uint8('abc'));
i_write_test_bytes(pathB, uint8('abd'));

[hashA, bytesA] = compute_file_sha256(pathA);
[hashB, bytesB] = compute_file_sha256(pathB);
verifyEqual(testCase, bytesA, 3);
verifyEqual(testCase, bytesB, 3);
verifyEqual(testCase, hashA, ...
    "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
verifyNotEqual(testCase, hashA, hashB);
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
    "build_run08_embedding_visualization_audits.m", ...
    "build_condition_blind_consensus_neighborhood.m"];
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
    string(fullfile(repoRoot, 'graph', 'build_condition_blind_consensus_neighborhood.m'))
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

verifyTrue(testCase, contains(txt, "rare_stratum_labels_used_for_edge_support"));
verifyTrue(testCase, contains(txt, "event_channels_used_for_consensus"));
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


function params = i_consensus_test_params()
params = struct();
repoRoot = fileparts(fileparts(mfilename('fullpath')));
params.config_path = string(fullfile(repoRoot, 'config', 'multiscale_graph_config.csv'));
params.config_table = table((1:height(readtable(params.config_path)))', ...
    'VariableNames', {'config_row'});
params.graph_n_pcs = 24;
params.k_neighbors = 5;
params.knn_block_size = 32;
params.audit_random_seed = 108;
params.consensus_replicates = 12;
params.consensus_candidate_k = 10;
params.consensus_k_values = [5 10];
params.consensus_support_thresholds = [0.25 0.5];
params.consensus_core_n_pcs = 4;
params.consensus_dimension_fraction = 0.75;
params.consensus_prefix_short_n_pcs = 8;
params.consensus_prefix_medium_n_pcs = 16;
params.consensus_min_co_inclusion_replicates = 1;
params.consensus_stable_min_degree = 2;
params.consensus_rare_min_degree = 1;
params.consensus_gate_min_stable_node_fraction = 0.1;
params.consensus_gate_min_largest_component_fraction = 0.1;
params.consensus_gate_min_rare_stable_fraction = 0.1;
params.consensus_gate_max_same_session_null_ratio = 100;
params.consensus_gate_max_same_scale_fraction_delta = 1;
params.consensus_min_passing_thresholds = 1;
params.consensus_write_edge_support = true;
params.figure_export_png = false;
params.figure_export_pdf = false;
params.figure_dpi = 100;
params.figure_font_name = "Arial";
params.figure_font_size = 8;
end

function i_remove_dirs(paths)
for pathText = paths
    if isfolder(pathText)
        rmdir(pathText, 's');
    end
end
end

function i_write_test_bytes(pathText, bytes)
fid = fopen(pathText, 'wb');
assert(fid >= 0, 'test_16_run08_condition_blind_graph:OpenFailed', ...
    'Could not create SHA-256 regression fixture.');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fwrite(fid, bytes, 'uint8');
end
