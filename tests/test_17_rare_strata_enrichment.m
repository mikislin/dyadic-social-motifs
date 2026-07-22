function tests = test_17_rare_strata_enrichment
tests = functiontests(localfunctions);
end

function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));
end

function testRareEnrichedModeRoutesToSeparateRoots(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
envNames = ["RUN06_ANCHOR_MANIFEST_MODE","RUN06_CHUNK_RUN_MODE", ...
    "RUN06_CHUNK_OUTPUT_DIR","RUN07_ANCHOR_MANIFEST_MODE", ...
    "RUN07_EMBEDDING_RUN_MODE","RUN07_EMBEDDING_OUTPUT_DIR", ...
    "RUN07_CHUNK_INPUT_DIR","RUN08_ANCHOR_MANIFEST_MODE", ...
    "RUN08_GRAPH_RUN_MODE","RUN08_GRAPH_OUTPUT_DIR","RUN08_EMBEDDING_INPUT_DIR"];
oldValues = i_capture_env(envNames);
cleanup = onCleanup(@() i_restore_env(envNames, oldValues)); %#ok<NASGU>
i_clear_env(envNames);

setenv('RUN06_ANCHOR_MANIFEST_MODE', 'rare_enriched');
setenv('RUN06_CHUNK_RUN_MODE', 'full');
setenv('RUN07_ANCHOR_MANIFEST_MODE', 'rare_enriched');
setenv('RUN07_EMBEDDING_RUN_MODE', 'full');
setenv('RUN08_ANCHOR_MANIFEST_MODE', 'rare_enriched');
setenv('RUN08_GRAPH_RUN_MODE', 'full');

p6 = load_multiscale_chunking_config(fullfile(repoRoot, 'config', 'multiscale_chunking_config.csv'));
p7 = load_multiscale_embedding_config(fullfile(repoRoot, 'config', 'multiscale_embedding_config.csv'));
p8 = load_multiscale_graph_config(fullfile(repoRoot, 'config', 'multiscale_graph_config.csv'));

verifyEqual(testCase, p6.output_dir, "derived/chunks_motif_discovery_rare_enriched");
verifyEqual(testCase, p7.chunk_input_dir, "derived/chunks_motif_discovery_rare_enriched");
verifyEqual(testCase, p7.output_dir, "derived/embedding_motif_discovery_rare_enriched");
verifyFalse(testCase, p7.allow_reviewed_snapshot_fallback);
verifyFalse(testCase, p7.use_anchor_weights_for_pca);
verifyEqual(testCase, p8.embedding_input_dir, "derived/embedding_motif_discovery_rare_enriched");
verifyEqual(testCase, p8.output_dir, "derived/graph_motif_discovery_rare_enriched");
verifyNotEqual(testCase, p6.output_dir, p6.baseline_chunk_input_dir);
verifyNotEqual(testCase, p8.output_dir, p8.baseline_graph_input_dir);
end

function testRareStrataIgnorePerturbedLabelColumns(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
graphRoot = fullfile(repoRoot, 'derived', 'graph_motif_discovery');
embeddingRoot = fullfile(repoRoot, 'derived', 'embedding_motif_discovery');
chunkRoot = fullfile(repoRoot, 'derived', 'chunks_motif_discovery');
required = [fullfile(graphRoot, "graph_node_manifest.csv"), ...
    fullfile(graphRoot, "graph_degree_audit.csv"), ...
    fullfile(embeddingRoot, "embedding_row_manifest.csv"), ...
    fullfile(embeddingRoot, "embedding_stability_audit.csv"), ...
    fullfile(chunkRoot, "primary_chunk_event_summary_audit.csv")];
if ~all(isfile(required))
    verifyTrue(testCase, contains(string(fileread(fullfile(repoRoot, 'graph', ...
        'define_condition_blind_rare_strata.m'))), "local_input_table"));
    return
end

Node = readtable(required(1), 'TextType', 'string');
Degree = readtable(required(2), 'TextType', 'string');
Rows = readtable(required(3), 'TextType', 'string');
Stability = readtable(required(4), 'TextType', 'string');
Event = readtable(required(5), 'TextType', 'string');
[D1, M1] = define_condition_blind_rare_strata(graphRoot, ...
    'nodeManifest', Node, 'degreeAudit', Degree, 'rowManifest', Rows, ...
    'eventSummary', Event, 'stabilityAudit', Stability, 'writeOutputs', false);

rng(17);
Node.condition_id = "condition_" + string(randperm(height(Node))');
Degree.cohort_id = "cohort_" + string(randperm(height(Degree))');
Rows.arena_label = "arena_" + string(randperm(height(Rows))');
Event.outcome_label = "outcome_" + string(randperm(height(Event))');
Stability.drug_label = "drug_" + string(randperm(height(Stability))');
[D2, M2] = define_condition_blind_rare_strata(graphRoot, ...
    'nodeManifest', Node, 'degreeAudit', Degree, 'rowManifest', Rows, ...
    'eventSummary', Event, 'stabilityAudit', Stability, 'writeOutputs', false);

definitionColumns = {'rare_stratum_id','scale_index','rare_stratum_threshold', ...
    'rare_stratum_secondary_threshold','rare_stratum_tertiary_threshold','n_member_nodes'};
membershipColumns = {'graph_node_id','scale_index','rare_stratum_id','rare_stratum_score'};
verifyEqual(testCase, D1(:, definitionColumns), D2(:, definitionColumns));
verifyEqual(testCase, M1(:, membershipColumns), M2(:, membershipColumns));
verifyTrue(testCase, all(D1.labels_used_for_rare_strata == "none"));
verifyFalse(testCase, any(D1.condition_used_for_rare_strata));
end

function testRareBuilderDoesNotReadProvenanceLabels(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
files = [fullfile(repoRoot, "graph", "define_condition_blind_rare_strata.m"), ...
    fullfile(repoRoot, "chunks", "build_rare_enriched_primary_anchor_manifest.m")];
txt = "";
for i = 1:numel(files)
    txt = txt + newline + string(fileread(files(i)));
end
forbiddenAccess = [".condition_id", ".condition_label", ".cohort_id", ...
    ".cohort_label", ".arena_label", ".drug", ".genotype", ".outcome"];
for token = forbiddenAccess
    verifyFalse(testCase, contains(lower(txt), token), ...
        "Rare-strata selection code reads a forbidden provenance label: " + token);
end
verifyTrue(testCase, contains(txt, "labels_used_for_anchor_selection"));
verifyTrue(testCase, contains(txt, "labels_used_for_rare_strata"));
end

function testSmokeExpandedBankSchemaAndWeightsIfPresent(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
root = fullfile(repoRoot, 'derived', 'chunks_motif_discovery_rare_enriched_smoke');
manifestPath = fullfile(root, 'expanded_primary_anchor_manifest.csv');
eventPath = fullfile(root, 'expanded_primary_chunk_event_summary_audit.csv');
if ~(isfile(manifestPath) && isfile(eventPath))
    verifyTrue(testCase, contains(string(fileread(fullfile(repoRoot, 'chunks', ...
        'build_rare_enriched_primary_anchor_manifest.m'))), "final_inclusion_probability"));
    return
end

A = readtable(manifestPath, 'TextType', 'string');
E = readtable(eventPath, 'TextType', 'string');
required = ["anchor_stage","rare_stratum_id","rare_stratum_rule", ...
    "rare_stratum_score","base_inclusion_probability", ...
    "rare_inclusion_probability","final_inclusion_probability", ...
    "pca_training_weight","graph_training_weight", ...
    "audit_inverse_probability_weight","selection_phase", ...
    "quota_requested_stratum_id","final_assigned_rare_stratum_id", ...
    "rare_strata_membership_ids","rare_strata_membership_count", ...
    "selection_composite_score","fill_composite_score","fill_reason", ...
    "audit_weight_interpretation"];
verifyTrue(testCase, all(ismember(required, string(A.Properties.VariableNames))));
key = string(A.scale_index) + "_" + string(A.raw_index) + "_" + string(A.anchor_frame);
verifyEqual(testCase, numel(unique(key)), height(A));
for name = ["base_inclusion_probability","rare_inclusion_probability", ...
        "final_inclusion_probability","pca_training_weight", ...
        "graph_training_weight","audit_inverse_probability_weight"]
    verifyTrue(testCase, all(isfinite(double(A.(name)))), name);
end
verifyTrue(testCase, all(A.final_inclusion_probability > 0 & A.final_inclusion_probability <= 1));
verifyEqual(testCase, height(E), height(A));
verifyEqual(testCase, numel(unique(A.scale_index)), 13);
[scaleGroup, ~] = findgroups(A.scale_index);
scaleCounts = splitapply(@numel, A.anchor_frame, scaleGroup);
verifyTrue(testCase, all(scaleCounts == 256));
verifyGreaterThan(testCase, nnz(E.contact_dwell_fraction > 0), 25);
verifyGreaterThan(testCase, nnz(E.contact_transition_count > 0), 22);
verifyTrue(testCase, all(E.labels_used_for_event_summary == "none"));

base = A.anchor_stage == "base_time_even";
rare = A.anchor_stage == "rare_strata_enriched";
quota = rare & A.selection_phase == "quota";
fill = rare & A.selection_phase == "fill";
verifyTrue(testCase, all(A.selection_phase(base) == "retained_base"));
verifyTrue(testCase, all(quota | fill | base));
verifyTrue(testCase, all(A.quota_requested_stratum_id(quota) ~= "none"));
verifyTrue(testCase, all(A.quota_requested_stratum_id(fill) == "none"));
verifyTrue(testCase, all(isfinite(A.fill_composite_score(fill))));
verifyTrue(testCase, all(isnan(A.fill_composite_score(~fill))));
verifyEqual(testCase, A.final_assigned_rare_stratum_id, A.rare_stratum_id);
verifyTrue(testCase, all(A.rare_strata_membership_count(rare) >= 1));
verifyTrue(testCase, all(A.audit_weight_interpretation == ...
    "assigned_stratum_enrichment_sensitivity_not_exact_overlap_adjusted_probability"));

samplingPath = fullfile(root, 'rare_strata_sampling_summary.csv');
verifyTrue(testCase, isfile(samplingPath));
Q = readtable(samplingPath, 'TextType', 'string');
quotaColumns = ["n_selectable_pool_before_base_exclusion", ...
    "n_selectable_after_base_exclusion", ...
    "n_selectable_after_prior_quota_assignments", ...
    "n_excluded_as_locked_base","n_depleted_by_prior_quota_assignments", ...
    "n_selected_quota_stage","n_selected_fill_stage", ...
    "n_selected_exclusive_assignment","n_selected_any_membership", ...
    "n_selected_with_multiple_memberships", ...
    "quota_shortfall_after_quota_stage","quota_shortfall_reason"];
verifyTrue(testCase, all(ismember(quotaColumns, string(Q.Properties.VariableNames))));
verifyEqual(testCase, sum(Q.n_selected_exclusive_assignment), nnz(rare));
verifyEqual(testCase, sum(Q.n_selected_quota_stage) + sum(Q.n_selected_fill_stage), nnz(rare));
verifyEqual(testCase, Q.n_selected_exclusive_assignment, Q.n_selected_rare_anchors);
verifyTrue(testCase, all(Q.n_selected_any_membership >= Q.n_selected_exclusive_assignment));
end

function testSmokeEmbeddingKeepsEventDimensionsSeparateIfPresent(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
root = fullfile(repoRoot, 'derived', 'embedding_motif_discovery_rare_enriched_smoke');
manifestPath = fullfile(root, 'embedding_matrix_manifest.csv');
pcaPath = fullfile(root, 'embedding_pca_by_scale.csv');
if ~(isfile(manifestPath) && isfile(pcaPath))
    verifyTrue(testCase, contains(string(fileread(fullfile(repoRoot, 'embedding', ...
        'run_condition_blind_embedding_build.m'))), "n_event_summary_dims"));
    return
end
M = readtable(manifestPath, 'TextType', 'string');
M = M(M.matrix_role == "per_scale_pca_input", :);
P = readtable(pcaPath, 'TextType', 'string');
verifyTrue(testCase, all(M.n_event_summary_dims == 21));
verifyTrue(testCase, all(M.n_columns == M.run06_expected_multiresolution_dims + M.n_event_summary_dims));
verifyTrue(testCase, all(P.weights_used_for_pca == 0));
verifyTrue(testCase, all(P.n_rare_enriched_rows > 0));
verifyTrue(testCase, all(P.run06_n_temporal_bins(P.hierarchical_role == "motif") == 12));
verifyTrue(testCase, all(P.run06_n_dct_coeffs(P.hierarchical_role == "motif") == 8));
end

function testGraphEdgesIgnoreRareStratumLabels(testCase)
rng(19);
n = 40;
X = randn(n, 6);
Meta = table();
Meta.embedding_row_id = (1:n)';
Meta.scale_index = repelem((1:4)', 10);
Meta.rare_stratum_id = repmat(["contact_present";"graph_periphery"], n / 2, 1);
Meta.rare_stratum_score = rand(n, 1);
params = struct('k_neighbors', 7, 'knn_block_size', 16);
G1 = build_condition_blind_knn_graph(X, Meta, params);
Meta.rare_stratum_id = flipud(Meta.rare_stratum_id);
Meta.rare_stratum_score = randn(n, 1) .* 100;
G2 = build_condition_blind_knn_graph(X, Meta, params);
verifyEqual(testCase, G1.Edges.source_node_id, G2.Edges.source_node_id);
verifyEqual(testCase, G1.Edges.target_node_id, G2.Edges.target_node_id);
verifyEqual(testCase, G1.Edges.neighbor_rank, G2.Edges.neighbor_rank);
verifyEqual(testCase, G1.Edges.neighbor_distance, G2.Edges.neighbor_distance, 'AbsTol', 1e-12);
end

function testRun08ConstructsGraphBeforeRareAudits(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
txt = string(fileread(fullfile(repoRoot, 'graph', 'run_condition_blind_motif_graph_build.m')));
graphPos = strfind(txt, "GraphMax = build_condition_blind_knn_graph(X, nodeManifest, buildParams)");
rarePos = strfind(txt, "define_condition_blind_rare_strata(outRoot");
verifyNotEmpty(testCase, graphPos);
verifyNotEmpty(testCase, rarePos);
verifyLessThan(testCase, graphPos(1), rarePos(1));
verifyTrue(testCase, contains(txt, "local_graph_score_matrix(scoreTable, params)"));
end

function values = i_capture_env(names)
values = strings(size(names));
for i = 1:numel(names)
    values(i) = string(getenv(names(i)));
end
end

function i_restore_env(names, values)
for i = 1:numel(names)
    setenv(char(names(i)), char(values(i)));
end
end

function i_clear_env(names)
for i = 1:numel(names)
    setenv(char(names(i)), '');
end
end
