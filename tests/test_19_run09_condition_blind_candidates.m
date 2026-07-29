function tests = test_19_run09_condition_blind_candidates
tests = functiontests(localfunctions);
end

function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));
end

function testConfigAndRegistryContract(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
envNames = ["RUN09_CANDIDATE_RUN_MODE", "RUN09_CANDIDATE_OUTPUT_DIR"];
oldValues = i_capture_env(envNames);
cleanup = onCleanup(@() i_restore_env(envNames, oldValues)); %#ok<NASGU>
setenv('RUN09_CANDIDATE_RUN_MODE', '');
setenv('RUN09_CANDIDATE_OUTPUT_DIR', '');
params = load_motif_candidate_discovery_config();

verifyEqual(testCase, params.run_mode, "smoke");
verifyEqual(testCase, params.output_dir, ...
    "derived/motif_candidate_discovery_smoke");
verifyEqual(testCase, params.input_dir, ...
    "derived/graph_motif_discovery_rare_enriched");
verifyEqual(testCase, params.expected_freeze_id, ...
    "run08_consensus_d48a1b618ae6_da79abad");
verifyEqual(testCase, params.expected_manifest_schema, ...
    "run08_handoff_manifest_v2_byte_exact");
verifyEqual(testCase, numel(params.cpm_seeds_full), 20);
verifyEqual(testCase, numel(params.active_cpm_density_multipliers), 11);
verifyEqual(testCase, params.primary_objective, "cpm");
verifyEqual(testCase, params.iteration_mode, ...
    "until_membership_unchanged");
verifyEqual(testCase, params.checkpoint_schema, ...
    "run09_replicate_checkpoint_v1");
verifyEqual(testCase, params.retention_min_pairwise_ari, 0.8);
verifyEqual(testCase, params.hierarchy_min_child_parent_fraction, 0.8);
verifyFalse(testCase, params.active_enable_consensus_partition);

setenv('RUN09_CANDIDATE_RUN_MODE', 'full');
setenv('RUN09_CANDIDATE_OUTPUT_DIR', '');
fullParams = load_motif_candidate_discovery_config();
verifyEqual(testCase, fullParams.output_dir, ...
    "derived/motif_candidate_discovery");
verifyEqual(testCase, numel(fullParams.active_cpm_seeds), 20);
verifyTrue(testCase, fullParams.active_enable_consensus_partition);
verifyEqual(testCase, numel(fullParams.active_cpm_density_multipliers), 11);
verifyFalse(testCase, ...
    fullParams.require_all_full_resolutions_pass_retention);
verifyTrue(testCase, ...
    fullParams.require_at_least_one_full_resolution_retained);
verifyEqual(testCase, ...
    fullParams.interpretation_eligibility_rule_version, ...
    "run09_graph_interpretation_v1");
verifyTrue(testCase, ...
    fullParams.interpretation_retain_all_candidate_ids);
verifyEqual(testCase, fullParams.postfreeze_merge_policy, ...
    "prohibited_preserve_frozen_membership");
verifyEqual(testCase, fullParams.postfreeze_split_policy, ...
    "prohibited_preserve_frozen_membership");
verifyEqual(testCase, fullParams.locked_candidate_freeze_id, ...
    "run09_candidates_2e58d214683b_8817d942_0e956a2c");
verifyEqual(testCase, fullParams.locked_candidate_membership_sha256, ...
    "2e58d214683badb0b86e746bf32a7944da664b67243e778ea2d2638249268442");
verifyEqual(testCase, fullParams.postfit_rule_version, ...
    "run09_postfit_interpretation_v1");
verifyEqual(testCase, fullParams.exemplar_selection_rule, ...
    "stability_plus_within_candidate_internal_weighted_degree_rank");
verifyEqual(testCase, fullParams.annotation_rule_version, ...
    "run09_provisional_social_annotation_v1");
verifyTrue(testCase, fullParams.annotation_eligible_candidates_only);
verifyEqual(testCase, fullParams.paper_figure_raster_format, "png");
verifyEqual(testCase, fullParams.paper_figure_vector_format, "pdf");
verifyEqual(testCase, fullParams.paper_figure_raster_dpi, 300);
verifyEqual(testCase, ...
    fullParams.exemplar_stability_weight + ...
    fullParams.exemplar_centrality_weight, 1, 'AbsTol', 1e-15);

registry = load_paper_step_registry(repoRoot);
row = registry(registry.step_id == "run_09", :);
verifyEqual(testCase, height(row), 1);
verifyEqual(testCase, row.script_path, ...
    "paper/run_09_discover_condition_blind_motif_candidates.m");
end

function testTopologyEligibilityAndObjectiveAuditAreGraphOnly(testCase)
params = load_motif_candidate_discovery_config();
params.run_mode = "full";
params.interpretation_small_candidate_max_nodes = 2;
Adapter = build_run09_frozen_graph_adapter(i_small_validation());
reference = uint32([ones(4, 1); 2 * ones(4, 1)]);
alternate = uint32([ones(3, 1); 2 * ones(5, 1)]);
candidateId = ["MC_M001000_C0001";"MC_M001000_C0002"];
resolutionId = repmat("CPM_M001000", 2, 1);
hierarchyLevel = ones(2, 1);
candidateLocalIndex = (1:2)';
nodeCount = 4 * ones(2, 1);
nodePrevalence = 0.5 * ones(2, 1);
meanNodeStability = 0.95 * ones(2, 1);
minimumNodeStability = 0.9 * ones(2, 1);
boundaryNodeFraction = 0.1 * ones(2, 1);
membershipStatus = repmat( ...
    "locked_graph_partition_pre_annotation", 2, 1);
Definition = table(candidateId, resolutionId, hierarchyLevel, ...
    candidateLocalIndex, nodeCount, nodePrevalence, ...
    meanNodeStability, minimumNodeStability, boundaryNodeFraction, ...
    membershipStatus, 'VariableNames', ...
    {'motif_candidate_id','resolution_id','hierarchy_level', ...
    'candidate_local_index','node_count','node_prevalence', ...
    'mean_node_stability','minimum_node_stability', ...
    'boundary_node_fraction','membership_status'});
Analysis = struct();
Analysis.reference_memberships = reference;
Analysis.resolution_ids = "CPM_M001000";
Analysis.reference_replicate_ids = "fixture_cpm_reference";
Analysis.candidate_definition = Definition;
Analysis.resolution_sensitivity = table("CPM_M001000", 0.1, ...
    'VariableNames', {'resolution_id','resolution'});

replicateId = ["mod_a";"mod_b"];
objective = repmat("modularity", 2, 1);
resolutionIndex = ones(2, 1);
resolution = ones(2, 1);
seed = [11;12];
quality = [0.4;0.3];
ReplicateAudit = table(replicateId, objective, resolutionIndex, ...
    resolution, seed, quality, 'VariableNames', ...
    {'replicate_id','objective','resolution_index','resolution', ...
    'seed','quality'});
Ensemble = struct('memberships', [reference alternate], ...
    'replicate_audit', ReplicateAudit);
comparison = compute_run09_partition_comparison(reference, alternate);
Analysis.partition_stability = table("modularity", 1, ...
    "mod_a", "mod_b", comparison.adjusted_rand_index, ...
    'VariableNames', {'objective','resolution_index', ...
    'replicate_id_a','replicate_id_b','adjusted_rand_index'});

Audits = build_run09_candidate_topology_audits( ...
    Adapter, Analysis, Ensemble, params);
T = Audits.candidate_topology;
verifyEqual(testCase, height(T), 2);
verifyEqual(testCase, T.internal_edge_count, 6 * ones(2, 1));
verifyEqual(testCase, T.cut_edge_count, ones(2, 1));
verifyEqual(testCase, T.internal_edge_density, ones(2, 1));
verifyEqual(testCase, T.induced_component_count, ones(2, 1));
verifyEqual(testCase, T.largest_component_fraction, ones(2, 1));
verifyEqual(testCase, T.weighted_conductance, ...
    repmat(0.01 / 12.01, 2, 1), 'AbsTol', 1e-15);
verifyTrue(testCase, all(T.eligible_for_behavioral_interpretation));
verifyTrue(testCase, all(T.candidate_id_retained));
verifyFalse(testCase, any(T.postfreeze_merge_allowed));
verifyFalse(testCase, any(T.postfreeze_split_allowed));
verifyEqual(testCase, height(Audits.objective_concordance), 3);
verifyEqual(testCase, sum( ...
    Audits.objective_concordance.analysis_role == ...
    "frozen_primary_reference"), 1);
verifyEqual(testCase, sum( ...
    Audits.objective_concordance.analysis_role == ...
    "audit_only_modularity_challenger"), 2);
verifyTrue(testCase, all(isfinite( ...
    Audits.objective_concordance{:, vartype('numeric')}), 'all'));
end

function testPartitionArtifactsAreFiniteAndHierarchyIsAudited(testCase)
params = load_motif_candidate_discovery_config();
params.active_retention_min_successful_replicates = 2;
params.retention_min_pairwise_ari = 0;
params.retention_max_normalized_vi = 1;
params.retention_max_candidate_count_relative_range = 1;
params.retention_max_singleton_node_fraction = 1;
params.retention_max_small_candidate_node_fraction = 1;
params.retention_min_mean_node_stability = 0;
params.hierarchy_min_child_parent_fraction = 0.5;
params.hierarchy_min_pair_weighted_child_purity = 0.5;
params.hierarchy_min_linked_node_fraction = 0.5;
Adapter = struct('node_ids', (1:8)');
Memberships = uint32([ ...
    1 1 1 1
    1 1 1 1
    1 1 1 1
    1 1 2 2
    2 2 2 2
    2 2 2 2
    2 2 3 3
    2 2 3 3]);
replicateId = ["a";"b";"c";"d"];
objective = repmat("cpm", 4, 1);
resolutionIndex = [1;1;2;2];
densityMultiplier = [1;1;4;4];
resolution = [0.1;0.1;0.4;0.4];
seed = [1;2;1;2];
quality = [0.5;0.5;0.4;0.4];
candidateCount = [2;2;3;3];
status = repmat("success", 4, 1);
converged = true(4, 1);
membershipHash = repmat("fixture", 4, 1);
ReplicateAudit = table(replicateId, objective, resolutionIndex, ...
    densityMultiplier, resolution, seed, quality, candidateCount, ...
    status, converged, membershipHash, ...
    'VariableNames', {'replicate_id','objective','resolution_index', ...
    'density_multiplier','resolution','seed','quality', ...
    'candidate_count','status','converged','membership_sha256'});
Ensemble = struct('memberships', Memberships, ...
    'replicate_audit', ReplicateAudit);

Analysis = analyze_run09_partition_ensemble(Adapter, Ensemble, params);
verifyEqual(testCase, height(Analysis.partition_stability), 2);
verifyEqual(testCase, height(Analysis.resolution_sensitivity), 2);
verifyEqual(testCase, height(Analysis.node_ambiguity), 16);
verifyEqual(testCase, height(Analysis.node_membership), 16);
verifyEqual(testCase, ...
    groupsummary(Analysis.node_membership, 'graph_node_id').GroupCount, ...
    2 * ones(8, 1));
verifyTrue(testCase, all(isfinite( ...
    Analysis.partition_stability.adjusted_rand_index)));
verifyTrue(testCase, all(isfinite( ...
    Analysis.node_ambiguity.ambiguity_score)));
verifyTrue(testCase, all( ...
    Analysis.hierarchy.child_parent_consistency >= 0 & ...
    Analysis.hierarchy.child_parent_consistency <= 1));
verifyTrue(testCase, all( ...
    Analysis.hierarchy.parent_link_retained <= ...
    Analysis.hierarchy.candidate_link_gate_pass));
end

function testAtomicCheckpointResume(testCase)
assumeNotEmpty(testCase, which('igraph.cluster'), ...
    'matlab-igraph is not installed.');
params = load_motif_candidate_discovery_config();
params.maximum_iterations = 20;
params.active_cpm_density_multipliers = 1;
params.active_cpm_seeds = [101 102];
params.enable_modularity_challenger = false;
outRoot = string(tempname);
mkdir(outRoot);
cleanup = onCleanup(@() i_remove_dir(outRoot)); %#ok<NASGU>
Bridge = prepare_run09_matlab_igraph_bridge(outRoot, params);
cleanupBridge = onCleanup(@() ...
    release_run09_matlab_igraph_bridge(Bridge)); %#ok<NASGU>
Adapter = build_run09_frozen_graph_adapter(i_small_validation());
opts = struct('checkpoint_dir', fullfile(outRoot, 'checkpoints'), ...
    'replicate_audit_path', fullfile(outRoot, 'audit.csv'));

First = execute_run09_partition_ensemble( ...
    Adapter, Bridge, params, opts);
Second = execute_run09_partition_ensemble( ...
    Adapter, Bridge, params, opts);
verifyEqual(testCase, First.memberships, Second.memberships);
verifyEqual(testCase, First.replicate_audit.execution_source, ...
    repmat("computed", 2, 1));
verifyEqual(testCase, Second.replicate_audit.execution_source, ...
    repmat("resumed_checkpoint", 2, 1));
verifyEqual(testCase, numel(dir(fullfile(opts.checkpoint_dir, '*.mat'))), 2);
verifyTrue(testCase, all(strlength( ...
    Second.replicate_audit.checkpoint_sha256) == 64));
verifyTrue(testCase, all(Second.replicate_audit.converged));

release_run09_matlab_igraph_bridge(Bridge);
clear cleanupBridge
end

function testAnnotationMetadataCannotChangeMembership(testCase)
params = load_motif_candidate_discovery_config();
candidateIds = ["MC_FIXTURE_C0001";"MC_FIXTURE_C0002"];
Topology = table(candidateIds, [4;4], [true;false], ...
    ["eligible_for_behavioral_interpretation"; ...
     "residual_small_candidate"], ...
    'VariableNames', {'motif_candidate_id','node_count', ...
    'eligible_for_behavioral_interpretation','interpretation_status'});

parts = strip(string(split(params.annotation_indicator_spec, ';')));
parts = parts(parts ~= "");
behaviorNames = strings(0, 1);
eventNames = strings(0, 1);
for i = 1:numel(parts)
    fields = string(split(parts(i), '|'));
    if fields(1) == "behavior"
        behaviorNames(end + 1, 1) = fields(2); %#ok<AGROW>
    else
        eventNames(end + 1, 1) = fields(2); %#ok<AGROW>
    end
end
Behavior = i_annotation_profile_fixture(candidateIds, ...
    behaviorNames, "behavior");
Event = i_annotation_profile_fixture(candidateIds, ...
    eventNames, "event");
Source = table(strings(0, 1), strings(0, 1), strings(0, 1), ...
    'VariableNames', {'source_id','source_path','sha256'});

membership = uint32([1;1;1;1;2;2;2;2]);
membershipBefore = membership;
membershipHash = compute_run09_membership_sha256(membership);
[A1, P1] = build_run09_provisional_candidate_annotations( ...
    Behavior, Event, Topology, Source, membershipHash, params);
paramsChanged = params;
paramsChanged.annotation_social_score_threshold = 2;
[A2, ~] = build_run09_provisional_candidate_annotations( ...
    Behavior, Event, Topology, Source, membershipHash, paramsChanged);

verifyEqual(testCase, membership, membershipBefore);
verifyEqual(testCase, ...
    compute_run09_membership_sha256(membership), membershipHash);
verifyNotEqual(testCase, A1.provisional_annotation(1), ...
    A2.provisional_annotation(1));
verifyEqual(testCase, A1.provisional_annotation(2), ...
    "not_interpreted_residual_or_unstable");
verifyTrue(testCase, all(A1.annotation_generated_after_candidate_freeze));
verifyFalse(testCase, any(A1.annotation_may_change_membership));
verifyFalse(testCase, ...
    any(A1.annotation_may_change_interpretation_eligibility));
verifyTrue(testCase, all(P1.recorded_after_candidate_freeze));
verifyFalse(testCase, any(P1.used_for_membership));
verifyFalse(testCase, any(P1.used_for_interpretation_eligibility));
end

function P = i_annotation_profile_fixture(candidateIds, featureNames, kind)
nFeature = numel(featureNames);
candidateId = repelem(candidateIds, nFeature);
featureName = repmat(featureNames, numel(candidateIds), 1);
z = [ones(nFeature, 1); -ones(nFeature, 1)];
P = table(candidateId, repmat(kind, numel(candidateId), 1), ...
    featureName, z, ...
    'VariableNames', {'motif_candidate_id','profile_type', ...
    'feature_name','standardized_mean_difference'});
end

function testFullArtifactsFreezeOnlyRetainedResolution(testCase)
params = load_motif_candidate_discovery_config();
params.run_mode = "full";
params.active_retention_min_successful_replicates = 2;
params.retention_min_pairwise_ari = 0;
params.retention_max_normalized_vi = 1;
params.retention_max_candidate_count_relative_range = 1;
params.retention_max_singleton_node_fraction = 1;
params.retention_max_small_candidate_node_fraction = 1;
params.retention_min_mean_node_stability = 0;
Adapter = struct('node_ids', (1:8)');
coarse = uint32([1;1;1;1;2;2;2;2]);
fine = uint32([1;1;1;2;2;2;3;3]);
Memberships = [coarse coarse fine fine];
replicateId = ["r1a";"r1b";"r2a";"r2b"];
objective = repmat("cpm", 4, 1);
resolutionIndex = [1;1;2;2];
densityMultiplier = [1;1;2;2];
resolution = [0.1;0.1;0.2;0.2];
seed = [1;2;1;2];
quality = [0.5;0.5;0.4;0.4];
candidateCount = [2;2;3;3];
status = repmat("success", 4, 1);
converged = true(4, 1);
membershipHash = [ ...
    repmat(compute_run09_membership_sha256(coarse), 2, 1)
    repmat(compute_run09_membership_sha256(fine), 2, 1)];
ReplicateAudit = table(replicateId, objective, resolutionIndex, ...
    densityMultiplier, resolution, seed, quality, candidateCount, ...
    status, converged, membershipHash, ...
    'VariableNames', {'replicate_id','objective','resolution_index', ...
    'density_multiplier','resolution','seed','quality', ...
    'candidate_count','status','converged','membership_sha256'});
Ensemble = struct('memberships', Memberships, ...
    'replicate_audit', ReplicateAudit);

Consensus = struct();
Consensus.memberships = [coarse fine];
Consensus.vote_support = ones(8, 2);
Consensus.split_node_stability_mean = [ones(8, 1), 0.5 * ones(8, 1)];
Consensus.split_node_stability_min = ...
    Consensus.split_node_stability_mean;
Consensus.split_mapped_assignment_rate = ...
    Consensus.split_node_stability_mean;
Consensus.resolution_ids = ["CPM_M001000";"CPM_M002000"];
Consensus.partition_stability = i_consensus_stability_fixture();
hierarchyLevel = repelem((1:2)', 3);
consensusVariant = repmat( ...
    ["full";"odd_seed_split";"even_seed_split"], 2, 1);
consensusStatus = repmat("success", 6, 1);
Consensus.audit = table(hierarchyLevel, consensusVariant, ...
    true(6, 1), consensusStatus, ...
    'VariableNames', {'hierarchy_level','consensus_variant', ...
    'converged','status'});

Analysis = analyze_run09_partition_ensemble( ...
    Adapter, Ensemble, params, Consensus);
verifyEqual(testCase, ...
    Analysis.resolution_sensitivity.retained_by_predeclared_rules, ...
    [true; false]);
verifyEqual(testCase, Analysis.n_retained_resolutions, 1);
verifyEqual(testCase, unique(Analysis.node_membership.resolution_id), ...
    "CPM_M001000");
verifyEqual(testCase, height(Analysis.node_membership), 8);
verifyEqual(testCase, unique(Analysis.node_membership.hierarchy_level), 1);
verifyEqual(testCase, height(Analysis.candidate_definition), 2);
verifyEqual(testCase, height(Analysis.hierarchy), 2);
verifyEqual(testCase, unique(Analysis.hierarchy.link_status), ...
    "retained_root_no_parent");
verifyEqual(testCase, height(Analysis.hierarchy_resolution_audit), 3);
verifyFalse(testCase, Analysis.multilevel_hierarchy_supported);
end

function T = i_consensus_stability_fixture()
objective = repmat("cpm_consensus_split", 6, 1);
resolutionIndex = repelem((1:2)', 3);
densityMultiplier = repelem([1;2], 3);
resolution = repelem([0.1;0.2], 3);
replicateIdA = repmat(["full";"full";"odd"], 2, 1);
replicateIdB = repmat(["odd";"even";"even"], 2, 1);
seedA = ones(6, 1);
seedB = ones(6, 1);
candidateCountA = repelem([2;3], 3);
candidateCountB = candidateCountA;
ari = [ones(3, 1); 0.5 * ones(3, 1)];
vi = [zeros(3, 1); 0.3 * ones(3, 1)];
nvi = vi;
comparisonType = repmat( ...
    ["full_vs_odd";"full_vs_even";"odd_vs_even"], 2, 1);
T = table(objective, resolutionIndex, densityMultiplier, resolution, ...
    replicateIdA, replicateIdB, seedA, seedB, candidateCountA, ...
    candidateCountB, ari, vi, nvi, comparisonType, ...
    'VariableNames', {'objective','resolution_index', ...
    'density_multiplier','resolution','replicate_id_a', ...
    'replicate_id_b','seed_a','seed_b','candidate_count_a', ...
    'candidate_count_b','adjusted_rand_index', ...
    'variation_of_information','normalized_variation_of_information', ...
    'comparison_type'});
end

function testLabelAlignedConsensusUsesFrozenGraphAndResumes(testCase)
assumeNotEmpty(testCase, which('igraph.cluster'), ...
    'matlab-igraph is not installed.');
params = load_motif_candidate_discovery_config();
params.run_mode = "full";
params.maximum_iterations = 20;
params.active_retention_min_successful_replicates = 4;
params.retention_min_pairwise_ari = 0;
params.retention_max_normalized_vi = 1;
params.retention_max_candidate_count_relative_range = 1;
params.retention_max_singleton_node_fraction = 1;
params.retention_max_small_candidate_node_fraction = 1;
params.retention_min_mean_node_stability = 0;
params.consensus_retention_min_split_ari = 0;
params.consensus_retention_max_normalized_vi = 1;
params.consensus_retention_min_mean_node_stability = 0;
params.consensus_retention_min_mean_vote_support = 0;
params.hierarchy_min_pair_weighted_child_purity = 0;
params.hierarchy_min_linked_node_fraction = 0;
outRoot = string(tempname);
mkdir(outRoot);
cleanup = onCleanup(@() i_remove_dir(outRoot)); %#ok<NASGU>
Bridge = prepare_run09_matlab_igraph_bridge(outRoot, params);
cleanupBridge = onCleanup(@() ...
    release_run09_matlab_igraph_bridge(Bridge)); %#ok<NASGU>
Adapter = build_run09_frozen_graph_adapter(i_small_validation());
base = uint32([ones(4, 1); 2 * ones(4, 1)]);
Memberships = repmat(base, 1, 4);
replicateId = compose("raw_%d", (1:4)');
objective = repmat("cpm", 4, 1);
resolutionIndex = ones(4, 1);
densityMultiplier = ones(4, 1);
resolution = repmat(0.05, 4, 1);
seed = (101:104)';
quality = zeros(4, 1);
candidateCount = 2 * ones(4, 1);
status = repmat("success", 4, 1);
converged = true(4, 1);
membershipHash = repmat( ...
    compute_run09_membership_sha256(base), 4, 1);
ReplicateAudit = table(replicateId, objective, resolutionIndex, ...
    densityMultiplier, resolution, seed, quality, candidateCount, ...
    status, converged, membershipHash, ...
    'VariableNames', {'replicate_id','objective','resolution_index', ...
    'density_multiplier','resolution','seed','quality', ...
    'candidate_count','status','converged','membership_sha256'});
Ensemble = struct('memberships', Memberships, ...
    'replicate_audit', ReplicateAudit);
opts = struct('checkpoint_dir', fullfile(outRoot, 'consensus'));

First = build_run09_label_aligned_consensus( ...
    Adapter, Bridge, Ensemble, params, opts);
Second = build_run09_label_aligned_consensus( ...
    Adapter, Bridge, Ensemble, params, opts);
verifyEqual(testCase, First.memberships, base);
verifyEqual(testCase, First.memberships, Second.memberships);
verifyEqual(testCase, height(First.audit), 3);
verifyEqual(testCase, height(First.partition_stability), 3);
verifyEqual(testCase, First.audit.execution_source, ...
    repmat("computed", 3, 1));
verifyEqual(testCase, Second.audit.execution_source, ...
    repmat("resumed_checkpoint", 3, 1));
verifyEqual(testCase, numel(dir(fullfile(opts.checkpoint_dir, '*.mat'))), 3);
Analysis = analyze_run09_partition_ensemble( ...
    Adapter, Ensemble, params, Second);
verifyEqual(testCase, Analysis.resolution_sensitivity.partition_source, ...
    "label_aligned_plurality_then_frozen_graph_refinement");
verifyTrue(testCase, ...
    Analysis.resolution_sensitivity.retained_by_predeclared_rules);
verifyEqual(testCase, height(Analysis.node_membership), 8);
verifyEqual(testCase, unique(Analysis.node_membership.hierarchy_level), 1);
verifyEqual(testCase, unique(Analysis.node_membership.membership_status), ...
    "locked_graph_partition_pre_annotation");
verifyTrue(testCase, all(isfinite( ...
    Analysis.node_ambiguity.consensus_vote_support)));
verifyEqual(testCase, height(Analysis.hierarchy), 2);
verifyEqual(testCase, unique(Analysis.hierarchy.link_status), ...
    "retained_root_no_parent");
verifyTrue(testCase, Analysis.has_retained_resolution);
verifyEqual(testCase, Analysis.n_retained_resolutions, 1);
verifyFalse(testCase, Analysis.multilevel_hierarchy_supported);
verifyEqual(testCase, width(Analysis.hierarchy_resolution_audit), 20);

release_run09_matlab_igraph_bridge(Bridge);
clear cleanupBridge
end

function testPaperScriptThinAndFreezePrecedesBackend(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
paperText = string(fileread(fullfile(repoRoot, 'paper', ...
    'run_09_discover_condition_blind_motif_candidates.m')));
runnerText = string(fileread(fullfile(repoRoot, 'clustering', ...
    'run_condition_blind_motif_candidate_discovery.m')));
verifyFalse(testCase, contains(paperText, newline + "function "));
verifyTrue(testCase, contains(paperText, ...
    "run_condition_blind_motif_candidate_discovery"));
freezeAt = strfind(runnerText, "validate_run09_frozen_handoff");
backendAt = strfind(runnerText, "prepare_run09_matlab_igraph_bridge");
verifyNotEmpty(testCase, freezeAt);
verifyNotEmpty(testCase, backendAt);
verifyLessThan(testCase, freezeAt(1), backendAt(1));

allRun09Source = "";
files = dir(fullfile(repoRoot, 'clustering', '*.m'));
for i = 1:numel(files)
    if startsWith(string(files(i).name), "run09") || ...
            contains(string(files(i).name), "run09_") || ...
            contains(string(files(i).name), "motif_candidate_discovery") || ...
            string(files(i).name) == ...
            "run_weighted_leiden_to_convergence.m"
        allRun09Source = allRun09Source + newline + ...
            string(fileread(fullfile(files(i).folder, files(i).name)));
    end
end
verifyFalse(testCase, contains(allRun09Source, "igraph.cluster("));
consensusText = string(fileread(fullfile(repoRoot, 'clustering', ...
    'build_run09_label_aligned_consensus.m')));
verifyFalse(testCase, contains(consensusText, "Adapter.weights ="));
verifyFalse(testCase, contains(consensusText, "Adapter.graph.Edges.Weight ="));
verifyFalse(testCase, contains(consensusText, "graph("));
verifyTrue(testCase, contains(consensusText, ...
    "run_weighted_leiden_to_convergence(Adapter, Bridge"));
end

function testProductionFreezeAndAdapterPass(testCase)
repoRoot = string(fileparts(fileparts(mfilename('fullpath'))));
params = load_motif_candidate_discovery_config();
auditPath = string(tempname) + ".csv";
cleanup = onCleanup(@() i_delete_file(auditPath)); %#ok<NASGU>
Validation = validate_run09_frozen_handoff(repoRoot, params, auditPath);
Adapter = build_run09_frozen_graph_adapter(Validation);

verifyTrue(testCase, Validation.pass);
verifyEqual(testCase, numel(Adapter.node_ids), 50548);
verifyEqual(testCase, numel(Adapter.weights), 760383);
verifyEqual(testCase, min(Adapter.degree), 10);
verifyTrue(testCase, all(Adapter.audit.pass));
verifyEqual(testCase, Adapter.backend_edge_columns, ...
    ["source_index", "target_index", "consensus_edge_weight"]);
end

function testFreezeIdMismatchFailsClosed(testCase)
repoRoot = string(fileparts(fileparts(mfilename('fullpath'))));
params = load_motif_candidate_discovery_config();
params.expected_freeze_id = "run08_consensus_wrong";
auditPath = string(tempname) + ".csv";
cleanup = onCleanup(@() i_delete_file(auditPath)); %#ok<NASGU>
verifyError(testCase, @() validate_run09_frozen_handoff( ...
    repoRoot, params, auditPath), ...
    'validate_run09_frozen_handoff:IdentityMismatch');
T = readtable(auditPath, 'TextType', 'string');
verifyEqual(testCase, T.pass(T.check_id == "manifest_freeze_id"), 0);
verifyEqual(testCase, T.pass(T.check_id == ...
    "structural_validation_executed"), 0);
end

function testModifiedFrozenWeightFailsByteGate(testCase)
repoRoot = string(fileparts(fileparts(mfilename('fullpath'))));
params = load_motif_candidate_discovery_config();
fixtureRoot = string(tempname);
mkdir(fixtureRoot);
cleanup = onCleanup(@() i_remove_dir(fixtureRoot)); %#ok<NASGU>
[params, edgePath] = i_make_identity_fixture(params, fixtureRoot);

E = readtable(edgePath);
E.consensus_edge_weight(1) = 2;
writetable(E, edgePath);
auditPath = fullfile(fixtureRoot, 'validation.csv');
verifyError(testCase, @() validate_run09_frozen_handoff( ...
    repoRoot, params, auditPath), ...
    'validate_run09_frozen_handoff:IdentityMismatch');
T = readtable(auditPath, 'TextType', 'string');
verifyEqual(testCase, T.pass(T.check_id == ...
    "actual_hash_run08_to_run09_edge_list"), 0);
end

function testAdapterIgnoresForbiddenMetadata(testCase)
Validation = i_small_validation();
Validation.nodes.condition_id = ["a";"b";"c";"d";"e";"f";"g";"h"];
Validation.nodes.session_index = (101:108)';
Validation.nodes.umap_1 = randn(8, 1);
Validation.edges.annotation = repmat("postfit", height(Validation.edges), 1);
Adapter1 = build_run09_frozen_graph_adapter(Validation);

Validation.nodes.condition_id = flipud(Validation.nodes.condition_id);
Validation.nodes.session_index = flipud(Validation.nodes.session_index);
Validation.nodes.umap_1 = 100 * randn(8, 1);
Validation.edges.annotation = repmat("changed", height(Validation.edges), 1);
Adapter2 = build_run09_frozen_graph_adapter(Validation);

verifyEqual(testCase, Adapter1.source_index, Adapter2.source_index);
verifyEqual(testCase, Adapter1.target_index, Adapter2.target_index);
verifyEqual(testCase, Adapter1.weights, Adapter2.weights);
verifyEqual(testCase, string(Adapter1.graph.Edges.Properties.VariableNames), ...
    ["EndNodes", "Weight"]);
verifyEmpty(testCase, Adapter1.graph.Nodes.Properties.VariableNames);
end

function testPinnedNativeBridgeIsWeightedAndDeterministic(testCase)
assumeNotEmpty(testCase, which('igraph.cluster'), ...
    'matlab-igraph is not installed.');
params = load_motif_candidate_discovery_config();
params.maximum_iterations = 20;
outRoot = string(tempname);
mkdir(outRoot);
Bridge = prepare_run09_matlab_igraph_bridge(outRoot, params);
cleanupBridge = onCleanup(@() ...
    release_run09_matlab_igraph_bridge(Bridge)); %#ok<NASGU>
Adapter = build_run09_frozen_graph_adapter(i_small_validation());

R1 = run_weighted_leiden_to_convergence( ...
    Adapter, Bridge, "cpm", 0.05, 240724, params);
R2 = run_weighted_leiden_to_convergence( ...
    Adapter, Bridge, "cpm", 0.05, 240724, params);
verifyTrue(testCase, R1.converged);
verifyEqual(testCase, R1.membership, ...
    uint32([ones(4, 1); 2 * ones(4, 1)]));
verifyEqual(testCase, R1.membership, R2.membership);
verifyEqual(testCase, R1.membership_hash, R2.membership_hash);
verifyEqual(testCase, R1.quality, R2.quality, 'AbsTol', 1e-15);

release_run09_matlab_igraph_bridge(Bridge);
clear cleanupBridge
i_remove_dir(outRoot);
end

function Validation = i_small_validation()
s = [1 1 1 2 2 3 5 5 5 6 6 7 4]';
t = [2 3 4 3 4 4 6 7 8 7 8 8 5]';
w = [ones(12, 1); 0.01];
n = 8;
degree = accumarray([s; t], 1, [n 1]);
weightedDegree = accumarray([s; t], [w; w], [n 1]);
Nodes = table((1:n)', ones(n, 1), ones(n, 1), degree, weightedDegree, ...
    'VariableNames', {'graph_node_id', 'scale_index', 'chunk_sec', ...
    'consensus_stable_induced_degree', 'consensus_weighted_degree'});
Edges = table(s, t, w, 'VariableNames', ...
    {'source_node_id', 'target_node_id', 'consensus_edge_weight'});
Validation = struct('pass', true, 'freeze_id', "fixture", ...
    'nodes', Nodes, 'edges', Edges, 'degree', degree, ...
    'weighted_degree', weightedDegree);
end

function [params, edgePath] = i_make_identity_fixture(params, root)
nodePath = fullfile(root, params.node_input_file);
edgePath = fullfile(root, params.edge_input_file);
freezePath = fullfile(root, params.freeze_config_file);
topologyPath = fullfile(root, params.topology_summary_file);
manifestPath = fullfile(root, params.handoff_manifest_file);

writetable(table((1:2)', 'VariableNames', {'graph_node_id'}), nodePath);
writetable(table(1, 2, 1, 'VariableNames', ...
    {'source_node_id','target_node_id','consensus_edge_weight'}), edgePath);
writetable(table("handoff_ready", "true", "fixture", ...
    'VariableNames', {'parameter','effective_value','freeze_role'}), freezePath);
writetable(table(1, 'VariableNames', {'placeholder'}), topologyPath);

ids = [params.edge_input_file; params.node_input_file; ...
    params.freeze_config_file; params.topology_summary_file];
paths = [edgePath; nodePath; freezePath; topologyPath];
hashes = strings(4, 1);
bytes = zeros(4, 1);
rows = [1; 2; 1; 1];
for i = 1:4
    hashes(i) = compute_file_sha256(paths(i));
    info = dir(paths(i));
    bytes(i) = info.bytes;
end
freezeId = "run08_consensus_" + extractBefore(hashes(1), 13) + ...
    "_" + extractBefore(hashes(3), 9);
Manifest = table(ids, paths, hashes, bytes, rows, true(4, 1), ...
    repmat("present_hashed_byte_exact", 4, 1), ...
    repmat(freezeId, 4, 1), true(4, 1), ...
    repmat("ready_frozen_consensus_graph", 4, 1), false(4, 1), ...
    repmat("SHA-256", 4, 1), repmat("exact_file_bytes", 4, 1), ...
    repmat(params.expected_manifest_schema, 4, 1), ...
    'VariableNames', {'artifact_id','artifact_path','sha256','file_bytes', ...
    'n_rows','required_for_run09','artifact_status','consensus_freeze_id', ...
    'handoff_ready','handoff_status','run09_may_modify_edges', ...
    'hash_algorithm','hash_scope','manifest_schema_version'});
writetable(Manifest, manifestPath);

params.input_dir = root;
params.expected_edge_hash = hashes(1);
params.expected_node_hash = hashes(2);
params.expected_freeze_config_hash = hashes(3);
params.expected_topology_hash = hashes(4);
params.expected_freeze_id = freezeId;
end

function oldValues = i_capture_env(names)
oldValues = strings(size(names));
for i = 1:numel(names)
    oldValues(i) = string(getenv(names(i)));
end
end

function i_restore_env(names, values)
for i = 1:numel(names)
    setenv(names(i), values(i));
end
end

function i_delete_file(pathText)
if isfile(pathText)
    delete(pathText);
end
end

function i_remove_dir(pathText)
if isfolder(pathText)
    rmdir(pathText, 's');
end
end
