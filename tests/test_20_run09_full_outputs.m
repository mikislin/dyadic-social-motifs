function tests = test_20_run09_full_outputs
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repoRoot = string(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(repoRoot));
outRoot = fullfile(repoRoot, 'derived', 'motif_candidate_discovery');
manifestPath = fullfile(outRoot, 'run09_output_manifest.csv');
assumeTrue(testCase, isfile(manifestPath), ...
    'Full run_09 outputs are not present in this checkout.');
testCase.TestData.repoRoot = repoRoot;
testCase.TestData.outRoot = outRoot;
end

function testFrozenInputAndEnsembleCompletion(testCase)
outRoot = testCase.TestData.outRoot;
Freeze = i_read(outRoot, 'run09_input_freeze_validation.csv');
verifyTrue(testCase, all(logical(Freeze.pass)));

Raw = i_read(outRoot, 'run09_partition_replicate_audit.csv');
verifyEqual(testCase, height(Raw), 235);
verifyEqual(testCase, sum(Raw.objective == "cpm"), 220);
verifyEqual(testCase, sum(Raw.objective == "modularity"), 15);
verifyTrue(testCase, all(Raw.status == "success"));
verifyTrue(testCase, all(logical(Raw.converged)));
verifyTrue(testCase, all(isfinite(Raw.resolution)));
verifyTrue(testCase, all(isfinite(Raw.quality)));
verifyTrue(testCase, all(isfinite(Raw.runtime_seconds)));
verifyTrue(testCase, all(Raw.candidate_count >= 1));
verifyTrue(testCase, all(strlength(Raw.membership_sha256) == 64));
verifyTrue(testCase, all(strlength(Raw.checkpoint_sha256) == 64));
verifyEqual(testCase, numel(unique(Raw.replicate_id)), 235);

Consensus = i_read(outRoot, 'run09_consensus_partition_audit.csv');
verifyEqual(testCase, height(Consensus), 33);
verifyEqual(testCase, numel(unique(Consensus.resolution_id)), 11);
verifyTrue(testCase, all( ...
    groupsummary(Consensus, 'resolution_id').GroupCount == 3));
verifyEqual(testCase, sort(unique(Consensus.consensus_variant)), ...
    sort(["full"; "odd_seed_split"; "even_seed_split"]));
verifyTrue(testCase, all(Consensus.status == "success"));
verifyTrue(testCase, all(logical(Consensus.converged)));
verifyTrue(testCase, all(isfinite(Consensus.quality)));
verifyTrue(testCase, all(isfinite(Consensus.mean_vote_support)));
verifyTrue(testCase, all(strlength(Consensus.membership_sha256) == 64));
end

function testOnlyPredeclaredStableResolutionWasFrozen(testCase)
outRoot = testCase.TestData.outRoot;
params = i_full_params();
Sensitivity = i_read( ...
    outRoot, 'motif_candidate_resolution_sensitivity.csv');
verifyEqual(testCase, height(Sensitivity), 11);
retained = logical(Sensitivity.retained_by_predeclared_rules);
verifyEqual(testCase, sum(retained), 1);
verifyEqual(testCase, Sensitivity.resolution_id(retained), ...
    "CPM_M000750");
verifyGreaterThanOrEqual(testCase, ...
    Sensitivity.minimum_consensus_split_ari(retained), ...
    params.consensus_retention_min_split_ari);
verifyLessThanOrEqual(testCase, ...
    Sensitivity.maximum_consensus_split_normalized_vi(retained), ...
    params.consensus_retention_max_normalized_vi);
verifyGreaterThanOrEqual(testCase, ...
    Sensitivity.mean_consensus_split_node_stability(retained), ...
    params.consensus_retention_min_mean_node_stability);
verifyGreaterThanOrEqual(testCase, ...
    Sensitivity.mean_consensus_vote_support(retained), ...
    params.consensus_retention_min_mean_vote_support);
verifyTrue(testCase, logical( ...
    Sensitivity.gate_consensus_convergence(retained)));
end

function testFrozenMembershipDefinitionAndAmbiguity(testCase)
outRoot = testCase.TestData.outRoot;
Membership = i_read(outRoot, 'motif_candidate_node_membership.csv');
Definition = i_read(outRoot, 'motif_candidate_definition.csv');
Ambiguity = i_read(outRoot, 'motif_candidate_node_ambiguity.csv');

verifyEqual(testCase, height(Membership), 50548);
verifyEqual(testCase, numel(unique(Membership.graph_node_id)), 50548);
verifyEqual(testCase, unique(Membership.resolution_id), ...
    "CPM_M000750");
verifyEqual(testCase, unique(Membership.hierarchy_level), 1);
verifyEqual(testCase, unique(Membership.membership_status), ...
    "locked_graph_partition_pre_annotation");
verifyEqual(testCase, numel(unique(Membership.motif_candidate_id)), 42);

verifyEqual(testCase, height(Definition), 42);
verifyEqual(testCase, sum(Definition.node_count), 50548);
verifyEqual(testCase, sum(Definition.node_prevalence), 1, ...
    'AbsTol', 1e-12);
verifyTrue(testCase, all(isfinite( ...
    Definition{:, vartype('numeric')}), 'all'));
verifyTrue(testCase, all(Definition.node_count >= 1));
verifyEqual(testCase, sort(Definition.motif_candidate_id), ...
    sort(unique(Membership.motif_candidate_id)));

verifyEqual(testCase, height(Ambiguity), 50548);
verifyEqual(testCase, sort(Ambiguity.graph_node_id), ...
    sort(Membership.graph_node_id));
scoreNames = {'coassignment_jaccard_mean', ...
    'coassignment_jaccard_min', 'mapped_assignment_rate', ...
    'consensus_vote_support', 'ambiguity_score'};
scores = Ambiguity{:, scoreNames};
verifyTrue(testCase, all(isfinite(scores), 'all'));
verifyTrue(testCase, all(scores >= 0 & scores <= 1, 'all'));
verifyEqual(testCase, Ambiguity.motif_candidate_id, ...
    Membership.motif_candidate_id);
end

function testFlatHierarchyWasNotForced(testCase)
outRoot = testCase.TestData.outRoot;
Hierarchy = i_read(outRoot, 'motif_candidate_hierarchy.csv');
Audit = i_read( ...
    outRoot, 'motif_candidate_hierarchy_resolution_audit.csv');

verifyEqual(testCase, height(Hierarchy), 42);
verifyEqual(testCase, unique(Hierarchy.fine_resolution_id), ...
    "CPM_M000750");
verifyEqual(testCase, unique(Hierarchy.link_status), ...
    "retained_root_no_parent");
verifyTrue(testCase, all(ismissing(Hierarchy.coarse_resolution_id)));
verifyTrue(testCase, ...
    all(ismissing(Hierarchy.hierarchy_parent_candidate_id)));
verifyFalse(testCase, any(logical(Hierarchy.parent_link_retained)));
verifyFalse(testCase, any(logical(Hierarchy.pair_hierarchy_supported)));
verifyEqual(testCase, numel(unique( ...
    Hierarchy.child_motif_candidate_id)), 42);

verifyGreaterThan(testCase, height(Audit), height(Hierarchy));
verifyFalse(testCase, any(logical(Audit.pair_hierarchy_supported)));
verifyFalse(testCase, any(logical(Audit.parent_link_retained)));
verifyTrue(testCase, all(isfinite(Audit.child_parent_consistency)));
verifyTrue(testCase, all(Audit.child_parent_consistency >= 0 & ...
    Audit.child_parent_consistency <= 1));
end

function testMembershipMatchesConsensusBinaryExactly(testCase)
outRoot = testCase.TestData.outRoot;
Membership = i_read(outRoot, 'motif_candidate_node_membership.csv');
S = load(fullfile(outRoot, 'run09_consensus_partitions.mat'), ...
    'memberships', 'resolution_ids');
column = find(S.resolution_ids == "CPM_M000750");
verifyEqual(testCase, numel(column), 1);
labels = uint32(S.memberships(:, column));
verifyEqual(testCase, height(Membership), numel(labels));
expectedIds = compose("MC_M000750_C%04d", labels);
verifyEqual(testCase, Membership.motif_candidate_id, expectedIds);
expectedHash = compute_run09_membership_sha256(labels);
verifyEqual(testCase, ...
    unique(Membership.representative_membership_sha256), expectedHash);
end

function testTopologyEligibilityAndObjectiveComparison(testCase)
outRoot = testCase.TestData.outRoot;
params = i_full_params();
Topology = i_read(outRoot, params.candidate_topology_audit_file);
Definition = i_read(outRoot, params.candidate_definition_file);
Objective = i_read(outRoot, params.objective_concordance_audit_file);

verifyEqual(testCase, height(Topology), 42);
verifyEqual(testCase, sort(Topology.motif_candidate_id), ...
    sort(Definition.motif_candidate_id));
verifyEqual(testCase, sum(Topology.node_count), 50548);
verifyEqual(testCase, ...
    sum(Topology.internal_edge_count) + ...
    sum(Topology.cut_edge_count) / 2, 760383);
verifyEqual(testCase, sum(Topology.unweighted_volume), ...
    2 * 760383);
verifyEqual(testCase, sum(Topology.node_prevalence), 1, ...
    'AbsTol', 1e-12);
verifyTrue(testCase, all(isfinite( ...
    Topology{:, vartype('numeric')}), 'all'));
verifyTrue(testCase, all(Topology.weighted_conductance >= 0 & ...
    Topology.weighted_conductance <= 1));
verifyTrue(testCase, all(Topology.internal_edge_density >= 0 & ...
    Topology.internal_edge_density <= 1));
verifyTrue(testCase, all(Topology.largest_component_fraction > 0 & ...
    Topology.largest_component_fraction <= 1));
verifyTrue(testCase, all(Topology.candidate_id_retained));
verifyFalse(testCase, any(Topology.postfreeze_merge_allowed));
verifyFalse(testCase, any(Topology.postfreeze_split_allowed));
expectedEligibility = ~( ...
    Topology.is_residual_small_candidate | ...
    Topology.is_stability_weak_candidate | ...
    Topology.is_boundary_heavy_candidate | ...
    Topology.is_topologically_weak_candidate);
verifyEqual(testCase, ...
    logical(Topology.eligible_for_behavioral_interpretation), ...
    logical(expectedEligibility));
verifyEqual(testCase, ...
    unique(Topology.interpretation_eligibility_rule_version), ...
    params.interpretation_eligibility_rule_version);

verifyEqual(testCase, height(Objective), 16);
reference = Objective.analysis_role == "frozen_primary_reference";
challenger = Objective.analysis_role == ...
    "audit_only_modularity_challenger";
verifyEqual(testCase, sum(reference), 1);
verifyEqual(testCase, sum(challenger), 15);
verifyEqual(testCase, Objective.objective(reference), "cpm");
verifyTrue(testCase, all(Objective.objective(challenger) == ...
    "modularity"));
verifyEqual(testCase, ...
    Objective.adjusted_rand_index_to_retained_cpm(reference), 1);
verifyEqual(testCase, ...
    Objective.variation_of_information_to_retained_cpm(reference), 0);
verifyEqual(testCase, ...
    Objective.normalized_variation_of_information_to_retained_cpm( ...
    reference), 0);
verifyEqual(testCase, sum( ...
    Objective.is_within_resolution_ari_medoid(challenger)), 3);
verifyTrue(testCase, all(isfinite( ...
    Objective{:, vartype('numeric')}), 'all'));
verifyTrue(testCase, all( ...
    Objective.internal_weight_fraction >= 0 & ...
    Objective.internal_weight_fraction <= 1));
verifyTrue(testCase, all( ...
    Objective.node_weighted_mean_conductance >= 0 & ...
    Objective.node_weighted_mean_conductance <= 1));
end

function testPostfitCompositionProfilesExemplarsAndAnnotations(testCase)
outRoot = testCase.TestData.outRoot;
params = i_full_params();
Topology = i_read(outRoot, params.candidate_topology_audit_file);
Definition = i_read(outRoot, params.candidate_definition_file);
Scale = i_read(outRoot, params.scale_composition_audit_file);
Session = i_read(outRoot, params.session_composition_audit_file);
Behavior = i_read(outRoot, params.behavioral_profile_file);
Event = i_read(outRoot, params.event_profile_file);
Exemplar = i_read(outRoot, params.exemplar_manifest_file);
Annotation = i_read(outRoot, params.annotation_file);
Provenance = i_read(outRoot, params.annotation_provenance_file);
Figures = i_read(outRoot, params.qc_figure_manifest_file);
PaperFigures = i_read(outRoot, params.paper_figure_manifest_file);

verifyEqual(testCase, height(Scale), 42 * 13);
verifyEqual(testCase, height(Session), 42 * 197);
verifyEqual(testCase, height(Behavior), 42 * 39);
verifyEqual(testCase, height(Event), 42 * 21);
verifyEqual(testCase, height(Annotation), 42);
verifyEqual(testCase, height(Figures), 6);
verifyEqual(testCase, height(PaperFigures), 6);
verifyEqual(testCase, sort(unique(Scale.motif_candidate_id)), ...
    sort(Definition.motif_candidate_id));
verifyEqual(testCase, sort(unique(Session.motif_candidate_id)), ...
    sort(Definition.motif_candidate_id));

ScaleSum = groupsummary(Scale, 'motif_candidate_id', 'sum', ...
    'node_count');
SessionSum = groupsummary(Session, 'motif_candidate_id', 'sum', ...
    'node_count');
[~, dLocScale] = ismember(ScaleSum.motif_candidate_id, ...
    Definition.motif_candidate_id);
[~, dLocSession] = ismember(SessionSum.motif_candidate_id, ...
    Definition.motif_candidate_id);
verifyEqual(testCase, ScaleSum.sum_node_count, ...
    Definition.node_count(dLocScale));
verifyEqual(testCase, SessionSum.sum_node_count, ...
    Definition.node_count(dLocSession));
verifyEqual(testCase, groupsummary(Scale, ...
    'motif_candidate_id', 'sum', ...
    'candidate_fraction').sum_candidate_fraction, ...
    ones(42, 1), 'AbsTol', 1e-12);
verifyEqual(testCase, groupsummary(Session, ...
    'motif_candidate_id', 'sum', ...
    'candidate_fraction').sum_candidate_fraction, ...
    ones(42, 1), 'AbsTol', 1e-12);
verifyTrue(testCase, all(isfinite( ...
    Scale{:, vartype('numeric')}), 'all'));
verifyTrue(testCase, all(isfinite( ...
    Session{:, vartype('numeric')}), 'all'));
verifyFalse(testCase, any(logical(Scale.used_for_membership)));
verifyFalse(testCase, ...
    any(logical(Scale.used_for_interpretation_eligibility)));
verifyFalse(testCase, any(logical(Session.used_for_membership)));
verifyFalse(testCase, ...
    any(logical(Session.used_for_interpretation_eligibility)));

verifyTrue(testCase, all(isfinite( ...
    Behavior{:, vartype('numeric')}), 'all'));
verifyTrue(testCase, all(isfinite( ...
    Event{:, vartype('numeric')}), 'all'));
verifyEqual(testCase, numel(unique(Behavior.feature_name)), 39);
verifyEqual(testCase, numel(unique(Event.feature_name)), 21);
verifyTrue(testCase, all(Behavior.experimental_labels_used == "none"));
verifyTrue(testCase, all(Event.experimental_labels_used == "none"));
verifyFalse(testCase, any(logical(Behavior.used_for_membership)));
verifyFalse(testCase, ...
    any(logical(Behavior.used_for_interpretation_eligibility)));

verifyEqual(testCase, numel(unique(Exemplar.graph_node_id)), ...
    height(Exemplar));
verifyLessThanOrEqual(testCase, max(Exemplar.exemplar_rank), ...
    params.exemplar_count_per_candidate);
verifyEqual(testCase, Exemplar.exemplar_score, ...
    params.exemplar_stability_weight .* Exemplar.node_stability + ...
    params.exemplar_centrality_weight .* ...
    Exemplar.centrality_percentile, 'AbsTol', 1e-12);
verifyFalse(testCase, any(logical(Exemplar.used_for_membership)));
verifyFalse(testCase, any(logical( ...
    Exemplar.used_for_interpretation_eligibility)));
verifyTrue(testCase, all(contains( ...
    Exemplar.postfit_metadata_not_used_for_selection, ...
    "attached_after_selection")));

[tf, topologyLoc] = ismember(Annotation.motif_candidate_id, ...
    Topology.motif_candidate_id);
verifyTrue(testCase, all(tf));
verifyEqual(testCase, logical( ...
    Annotation.eligible_for_behavioral_interpretation), ...
    logical(Topology.eligible_for_behavioral_interpretation( ...
    topologyLoc)));
withheld = ~logical( ...
    Annotation.eligible_for_behavioral_interpretation);
verifyTrue(testCase, all(Annotation.provisional_annotation(withheld) == ...
    "not_interpreted_residual_or_unstable"));
verifyFalse(testCase, any(logical( ...
    Annotation.annotation_may_change_membership)));
verifyFalse(testCase, any(logical( ...
    Annotation.annotation_may_change_interpretation_eligibility)));
verifyEqual(testCase, unique(Annotation.frozen_membership_sha256), ...
    params.locked_candidate_membership_sha256);
verifyTrue(testCase, all(Provenance.experimental_labels_used == "none"));
verifyFalse(testCase, any(logical(Provenance.used_for_membership)));
verifyFalse(testCase, any(logical( ...
    Provenance.used_for_interpretation_eligibility)));
sourceRows = Provenance.record_type == "source_artifact";
for i = find(sourceRows)'
    verifyTrue(testCase, isfile(Provenance.source_path(i)));
    verifyEqual(testCase, compute_file_sha256( ...
        Provenance.source_path(i)), Provenance.sha256(i));
end

for i = 1:height(Figures)
    verifyTrue(testCase, isfile(Figures.figure_path(i)));
    verifyEqual(testCase, ...
        compute_file_sha256(Figures.figure_path(i)), ...
        Figures.sha256(i));
    info = dir(Figures.figure_path(i));
    verifyEqual(testCase, double(info.bytes), Figures.file_bytes(i));
end
verifyTrue(testCase, all(logical( ...
    Figures.generated_after_candidate_freeze)));
verifyFalse(testCase, any(logical(Figures.used_for_membership)));
verifyFalse(testCase, any(logical( ...
    Figures.used_for_interpretation_eligibility)));
verifyFalse(testCase, any(logical(Figures.used_for_annotation)));
verifyEqual(testCase, ...
    numel(unique(PaperFigures.recommended_paper_role)), 3);
for i = 1:height(PaperFigures)
    verifyTrue(testCase, isfile(PaperFigures.figure_path(i)));
    verifyEqual(testCase, ...
        compute_file_sha256(PaperFigures.figure_path(i)), ...
        PaperFigures.sha256(i));
    info = dir(PaperFigures.figure_path(i));
    verifyEqual(testCase, double(info.bytes), ...
        PaperFigures.file_bytes(i));
end
verifyTrue(testCase, all(logical( ...
    PaperFigures.generated_after_candidate_freeze)));
verifyFalse(testCase, any(logical(PaperFigures.used_for_membership)));
verifyFalse(testCase, any(logical( ...
    PaperFigures.used_for_interpretation_eligibility)));
verifyFalse(testCase, any(logical(PaperFigures.used_for_annotation)));
end

function testManifestRehashesEveryRecordedFile(testCase)
repoRoot = testCase.TestData.repoRoot;
outRoot = testCase.TestData.outRoot;
params = i_full_params();
Manifest = i_read(outRoot, 'run09_output_manifest.csv');
fileRows = Manifest.record_type == "output_artifact" | ...
    Manifest.record_type == "frozen_input" | ...
    Manifest.record_type == "source_provenance" | ...
    Manifest.record_type == "qc_figure" | ...
    Manifest.record_type == "paper_figure";
Files = Manifest(fileRows, :);
verifyTrue(testCase, all(logical(Manifest.no_experimental_labels_used)));
verifyEqual(testCase, sum(logical( ...
    Manifest.post_fit_interpretation)), 21);
verifyEqual(testCase, numel(unique(Manifest.run08_freeze_id)), 1);
verifyEqual(testCase, unique(Manifest.run08_freeze_id), ...
    params.expected_freeze_id);
verifyEqual(testCase, numel(unique(Manifest.candidate_freeze_id)), 1);
verifyEqual(testCase, unique(Manifest.candidate_freeze_id), ...
    params.locked_candidate_freeze_id);
verifyEqual(testCase, numel(unique(Manifest.config_sha256)), 1);
verifyEqual(testCase, unique(Manifest.config_sha256), ...
    compute_file_sha256(params.config_path));

for i = 1:height(Files)
    pathText = Files.artifact_path(i);
    verifyTrue(testCase, isfile(pathText), ...
        "Missing manifest file: " + pathText);
    verifyEqual(testCase, compute_file_sha256(pathText), ...
        Files.sha256(i), "Hash mismatch: " + pathText);
    info = dir(pathText);
    verifyEqual(testCase, double(info.bytes), Files.file_bytes(i), ...
        "Byte-count mismatch: " + pathText);
    if endsWith(lower(pathText), ".csv") && ...
            isfinite(Files.row_count(i))
        verifyEqual(testCase, height(readtable(pathText, ...
            'Delimiter', ',')), ...
            Files.row_count(i), "Row-count mismatch: " + pathText);
    end
end

Output = Manifest(Manifest.record_type == "output_artifact", :);
verifyEqual(testCase, sum(logical(Output.defines_membership)), 3);
verifyEqual(testCase, sort(Output.artifact_id( ...
    logical(Output.defines_membership))), sort([ ...
    params.node_membership_file
    params.candidate_definition_file
    params.hierarchy_file]));
verifyEqual(testCase, sum(Output.artifact_id == ...
    params.candidate_topology_audit_file), 1);
verifyEqual(testCase, sum(Output.artifact_id == ...
    params.objective_concordance_audit_file), 1);
postfitIds = [ ...
    params.scale_composition_audit_file
    params.session_composition_audit_file
    params.behavioral_profile_file
    params.event_profile_file
    params.exemplar_manifest_file
    params.annotation_file
    params.annotation_provenance_file
    params.qc_figure_manifest_file
    params.paper_figure_manifest_file];
verifyEqual(testCase, sum(logical(Output.post_fit_interpretation)), 9);
verifyEqual(testCase, sort(Output.artifact_id( ...
    logical(Output.post_fit_interpretation))), sort(postfitIds));
verifyFalse(testCase, logical(Output.defines_membership( ...
    Output.artifact_id == params.candidate_topology_audit_file)));
verifyFalse(testCase, logical(Output.defines_membership( ...
    Output.artifact_id == params.objective_concordance_audit_file)));

Input = Manifest(Manifest.record_type == "frozen_input", :);
expectedInputIds = [params.node_input_file; params.edge_input_file; ...
    params.freeze_config_file; params.topology_summary_file];
expectedInputHashes = [params.expected_node_hash; ...
    params.expected_edge_hash; params.expected_freeze_config_hash; ...
    params.expected_topology_hash];
for i = 1:numel(expectedInputIds)
    row = Input.artifact_id == expectedInputIds(i);
    verifyEqual(testCase, sum(row), 1);
    verifyEqual(testCase, Input.sha256(row), expectedInputHashes(i));
end

contract = @(id) Manifest.details( ...
    Manifest.record_type == "run_contract" & ...
    Manifest.artifact_id == id);
verifyEqual(testCase, contract("forbidden_labels_used"), "none");
verifyEqual(testCase, contract("retained_resolution_count"), "1");
verifyEqual(testCase, contract("retained_resolution_ids"), ...
    "CPM_M000750");
verifyEqual(testCase, contract("candidate_freeze_status"), ...
    "frozen_graph_partition_pre_annotation");
verifyEqual(testCase, contract("postfreeze_membership_policy"), ...
    "retain_all_ids=true;merge=prohibited;split=prohibited");
verifyEqual(testCase, contract("objective_challenger_role"), ...
    "weighted_modularity_audit_only_cannot_select_or_replace_membership");
verifyEqual(testCase, contract("scale_session_audit_role"), ...
    "post_fit_only_cannot_change_membership_or_interpretation_eligibility");
verifyEqual(testCase, contract("annotation_membership_invariance"), ...
    "annotation_metadata_and_labels_are_downstream_outputs;membership_hash_remains_pinned");
verifyEqual(testCase, contract("candidate_membership_freeze_hash"), ...
    params.locked_candidate_membership_sha256);
verifyEqual(testCase, contract("partition_input_columns"), ...
    "source_index;target_index;consensus_edge_weight");

membershipHash = compute_file_sha256(fullfile( ...
    repoRoot, 'derived', 'motif_candidate_discovery', ...
    params.node_membership_file));
verifyEqual(testCase, membershipHash, ...
    params.locked_candidate_membership_sha256);
end

function params = i_full_params()
oldMode = string(getenv('RUN09_CANDIDATE_RUN_MODE'));
oldOut = string(getenv('RUN09_CANDIDATE_OUTPUT_DIR'));
cleanup = onCleanup(@() i_restore_env(oldMode, oldOut));
setenv('RUN09_CANDIDATE_RUN_MODE', 'full');
setenv('RUN09_CANDIDATE_OUTPUT_DIR', '');
params = load_motif_candidate_discovery_config();
end

function T = i_read(outRoot, fileName)
T = readtable(fullfile(outRoot, fileName), ...
    'Delimiter', ',', 'TextType', 'string');
end

function i_restore_env(modeValue, outValue)
setenv('RUN09_CANDIDATE_RUN_MODE', modeValue);
setenv('RUN09_CANDIDATE_OUTPUT_DIR', outValue);
end
