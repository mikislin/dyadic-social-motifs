function [Manifest, candidateFreezeId] = build_run09_output_manifest( ...
        repoRoot, outRoot, params, Validation)
%BUILD_RUN09_OUTPUT_MANIFEST Hash run_09 inputs sources contracts and outputs.

repoRoot = string(repoRoot);
outRoot = string(outRoot);
configHash = compute_file_sha256(params.config_path);
sourcePaths = i_source_paths(repoRoot);
sourceHashes = strings(numel(sourcePaths), 1);
for i = 1:numel(sourcePaths)
    sourceHashes(i) = compute_file_sha256(sourcePaths(i));
end
[sourcePaths, order] = sort(sourcePaths);
sourceHashes = sourceHashes(order);
sourcePayload = strjoin(sourcePaths + "=" + sourceHashes, newline);
sourceBundleHash = compute_run09_sha256_text(sourcePayload);

outputNames = i_output_names(params);
outputRoles = i_output_roles(params);
postFitFlags = ismember(outputNames, [ ...
    params.scale_composition_audit_file
    params.session_composition_audit_file
    params.behavioral_profile_file
    params.event_profile_file
    params.exemplar_manifest_file
    params.annotation_file
    params.annotation_provenance_file
    params.qc_figure_manifest_file
    params.paper_figure_manifest_file]);
definesMembership = ismember(outputNames, [ ...
    params.node_membership_file
    params.candidate_definition_file
    params.hierarchy_file]);
membershipPath = fullfile(outRoot, params.node_membership_file);
membershipHash = compute_file_sha256(membershipPath);
if params.run_mode == "full"
    assert(membershipHash == params.locked_candidate_membership_sha256, ...
        'build_run09_output_manifest:CandidateMembershipFreezeMismatch', ...
        ['The frozen candidate membership bytes changed after the issued ' ...
         'candidate freeze. Audit-only work may not rename or rewrite it.']);
    candidateFreezeId = params.locked_candidate_freeze_id;
else
    candidateFreezeId = "not_frozen_smoke";
    definesMembership(:) = false;
end

Manifest = i_empty_manifest();
for i = 1:numel(outputNames)
    pathText = fullfile(outRoot, outputNames(i));
    assert(isfile(pathText), ...
        'build_run09_output_manifest:MissingOutput', ...
        'Required generated run_09 output is missing: %s', pathText);
    [hash, bytes] = compute_file_sha256(pathText);
    rows = i_csv_rows(pathText);
    Manifest = [Manifest; i_row("output_artifact", outputNames(i), ... %#ok<AGROW>
        pathText, outputRoles(i), definesMembership(i), postFitFlags(i), ...
        "present_hashed_byte_exact", rows, bytes, hash, ...
        Validation.freeze_id, candidateFreezeId, configHash, ...
        sourceBundleHash, true, ...
        i_output_details(definesMembership(i), postFitFlags(i)))];
end

Manifest = i_append_figure_manifest(Manifest, outRoot, ...
    params.qc_figure_manifest_file, "qc_figure", Validation.freeze_id, ...
    candidateFreezeId, configHash, sourceBundleHash);
Manifest = i_append_figure_manifest(Manifest, outRoot, ...
    params.paper_figure_manifest_file, "paper_figure", ...
    Validation.freeze_id, candidateFreezeId, configHash, ...
    sourceBundleHash);

inputNames = [params.node_input_file; params.edge_input_file; ...
    params.freeze_config_file; params.topology_summary_file];
inputRoles = [ ...
    "frozen_node_identity_and_provenance"
    "frozen_weighted_graph_edges"
    "frozen_graph_contract"
    "frozen_topology_readiness_audit"];
for i = 1:numel(inputNames)
    pathText = fullfile(resolve_repo_path(repoRoot, params.input_dir), ...
        inputNames(i));
    [hash, bytes] = compute_file_sha256(pathText);
    rows = i_csv_rows(pathText);
    Manifest = [Manifest; i_row("frozen_input", inputNames(i), ... %#ok<AGROW>
        pathText, inputRoles(i), false, false, ...
        "validated_read_only_byte_exact", rows, bytes, hash, ...
        Validation.freeze_id, candidateFreezeId, configHash, ...
        sourceBundleHash, true, ...
        "Run_09 did not modify rebuild prune supplement or reweight this input.")];
end

for i = 1:numel(sourcePaths)
    [~, stem, ext] = fileparts(sourcePaths(i));
    Manifest = [Manifest; i_row("source_provenance", ... %#ok<AGROW>
        string(stem) + string(ext), sourcePaths(i), ...
        "run09_implementation_source", false, false, ...
        "present_hashed_byte_exact", NaN, ...
        i_file_bytes(sourcePaths(i)), sourceHashes(i), ...
        Validation.freeze_id, candidateFreezeId, configHash, ...
        sourceBundleHash, true, ...
        "Exact source contributing to the locked run_09 implementation.")];
end

Sensitivity = readtable(fullfile(outRoot, ...
    params.resolution_sensitivity_file), 'TextType', 'string');
retainedIds = Sensitivity.resolution_id( ...
    logical(Sensitivity.retained_by_predeclared_rules));
contractIds = [ ...
    "algorithm_and_backend"
    "primary_resolution_grid"
    "primary_random_seeds"
    "challenger_resolution_grid"
    "challenger_random_seeds"
    "resolution_retention_rules"
    "hierarchy_rules"
    "partition_input_columns"
    "forbidden_labels_used"
    "arena_session_usage"
    "annotation_dependency"
    "interpretation_eligibility_rules"
    "postfreeze_membership_policy"
    "objective_challenger_role"
    "scale_session_audit_role"
    "behavior_event_profile_role"
    "exemplar_selection_rule"
    "provisional_annotation_rule"
    "annotation_membership_invariance"
    "candidate_membership_freeze_hash"
    "retained_resolution_count"
    "retained_resolution_ids"
    "candidate_freeze_status"];
contractValues = [ ...
    params.backend_id + "|weighted_leiden_cpm|" + ...
        params.consensus_partition_rule
    strjoin(compose('%.17g', params.active_cpm_density_multipliers), ";")
    strjoin(compose('%.0f', params.active_cpm_seeds), ";")
    strjoin(compose('%.17g', params.active_modularity_resolutions), ";")
    strjoin(compose('%.0f', params.active_modularity_seeds), ";")
    "split_ari>=" + compose('%.17g', ...
        params.consensus_retention_min_split_ari) + ...
        ";normalized_vi<=" + compose('%.17g', ...
        params.consensus_retention_max_normalized_vi) + ...
        ";mean_node_stability>=" + compose('%.17g', ...
        params.consensus_retention_min_mean_node_stability) + ...
        ";mean_vote_support>=" + compose('%.17g', ...
        params.consensus_retention_min_mean_vote_support)
    "child_parent_fraction>=" + compose('%.17g', ...
        params.hierarchy_min_child_parent_fraction) + ...
        ";pair_weighted_purity>=" + compose('%.17g', ...
        params.hierarchy_min_pair_weighted_child_purity) + ...
        ";linked_node_fraction>=" + compose('%.17g', ...
        params.hierarchy_min_linked_node_fraction)
    "source_index;target_index;consensus_edge_weight"
    "none"
    "post_fit_audit_only_not_used_for_selection"
    "none_candidate_ids_frozen_before_annotation"
    "small_nodes<=" + compose('%.0f', ...
        params.interpretation_small_candidate_max_nodes) + ...
        ";mean_node_stability>=" + compose('%.17g', ...
        params.interpretation_min_mean_node_stability) + ...
        ";boundary_fraction<" + compose('%.17g', ...
        params.interpretation_max_boundary_node_fraction) + ...
        ";weighted_conductance<=" + compose('%.17g', ...
        params.interpretation_max_weighted_conductance) + ...
        ";largest_component_fraction>=" + compose('%.17g', ...
        params.interpretation_min_largest_component_fraction) + ...
        ";internal_edges>=" + compose('%.0f', ...
        params.interpretation_min_internal_edge_count)
    "retain_all_ids=true;merge=prohibited;split=prohibited"
    "weighted_modularity_audit_only_cannot_select_or_replace_membership"
    "post_fit_only_cannot_change_membership_or_interpretation_eligibility"
    "condition_blind_post_fit_only_cannot_change_membership_or_interpretation_eligibility"
    params.exemplar_selection_rule + ...
        ";stability_weight=" + compose('%.17g', ...
        params.exemplar_stability_weight) + ...
        ";centrality_weight=" + compose('%.17g', ...
        params.exemplar_centrality_weight)
    params.annotation_rule_version + ...
        ";eligible_only=true;social_threshold=" + ...
        compose('%.17g', params.annotation_social_score_threshold) + ...
        ";non_social_threshold=" + ...
        compose('%.17g', params.annotation_non_social_score_threshold)
    "annotation_metadata_and_labels_are_downstream_outputs;membership_hash_remains_pinned"
    params.locked_candidate_membership_sha256
    string(numel(retainedIds))
    strjoin(retainedIds, ";")
    i_freeze_status(params)];
for i = 1:numel(contractIds)
    Manifest = [Manifest; i_row("run_contract", contractIds(i), ... %#ok<AGROW>
        "", "reviewer_audit_contract", false, false, "locked", ...
        NaN, NaN, "", Validation.freeze_id, candidateFreezeId, ...
        configHash, sourceBundleHash, true, contractValues(i))];
end

function value = i_freeze_status(params)
if params.run_mode == "full"
    value = "frozen_graph_partition_pre_annotation";
else
    value = "provisional_smoke_not_frozen";
end
end
end

function paths = i_source_paths(repoRoot)
paths = [
    fullfile(repoRoot, 'config', 'motif_candidate_discovery_config.csv')
    fullfile(repoRoot, 'paper', ...
        'run_09_discover_condition_blind_motif_candidates.m')
    fullfile(repoRoot, 'clustering', ...
        'run_condition_blind_motif_candidate_discovery.m')
    fullfile(repoRoot, 'clustering', ...
        'refresh_run09_postfit_outputs.m')
    fullfile(repoRoot, 'clustering', ...
        'load_motif_candidate_discovery_config.m')
    fullfile(repoRoot, 'clustering', 'validate_run09_frozen_handoff.m')
    fullfile(repoRoot, 'clustering', ...
        'build_run09_frozen_graph_adapter.m')
    fullfile(repoRoot, 'clustering', ...
        'prepare_run09_matlab_igraph_bridge.m')
    fullfile(repoRoot, 'clustering', ...
        'validate_matlab_igraph_backend.m')
    fullfile(repoRoot, 'clustering', ...
        'release_run09_matlab_igraph_bridge.m')
    fullfile(repoRoot, 'clustering', ...
        'run_weighted_leiden_to_convergence.m')
    fullfile(repoRoot, 'clustering', ...
        'execute_run09_partition_ensemble.m')
    fullfile(repoRoot, 'clustering', ...
        'build_run09_label_aligned_consensus.m')
    fullfile(repoRoot, 'clustering', ...
        'analyze_run09_partition_ensemble.m')
    fullfile(repoRoot, 'clustering', ...
        'build_run09_candidate_topology_audits.m')
    fullfile(repoRoot, 'clustering', ...
        'build_run09_postfit_candidate_audits.m')
    fullfile(repoRoot, 'clustering', ...
        'build_run09_provisional_candidate_annotations.m')
    fullfile(repoRoot, 'clustering', ...
        'make_run09_candidate_qc_figures.m')
    fullfile(repoRoot, 'clustering', ...
        'make_run09_candidate_paper_figures.m')
    fullfile(repoRoot, 'clustering', ...
        'canonicalize_run09_membership.m')
    fullfile(repoRoot, 'clustering', ...
        'compute_run09_partition_quality.m')
    fullfile(repoRoot, 'clustering', ...
        'compute_run09_partition_comparison.m')
    fullfile(repoRoot, 'clustering', ...
        'compute_run09_membership_sha256.m')
    fullfile(repoRoot, 'clustering', ...
        'compute_run09_sha256_text.m')
    fullfile(repoRoot, 'clustering', ...
        'save_run09_mat_atomic.m')
    fullfile(repoRoot, 'clustering', ...
        'write_run09_table_atomic.m')
    fullfile(repoRoot, 'io', 'write_table_atomic.m')
    fullfile(repoRoot, 'graph', 'compute_file_sha256.m')
    fullfile(repoRoot, 'clustering', ...
        'build_run09_output_manifest.m')
    ];
end

function names = i_output_names(params)
names = [
    params.parameter_audit_file
    params.freeze_validation_file
    params.adapter_audit_file
    params.node_index_file
    params.backend_audit_file
    params.replicate_audit_file
    params.membership_mat_file];
if params.run_mode == "full"
    names = [names
        params.consensus_audit_file
        params.consensus_mat_file]; %#ok<AGROW>
end
names = [names
    params.partition_stability_file
    params.resolution_sensitivity_file
    params.node_ambiguity_file
    params.hierarchy_file
    params.hierarchy_resolution_audit_file
    params.candidate_definition_file
    params.node_membership_file
    params.candidate_topology_audit_file
    params.objective_concordance_audit_file
    params.scale_composition_audit_file
    params.session_composition_audit_file
    params.behavioral_profile_file
    params.event_profile_file
    params.exemplar_manifest_file
    params.annotation_file
    params.annotation_provenance_file
    params.qc_figure_manifest_file
    params.paper_figure_manifest_file
    ];
end

function roles = i_output_roles(params)
roles = [
    "effective_parameter_audit"
    "frozen_input_validation"
    "strict_graph_adapter_audit"
    "node_identity_index_audit"
    "weighted_backend_validation"
    "raw_partition_replicate_audit"
    "raw_partition_binary_audit"];
if params.run_mode == "full"
    roles = [roles
        "consensus_partition_audit"
        "consensus_partition_binary_audit"]; %#ok<AGROW>
end
roles = [roles
    "partition_stability_audit"
    "resolution_selection_audit"
    "retained_node_ambiguity_audit"
    "defines_frozen_candidate_hierarchy"
    "all_resolution_hierarchy_audit"
    "defines_frozen_candidate_ids"
    "defines_frozen_candidate_membership"
    "postfreeze_graph_topology_and_interpretation_eligibility_audit"
    "audit_only_cpm_modularity_concordance_and_topology"
    "postfit_scale_composition_audit"
    "postfit_session_composition_audit"
    "postfit_condition_blind_behavioral_profile"
    "postfit_condition_blind_event_profile"
    "postfit_graph_stability_centrality_exemplars"
    "postfit_provisional_social_non_social_annotation"
    "postfit_annotation_rule_and_source_provenance"
    "postfit_qc_figure_hash_manifest"
    "postfit_paper_figure_hash_manifest"
    ];
end

function Manifest = i_append_figure_manifest(Manifest, outRoot, ...
        manifestFile, recordType, freezeId, candidateFreezeId, ...
        configHash, sourceBundleHash)
figureManifestPath = fullfile(outRoot, manifestFile);
FigureManifest = readtable(figureManifestPath, ...
    'Delimiter', ',', 'TextType', 'string');
for i = 1:height(FigureManifest)
    pathText = FigureManifest.figure_path(i);
    assert(isfile(pathText), ...
        'build_run09_output_manifest:MissingFigure', ...
        'Figure recorded in %s is missing: %s', ...
        manifestFile, pathText);
    [hash, bytes] = compute_file_sha256(pathText);
    assert(hash == FigureManifest.sha256(i) && ...
        bytes == FigureManifest.file_bytes(i), ...
        'build_run09_output_manifest:FigureHashMismatch', ...
        'Figure hash or byte count differs from %s.', manifestFile);
    Manifest = [Manifest; i_row(recordType, ... %#ok<AGROW>
        FigureManifest.figure_id(i), pathText, ...
        FigureManifest.scientific_role(i), false, true, ...
        "present_hashed_byte_exact", NaN, bytes, hash, ...
        freezeId, candidateFreezeId, configHash, ...
        sourceBundleHash, true, ...
        ['Generated after freeze and never used for membership, ' ...
         'eligibility, exemplar selection, or annotation.'])];
end
end

function row = i_row(recordType, artifactId, artifactPath, artifactRole, ...
        definesMembership, postFit, status, rowCount, bytes, hash, ...
        freezeId, candidateFreezeId, configHash, sourceBundleHash, ...
        noLabels, details)
row = table(string(recordType), string(artifactId), string(artifactPath), ...
    string(artifactRole), logical(definesMembership), logical(postFit), ...
    string(status), double(rowCount), double(bytes), string(hash), ...
    string(freezeId), string(candidateFreezeId), string(configHash), ...
    string(sourceBundleHash), logical(noLabels), string(details), ...
    'VariableNames', i_manifest_names());
end

function T = i_empty_manifest()
T = table(strings(0, 1), strings(0, 1), strings(0, 1), ...
    strings(0, 1), false(0, 1), false(0, 1), strings(0, 1), ...
    nan(0, 1), nan(0, 1), strings(0, 1), strings(0, 1), ...
    strings(0, 1), strings(0, 1), strings(0, 1), false(0, 1), ...
    strings(0, 1), 'VariableNames', i_manifest_names());
end

function names = i_manifest_names()
names = {'record_type','artifact_id','artifact_path','artifact_role', ...
    'defines_membership','post_fit_interpretation','status','row_count', ...
    'file_bytes','sha256','run08_freeze_id','candidate_freeze_id', ...
    'config_sha256','source_bundle_sha256', ...
    'no_experimental_labels_used','details'};
end

function value = i_output_details(definesMembership, postFit)
if definesMembership
    value = "Frozen graph-derived candidate identity or hierarchy.";
elseif postFit
    value = "Post-fit interpretation or audit only; cannot revise membership or graph-only eligibility.";
else
    value = "Audit only and cannot revise frozen candidate membership.";
end
end

function rows = i_csv_rows(pathText)
if endsWith(lower(string(pathText)), ".csv")
    % Count parsed records rather than physical lines. Quoted CSV fields may
    % legally contain line feeds, so byte-level newline counts can overstate
    % the auditable table row count.
    rows = height(readtable(pathText, 'Delimiter', ','));
else
    rows = NaN;
end
end

function value = i_file_bytes(pathText)
info = dir(pathText);
value = info.bytes;
end
