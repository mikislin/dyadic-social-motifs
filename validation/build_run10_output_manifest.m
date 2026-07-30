function Manifest = build_run10_output_manifest( ...
        repoRoot, outRoot, params, Validation, Measurement, ...
        RatingAudit, QcFigures, PaperFigures)
%BUILD_RUN10_OUTPUT_MANIFEST Hash run_10 inputs, sources, outputs, figures.

repoRoot = string(repoRoot);
outRoot = string(outRoot);
configHash = compute_file_sha256(params.config_path);
Manifest = i_empty();

inputPaths = [ ...
    Validation.paths.membership
    Validation.paths.definition
    Validation.paths.hierarchy
    Validation.paths.topology
    Validation.paths.ambiguity
    Validation.paths.manifest
    Validation.paths.run09FreezeValidation
    string(resolve_repo_path(repoRoot,params.run08_node_manifest_file))
    string(resolve_repo_path(repoRoot,params.preprocess_qc_file))
    fullfile(resolve_repo_path(repoRoot,params.pre_run10_review_dir), ...
        'pre_run10_blinded_video_review_output_manifest.csv')];
inputRoles = [ ...
    "frozen_run09_candidate_membership"
    "frozen_run09_candidate_ids"
    "frozen_run09_flat_hierarchy"
    "frozen_run09_graph_only_eligibility"
    "frozen_run09_node_ambiguity_audit"
    "run09_authoritative_output_manifest"
    "run09_input_freeze_provenance"
    "safe_column_node_session_scale_anchor_provenance"
    "safe_column_pose_file_provenance"
    "hashed_blinded_review_preparation_manifest"];
for i = 1:numel(inputPaths)
    Manifest = [Manifest;i_file_row("frozen_or_audit_input", ...
        i_name(inputPaths(i)),inputPaths(i),inputRoles(i), ...
        params,configHash)]; %#ok<AGROW>
end

poseFiles = unique(string(Measurement.preprocess_output_file));
for i = 1:numel(poseFiles)
    absolute = string(resolve_repo_path(repoRoot,poseFiles(i)));
    Manifest = [Manifest;i_file_row("pose_measurement_input", ...
        i_name(absolute),absolute, ...
        "cleaned_pose_disjoint_landmark_measurement_input", ...
        params,configHash)]; %#ok<AGROW>
end

sourcePaths = i_source_paths(repoRoot,params);
sourceHashes = strings(numel(sourcePaths),1);
for i = 1:numel(sourcePaths)
    sourceHashes(i) = compute_file_sha256(sourcePaths(i));
end
sourceBundle = compute_run09_sha256_text(strjoin( ...
    sourcePaths+"|"+sourceHashes,newline));
for i = 1:numel(sourcePaths)
    row = i_file_row("source_provenance",i_name(sourcePaths(i)), ...
        sourcePaths(i),"run10_implementation_source",params,configHash);
    row.source_bundle_sha256 = sourceBundle;
    Manifest = [Manifest;row]; %#ok<AGROW>
end

outputNames = [ ...
    params.parameter_audit_file
    params.freeze_validation_file
    params.feature_registry_file
    params.overlap_audit_file
    params.sample_manifest_file
    params.fold_manifest_file
    params.measurement_file
    params.coherence_file
    params.distinctness_file
    params.cross_session_file
    params.cross_scale_file
    params.heldout_file
    params.confusion_file
    params.null_audit_file
    params.rating_packet_manifest_file
    params.rating_template_file
    params.rating_codebook_file
    params.rating_randomization_file
    params.rating_instructions_file
    params.rating_ingestion_audit_file
    params.blinded_review_summary_file
    params.rater_agreement_file
    params.graph_behavior_file
    params.status_file
    params.provenance_file
    params.qc_figure_manifest_file
    params.paper_figure_manifest_file];
outputRoles = [ ...
    "effective_parameter_audit"
    "candidate_freeze_validation"
    "independent_feature_lineage_registry"
    "discovery_feature_exclusion_audit"
    "value_blind_validation_sample_manifest"
    "grouped_session_fold_manifest"
    "automated_independent_behavioral_measurement"
    "automated_within_candidate_coherence"
    "automated_between_candidate_distinctness"
    "automated_cross_session_reproducibility"
    "audit_only_cross_scale_validation"
    "automated_grouped_heldout_discriminability"
    "automated_validation_confusion_matrix"
    "automated_blocked_null_audit"
    "blinded_human_review_packet"
    "blank_human_rating_template"
    "human_rating_codebook"
    "sealed_human_review_randomization"
    "human_review_instructions"
    "human_rating_ingestion_audit"
    "blinded_human_review_summary"
    "blinded_human_rater_agreement"
    "audit_only_graph_behavior_correspondence"
    "run10_validation_status"
    "run10_validation_provenance"
    "qc_figure_hash_manifest"
    "paper_figure_hash_manifest"];
for i = 1:numel(outputNames)
    pathText = fullfile(outRoot,outputNames(i));
    Manifest = [Manifest;i_file_row("output_artifact",outputNames(i), ...
        pathText,outputRoles(i),params,configHash)]; %#ok<AGROW>
end

Figures = [QcFigures;PaperFigures];
for i = 1:height(Figures)
    Manifest = [Manifest;i_file_row("figure",Figures.figure_id(i), ...
        Figures.figure_path(i),Figures.scientific_role(i), ...
        params,configHash)]; %#ok<AGROW>
end

contracts = [ ...
    "candidate_freeze_id",params.expected_candidate_freeze_id
    "membership_sha256",params.expected_membership_sha256
    "feature_panel_version",params.validation_feature_panel_version
    "feature_lineage_registry_sha256",compute_file_sha256( ...
        fullfile(outRoot,params.feature_registry_file))
    "classifier_and_parameters",params.classifier_id+"|"+ ...
        params.standardization_rule+"|fixed_hyperparameters"
    "null_construction","candidate_target_permuted_within_session_by_scale"
    "fold_assignment_sha256",compute_file_sha256( ...
        fullfile(outRoot,params.fold_manifest_file))
    "permutation_seed_and_count",string(params.permutation_seed)+"|"+ ...
        string(params.active_permutation_count)
    "validation_status_rule_version",params.validation_status_rule_version
    "rating_packet_sha256",compute_file_sha256( ...
        fullfile(outRoot,params.rating_packet_manifest_file))
    "rating_input_status",string(RatingAudit.ingestion_status(1))
    "rating_input_sha256",string(RatingAudit.rating_input_sha256(1))
    "experimental_labels_used","none"
    "membership_mutation_policy","prohibited"
    "run09_profiles_and_annotations_as_independent_evidence","prohibited"
    "same_dataset_external_replication_claim","prohibited"
    ];
for i = 1:size(contracts,1)
    Manifest = [Manifest;i_contract_row(contracts(i,1),contracts(i,2), ...
        params,configHash,sourceBundle)]; %#ok<AGROW>
end

Manifest.source_bundle_sha256( ...
    Manifest.source_bundle_sha256=="") = sourceBundle;
Manifest.matlab_version(:) = string(version);
v = ver;
toolboxText = strjoin(string({v.Name})+"="+string({v.Version}),";");
Manifest.toolbox_versions(:) = toolboxText;
assert(~any(Manifest.record_type=="output_artifact" & ...
    contains(lower(Manifest.artifact_id),"membership")), ...
    'build_run10_output_manifest:ReplacementMembershipOutput', ...
    'Run_10 must not write a membership artifact.');
end

function paths = i_source_paths(repoRoot,params)
files = [ ...
    dir(fullfile(repoRoot,'validation','*run10*.m'))
    dir(fullfile(repoRoot,'validation', ...
        '*motif_candidate_behavioral_validation*.m'))
    dir(fullfile(repoRoot,'validation','*behavioral_validation*.m'))];
paths = string(fullfile({files.folder},{files.name}))';
paths = unique(paths,'stable');
paths = [paths
    fullfile(repoRoot,'paper','run_10_validate_motif_candidates_behaviorally.m')
    string(params.config_path)
    fullfile(repoRoot,'docs','run10_behavioral_validation_design_decision.md')
    fullfile(repoRoot,'io','write_table_atomic.m')
    fullfile(repoRoot,'io','write_text_atomic.m')
    fullfile(repoRoot,'graph','compute_file_sha256.m')];
paths = unique(paths,'stable');
assert(all(isfile(paths)), ...
    'build_run10_output_manifest:MissingSource', ...
    'A contributing run_10 source file is missing.');
end

function row = i_file_row(recordType,id,pathText,role,params,configHash)
assert(isfile(pathText), ...
    'build_run10_output_manifest:MissingArtifact', ...
    'Required manifest artifact is missing: %s',pathText);
[hash,bytes] = compute_file_sha256(pathText);
rows = NaN;
if endsWith(lower(string(pathText)),'.csv')
    rows = height(readtable(pathText,'Delimiter',','));
end
row = table(string(recordType),string(id),string(pathText), ...
    string(role),"present_hashed_byte_exact",rows,bytes,hash, ...
    params.expected_candidate_freeze_id, ...
    params.expected_membership_sha256,configHash,"","", ...
    "","","",true,false,"none", ...
    'VariableNames',i_names());
end

function row = i_contract_row(id,value,params,configHash,sourceBundle)
row = table("run_contract",string(id),"","reviewer_audit_contract", ...
    "locked",NaN,NaN,"",params.expected_candidate_freeze_id, ...
    params.expected_membership_sha256,configHash,sourceBundle,"","","", ...
    string(value),true,false,"none",'VariableNames',i_names());
end

function T = i_empty()
T = table(strings(0,1),strings(0,1),strings(0,1),strings(0,1), ...
    strings(0,1),nan(0,1),nan(0,1),strings(0,1),strings(0,1), ...
    strings(0,1),strings(0,1),strings(0,1),strings(0,1), ...
    strings(0,1),strings(0,1),strings(0,1),false(0,1),false(0,1), ...
    strings(0,1),'VariableNames',i_names());
end

function names = i_names()
names = {'record_type','artifact_id','artifact_path','artifact_role', ...
    'status','row_count','file_bytes','sha256','candidate_freeze_id', ...
    'frozen_membership_sha256','config_sha256', ...
    'source_bundle_sha256','matlab_version','toolbox_versions','details', ...
    'contract_value','no_experimental_labels_used', ...
    'may_change_membership','experimental_labels_used'};
end

function value = i_name(pathText)
[~,stem,ext] = fileparts(pathText);
value = string(stem)+string(ext);
end
