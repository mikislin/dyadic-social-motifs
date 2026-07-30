function Provenance = build_run10_behavioral_validation_provenance( ...
        repoRoot, outRoot, params, Validation, RatingAudit)
%BUILD_RUN10_BEHAVIORAL_VALIDATION_PROVENANCE Reviewer-facing rule record.

repoRoot = string(repoRoot);
outRoot = string(outRoot);
ids = [ ...
    "candidate_freeze_id"
    "frozen_membership_sha256"
    "validation_feature_panel_version"
    "validation_status_rule_version"
    "rating_schema_version"
    "independence_definition"
    "candidate_membership_policy"
    "run09_eligibility_policy"
    "experimental_label_policy"
    "discovery_profile_policy"
    "session_usage"
    "scale_usage"
    "arena_usage"
    "human_rating_status"
    "same_dataset_replication_language"
    ];
values = [ ...
    params.expected_candidate_freeze_id
    params.expected_membership_sha256
    params.validation_feature_panel_version
    params.validation_status_rule_version
    params.rating_schema_version
    "disjoint raw-pose measurement channels plus blinded review and grouped-session generalization; not external cohort replication"
    "immutable; no merge split rename reassignment hierarchy or replacement membership output"
    "frozen graph-only input; run_10 status cannot change it"
    "condition cohort group drug genotype perturbation treatment outcome and response excluded"
    "run_06 events run_07 embedding PCA UMAP and run_09 profiles annotations are not independent validation evidence"
    "blocking sampling diversity and reproducibility only"
    "stratified null and postfreeze audit only"
    "not used"
    string(RatingAudit.ingestion_status(1))
    "independent behavioral validation within the same tracked dataset unless an external cohort is later supplied"
    ];
sourcePath = repmat(string(params.config_path),numel(ids),1);
sourcePath(ids=="candidate_freeze_id" | ...
    ids=="frozen_membership_sha256") = ...
    fullfile(resolve_repo_path(repoRoot,params.run09_output_dir), ...
    params.run09_manifest_file);
sourcePath(ids=="human_rating_status") = ...
    fullfile(outRoot,params.rating_ingestion_audit_file);
sha256 = strings(numel(ids),1);
for i = 1:numel(ids)
    if isfile(sourcePath(i))
        sha256(i) = compute_file_sha256(sourcePath(i));
    end
end
Provenance = table(repmat("run10_validation_contract",numel(ids),1), ...
    ids,values,sourcePath,sha256, ...
    repmat(Validation.candidate_freeze_id,numel(ids),1), ...
    repmat(Validation.membership_sha256,numel(ids),1), ...
    repmat(true,numel(ids),1),repmat(false,numel(ids),1), ...
    repmat("none",numel(ids),1), ...
    'VariableNames', {'record_type','provenance_id','value', ...
    'source_path','sha256','candidate_freeze_id', ...
    'frozen_membership_sha256','recorded_after_freeze_validation', ...
    'may_change_membership','experimental_labels_used'});
end
