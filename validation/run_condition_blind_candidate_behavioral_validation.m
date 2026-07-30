function outputs = run_condition_blind_candidate_behavioral_validation( ...
        repoRoot, opts)
%RUN_CONDITION_BLIND_CANDIDATE_BEHAVIORAL_VALIDATION Registered run_10.
%
% The candidate-freeze validator is the first data-reading operation.
% Experimental labels and run_05-run_09 discovery values never enter
% independent measurement, folds, nulls, models, or validation status.

if nargin<1 || strlength(string(repoRoot))==0
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
end
if nargin<2 || isempty(opts)
    opts = struct();
end
if ~isfield(opts,'configPath') || strlength(string(opts.configPath))==0
    opts.configPath = fullfile(repoRoot,'config', ...
        'motif_candidate_behavioral_validation_config.csv');
end

repoRoot = string(repoRoot);
cd(repoRoot);
addpath(genpath(repoRoot));
params = load_motif_candidate_behavioral_validation_config( ...
    opts.configPath);
outRoot = string(resolve_repo_path(repoRoot,params.output_dir));
i_make_dir(outRoot);
i_make_dir(fullfile(outRoot,'logs'));
diary(fullfile(outRoot,'logs', ...
    'run_10_validate_motif_candidates_behaviorally_latest.log'));
diaryCleanup = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('run_10_validate_motif_candidates_behaviorally\n');
fprintf('Run mode: %s | phase: %s\n',params.run_mode,params.phase);
fprintf('Output root: %s\n',outRoot);
paths = i_paths(outRoot,params);
write_table_atomic(i_parameter_audit(params),paths.parameterAudit);

% Mandatory fail-closed gate. No pose, video, rating, provisional
% annotation, or candidate-specific validation feature is read above it.
Validation = validate_run10_candidate_freeze( ...
    repoRoot,params,paths.freezeValidation);
fprintf('Freeze passed: %s | membership %s\n', ...
    Validation.candidate_freeze_id,Validation.membership_sha256);

[Registry,OverlapAudit] = ...
    build_run10_validation_feature_registry( ...
    repoRoot,params,Validation);
write_table_atomic(Registry,paths.featureRegistry);
write_table_atomic(OverlapAudit,paths.overlapAudit);

Sample = build_run10_validation_sample_manifest( ...
    repoRoot,params,Validation,Registry);
write_table_atomic(Sample,paths.sampleManifest);
Measurement = extract_run10_independent_behavioral_measures( ...
    repoRoot,params,Validation,Registry,Sample);
write_table_atomic(Measurement,paths.measurement);
FoldManifest = build_run10_blocked_validation_folds( ...
    Measurement,params);
write_table_atomic(FoldManifest,paths.foldManifest);

Packet = prepare_run10_blinded_rating_packet( ...
    repoRoot,outRoot,params,Validation);
write_table_atomic(Packet.manifest,paths.ratingPacketManifest);
write_table_atomic(Packet.template,paths.ratingTemplate);
write_table_atomic(Packet.codebook,paths.ratingCodebook);
write_table_atomic(Packet.randomization,paths.ratingRandomization);
write_text_atomic(Packet.instructions,paths.ratingInstructions);
[Ratings,RatingAudit] = ingest_run10_blinded_ratings( ...
    outRoot,params,Packet);
write_table_atomic(RatingAudit,paths.ratingIngestionAudit);

if params.phase=="prepare"
    fprintf(['Preparation complete. Automated candidate analysis was not ' ...
        'requested in phase=prepare.\n']);
    outputs = i_outputs(outRoot,paths,Validation,Measurement,RatingAudit);
    return
end

[Coherence,CoherenceNull] = ...
    analyze_run10_candidate_behavioral_coherence( ...
    Measurement,Registry,Validation,params);
write_table_atomic(Coherence,paths.coherence);
[Distinctness,Heldout,Confusion,DiscriminabilityNull] = ...
    analyze_run10_candidate_discriminability( ...
    Measurement,Registry,FoldManifest,Validation,params);
write_table_atomic(Distinctness,paths.distinctness);
write_table_atomic(Heldout,paths.heldout);
write_table_atomic(Confusion,paths.confusion);
NullAudit = [CoherenceNull;DiscriminabilityNull];
write_table_atomic(NullAudit,paths.nullAudit);
CrossSession = analyze_run10_cross_session_reproducibility( ...
    Measurement,Registry,Validation,params);
write_table_atomic(CrossSession,paths.crossSession);
CrossScale = analyze_run10_cross_scale_consistency( ...
    Measurement,Registry,Validation,params);
write_table_atomic(CrossScale,paths.crossScale);

[ReviewSummary,RaterAgreement] = analyze_run10_rater_agreement( ...
    Ratings,Validation,params);
write_table_atomic(ReviewSummary,paths.blindedReviewSummary);
write_table_atomic(RaterAgreement,paths.raterAgreement);
Status = assign_run10_behavioral_validation_status( ...
    Validation,Coherence,Heldout,CrossSession,CrossScale, ...
    ReviewSummary,params);
write_table_atomic(Status,paths.status);

% Only after automated metrics, ratings, and status are frozen on disk may
% graph-to-behavior quantities be compared. They remain audit-only.
GraphBehavior = analyze_run10_graph_behavior_correspondence( ...
    Validation,Coherence,ReviewSummary,params);
write_table_atomic(GraphBehavior,paths.graphBehavior);
Provenance = build_run10_behavioral_validation_provenance( ...
    repoRoot,outRoot,params,Validation,RatingAudit);
write_table_atomic(Provenance,paths.provenance);

[QcFigures,PaperFigures] = ...
    make_run10_behavioral_validation_figures( ...
    outRoot,Registry,Measurement,Coherence,Heldout,Confusion, ...
    CrossSession,CrossScale,ReviewSummary,Status,GraphBehavior,params);
write_table_atomic(QcFigures,paths.qcFigureManifest);
write_table_atomic(PaperFigures,paths.paperFigureManifest);
OutputManifest = build_run10_output_manifest( ...
    repoRoot,outRoot,params,Validation,Measurement,RatingAudit, ...
    QcFigures,PaperFigures);
write_table_atomic(OutputManifest,paths.outputManifest);

fprintf('Selected nodes: %d | analyzed nodes: %d\n', ...
    nnz(Sample.sample_selected), ...
    nnz(Measurement.eligible_for_automated_analysis));
fprintf('Eligible candidates: %d | ratings: %s\n', ...
    numel(Validation.eligible_candidate_ids), ...
    RatingAudit.ingestion_status(1));
fprintf('Output manifest: %s\n',paths.outputManifest);
outputs = i_outputs(outRoot,paths,Validation,Measurement,RatingAudit);
outputs.status_path = paths.status;
outputs.output_manifest_path = paths.outputManifest;
outputs.n_behaviorally_supported = nnz( ...
    Status.run10_validation_status=="behaviorally_supported");
end

function paths = i_paths(outRoot,params)
paths = struct();
paths.parameterAudit = fullfile(outRoot,params.parameter_audit_file);
paths.freezeValidation = fullfile(outRoot,params.freeze_validation_file);
paths.featureRegistry = fullfile(outRoot,params.feature_registry_file);
paths.overlapAudit = fullfile(outRoot,params.overlap_audit_file);
paths.sampleManifest = fullfile(outRoot,params.sample_manifest_file);
paths.foldManifest = fullfile(outRoot,params.fold_manifest_file);
paths.measurement = fullfile(outRoot,params.measurement_file);
paths.coherence = fullfile(outRoot,params.coherence_file);
paths.distinctness = fullfile(outRoot,params.distinctness_file);
paths.crossSession = fullfile(outRoot,params.cross_session_file);
paths.crossScale = fullfile(outRoot,params.cross_scale_file);
paths.heldout = fullfile(outRoot,params.heldout_file);
paths.confusion = fullfile(outRoot,params.confusion_file);
paths.nullAudit = fullfile(outRoot,params.null_audit_file);
paths.ratingPacketManifest = fullfile(outRoot, ...
    params.rating_packet_manifest_file);
paths.ratingTemplate = fullfile(outRoot,params.rating_template_file);
paths.ratingCodebook = fullfile(outRoot,params.rating_codebook_file);
paths.ratingRandomization = fullfile(outRoot, ...
    params.rating_randomization_file);
paths.ratingInstructions = fullfile(outRoot, ...
    params.rating_instructions_file);
paths.ratingIngestionAudit = fullfile(outRoot, ...
    params.rating_ingestion_audit_file);
paths.blindedReviewSummary = fullfile(outRoot, ...
    params.blinded_review_summary_file);
paths.raterAgreement = fullfile(outRoot,params.rater_agreement_file);
paths.graphBehavior = fullfile(outRoot,params.graph_behavior_file);
paths.status = fullfile(outRoot,params.status_file);
paths.provenance = fullfile(outRoot,params.provenance_file);
paths.qcFigureManifest = fullfile(outRoot, ...
    params.qc_figure_manifest_file);
paths.paperFigureManifest = fullfile(outRoot, ...
    params.paper_figure_manifest_file);
paths.outputManifest = fullfile(outRoot,params.output_manifest_file);
end

function T = i_parameter_audit(params)
C = params.config_table;
T = table(string(C.parameter),string(C.effective_value), ...
    string(C.type),string(C.env_override),logical(C.env_override_used), ...
    string(C.description), ...
    'VariableNames', {'parameter','effective_value','type', ...
    'env_override','env_override_used','description'});
T.used_for_validation_status = ismember(T.parameter,[ ...
    "minimum_nodes_primary_full","minimum_nodes_primary_smoke", ...
    "minimum_sessions_primary","coherence_min_relative_improvement", ...
    "coherence_max_fdr_q","heldout_min_candidate_recall", ...
    "heldout_min_macro_balanced_accuracy","heldout_null_quantile", ...
    "cross_session_min_split_profile_correlation", ...
    "cross_session_min_largest_session_removed_correlation", ...
    "cross_session_max_largest_session_fraction", ...
    "supported_required_primary_gates","partial_required_primary_gates", ...
    "validation_status_rule_version"]);
T.experimental_labels_used = repmat("none",height(T),1);
end

function outputs = i_outputs(outRoot,paths,Validation,Measurement,RatingAudit)
outputs = struct();
outputs.output_root = outRoot;
outputs.freeze_validation_path = paths.freezeValidation;
outputs.feature_registry_path = paths.featureRegistry;
outputs.sample_manifest_path = paths.sampleManifest;
outputs.measurement_path = paths.measurement;
outputs.fold_manifest_path = paths.foldManifest;
outputs.rating_packet_manifest_path = paths.ratingPacketManifest;
outputs.rating_ingestion_audit_path = paths.ratingIngestionAudit;
outputs.candidate_freeze_id = Validation.candidate_freeze_id;
outputs.membership_sha256 = Validation.membership_sha256;
outputs.n_measurements = height(Measurement);
outputs.rating_status = string(RatingAudit.ingestion_status(1));
end

function i_make_dir(pathText)
if ~isfolder(pathText)
    mkdir(pathText);
end
end
