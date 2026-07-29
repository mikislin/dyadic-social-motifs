function outputs = refresh_run09_postfit_outputs(repoRoot)
%REFRESH_RUN09_POSTFIT_OUTPUTS Regenerate only run_09 post-fit work.
%
% This entry point validates the run_08 freeze and the issued run_09
% membership hash, reads existing frozen partition artifacts, and generates
% only downstream audits, annotations, figures, and the output manifest.
% It never initializes or calls matlab-igraph.

if nargin < 1 || strlength(string(repoRoot)) == 0
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
end
repoRoot = string(repoRoot);
addpath(genpath(repoRoot));
params = load_motif_candidate_discovery_config();
assert(params.run_mode == "full", ...
    'refresh_run09_postfit_outputs:FullOnly', ...
    'Post-fit resume requires the issued full candidate freeze.');
outRoot = string(resolve_repo_path(repoRoot, params.output_dir));

Validation = validate_run09_frozen_handoff(repoRoot, params, ...
    fullfile(outRoot, params.freeze_validation_file));
Adapter = build_run09_frozen_graph_adapter(Validation);
Analysis = struct();
Analysis.node_membership = i_read(outRoot, ...
    params.node_membership_file);
Analysis.node_ambiguity = i_read(outRoot, ...
    params.node_ambiguity_file);
TopologyAudits = struct();
TopologyAudits.candidate_topology = i_read(outRoot, ...
    params.candidate_topology_audit_file);
TopologyAudits.objective_concordance = i_read(outRoot, ...
    params.objective_concordance_audit_file);

Postfit = build_run09_postfit_candidate_audits( ...
    repoRoot, outRoot, Validation, Adapter, Analysis, ...
    TopologyAudits, params);
i_write(Postfit.scale_composition, outRoot, ...
    params.scale_composition_audit_file);
i_write(Postfit.session_composition, outRoot, ...
    params.session_composition_audit_file);
i_write(Postfit.behavioral_profile, outRoot, ...
    params.behavioral_profile_file);
i_write(Postfit.event_profile, outRoot, params.event_profile_file);
i_write(Postfit.exemplar_manifest, outRoot, ...
    params.exemplar_manifest_file);
i_write(Postfit.annotation, outRoot, params.annotation_file);
i_write(Postfit.annotation_provenance, outRoot, ...
    params.annotation_provenance_file);

FigureManifest = make_run09_candidate_qc_figures( ...
    outRoot, Postfit, TopologyAudits, params);
i_write(FigureManifest, outRoot, params.qc_figure_manifest_file);
PaperFigureManifest = make_run09_candidate_paper_figures( ...
    outRoot, Postfit, TopologyAudits, params);
i_write(PaperFigureManifest, outRoot, ...
    params.paper_figure_manifest_file);
[Manifest, candidateFreezeId] = build_run09_output_manifest( ...
    repoRoot, outRoot, params, Validation);
i_write(Manifest, outRoot, params.output_manifest_file);

outputs = struct();
outputs.output_root = outRoot;
outputs.candidate_freeze_id = candidateFreezeId;
outputs.membership_sha256 = Postfit.membership_sha256;
outputs.n_candidates = numel(Postfit.candidate_ids);
outputs.n_scale_rows = height(Postfit.scale_composition);
outputs.n_session_rows = height(Postfit.session_composition);
outputs.n_behavior_rows = height(Postfit.behavioral_profile);
outputs.n_event_rows = height(Postfit.event_profile);
outputs.n_exemplars = height(Postfit.exemplar_manifest);
outputs.n_annotations = height(Postfit.annotation);
outputs.n_qc_figures = height(FigureManifest);
outputs.n_paper_figure_files = height(PaperFigureManifest);
end

function T = i_read(root, fileName)
T = readtable(fullfile(root, fileName), ...
    'Delimiter', ',', 'TextType', 'string');
end

function i_write(T, root, fileName)
write_run09_table_atomic(T, fullfile(root, fileName));
end
