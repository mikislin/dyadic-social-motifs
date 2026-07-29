function FigureManifest = make_run09_candidate_qc_figures( ...
        outRoot, Postfit, TopologyAudits, params)
%MAKE_RUN09_CANDIDATE_QC_FIGURES Create post-freeze audit figures.
%
% Figure appearance is never consumed by membership, eligibility, exemplar,
% or annotation logic.

outRoot = string(outRoot);
figureRoot = fullfile(outRoot, params.qc_figure_directory);
if ~isfolder(figureRoot)
    mkdir(figureRoot);
end

records = struct('id', {}, 'file', {}, 'role', {}, 'inputs', {});

Topology = TopologyAudits.candidate_topology;
fig = i_figure();
scatter(Topology.node_count, Topology.weighted_conductance, 48, ...
    double(Topology.eligible_for_behavioral_interpretation), 'filled');
set(gca, 'XScale', 'log');
xlabel('Candidate node count');
ylabel('Weighted conductance');
title('Frozen candidate topology and interpretation eligibility');
cb = colorbar;
cb.Ticks = [0 1];
cb.TickLabels = {'residual/unstable','eligible'};
grid on;
records(end + 1) = i_export(fig, figureRoot, ... %#ok<AGROW>
    "run09_qc_candidate_topology", ...
    "Frozen size-conductance topology with graph-only eligibility", ...
    "motif_candidate_topology_audit.csv", params);

fig = i_figure();
Scale = Postfit.scale_composition;
candidates = Postfit.candidate_ids(Postfit.eligible);
if isempty(candidates)
    candidates = Postfit.candidate_ids;
end
categories = unique(Scale.category_index, 'sorted');
Z = i_composition_matrix(Scale, candidates, categories);
imagesc(log2(max(Z, 1e-3)));
colormap(fig, parula);
colorbar;
xlabel('Scale index');
ylabel('Frozen candidate');
title('Post-fit scale enrichment (log_2 observed/overall)');
i_ticks(categories, candidates);
records(end + 1) = i_export(fig, figureRoot, ... %#ok<AGROW>
    "run09_qc_scale_composition", ...
    "Post-fit scale mixing; never used for membership or eligibility", ...
    "motif_candidate_scale_composition_audit.csv", params);

fig = i_figure();
Session = Postfit.session_composition;
[candidateSize, sessionEntropy] = i_candidate_composition_summary(Session, ...
    Postfit.candidate_ids);
scatter(candidateSize, sessionEntropy, 48, ...
    double(Postfit.eligible), 'filled');
set(gca, 'XScale', 'log');
xlabel('Candidate node count');
ylabel('Normalized session mixing entropy');
title('Post-fit session mixing audit');
ylim([0 1]);
grid on;
records(end + 1) = i_export(fig, figureRoot, ... %#ok<AGROW>
    "run09_qc_session_mixing", ...
    "Post-fit session mixing; never used for membership or eligibility", ...
    "motif_candidate_session_composition_audit.csv", params);

fig = i_figure();
[Z, featureNames] = i_profile_matrix(Postfit.behavioral_profile, ...
    candidates, params.qc_top_profile_features);
imagesc(Z);
colormap(fig, parula);
colorbar;
xlabel('Condition-blind behavioral feature');
ylabel('Frozen candidate');
title('Post-fit behavioral standardized mean differences');
i_profile_ticks(featureNames, candidates);
records(end + 1) = i_export(fig, figureRoot, ... %#ok<AGROW>
    "run09_qc_behavioral_profiles", ...
    "Condition-blind behavioral profiles generated after candidate freeze", ...
    "motif_candidate_behavioral_profile.csv", params);

fig = i_figure();
[Z, featureNames] = i_profile_matrix(Postfit.event_profile, ...
    candidates, params.qc_top_profile_features);
imagesc(Z);
colormap(fig, parula);
colorbar;
xlabel('Condition-blind event feature');
ylabel('Frozen candidate');
title('Post-fit event standardized mean differences');
i_profile_ticks(featureNames, candidates);
records(end + 1) = i_export(fig, figureRoot, ... %#ok<AGROW>
    "run09_qc_event_profiles", ...
    "Condition-blind event profiles generated after candidate freeze", ...
    "motif_candidate_event_profile.csv", params);

fig = i_figure();
Objective = TopologyAudits.objective_concordance;
challenger = Objective.analysis_role == ...
    "audit_only_modularity_challenger";
scatter(Objective.adjusted_rand_index_to_retained_cpm(challenger), ...
    Objective.node_weighted_mean_conductance(challenger), 50, ...
    Objective.resolution(challenger), 'filled');
xlabel('ARI to retained weighted CPM partition');
ylabel('Node-weighted mean conductance');
title('Audit-only weighted modularity challenger');
colorbar;
grid on;
records(end + 1) = i_export(fig, figureRoot, ... %#ok<AGROW>
    "run09_qc_objective_concordance", ...
    "CPM versus modularity concordance and topology comparison", ...
    "motif_candidate_objective_concordance_audit.csv", params);

n = numel(records);
figureId = strings(n, 1);
fileName = strings(n, 1);
figurePath = strings(n, 1);
scientificRole = strings(n, 1);
inputArtifacts = strings(n, 1);
fileBytes = zeros(n, 1);
sha256 = strings(n, 1);
for i = 1:n
    figureId(i) = records(i).id;
    fileName(i) = records(i).file;
    figurePath(i) = fullfile(figureRoot, records(i).file);
    scientificRole(i) = records(i).role;
    inputArtifacts(i) = records(i).inputs;
    info = dir(figurePath(i));
    fileBytes(i) = info.bytes;
    sha256(i) = compute_file_sha256(figurePath(i));
end
FigureManifest = table(figureId, fileName, figurePath, ...
    scientificRole, inputArtifacts, fileBytes, sha256, ...
    repmat(params.qc_figure_format, n, 1), ...
    repmat(true, n, 1), repmat(false, n, 1), ...
    repmat(false, n, 1), repmat(false, n, 1), ...
    repmat(string(Postfit.membership_sha256), n, 1), ...
    repmat("none", n, 1), ...
    'VariableNames', { ...
    'figure_id','file_name','figure_path','scientific_role', ...
    'input_artifacts','file_bytes','sha256','format', ...
    'generated_after_candidate_freeze','used_for_membership', ...
    'used_for_interpretation_eligibility','used_for_annotation', ...
    'frozen_membership_sha256','experimental_labels_used'});
end

function fig = i_figure()
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100 100 1200 760]);
axes(fig);
end

function record = i_export(fig, root, id, role, inputs, params)
fileName = id + "." + params.qc_figure_format;
pathText = fullfile(root, fileName);
exportgraphics(fig, pathText, 'Resolution', 160);
close(fig);
record = struct('id', id, 'file', fileName, ...
    'role', role, 'inputs', inputs);
end

function Z = i_composition_matrix(T, candidateIds, categories)
Z = zeros(numel(candidateIds), numel(categories));
for c = 1:numel(candidateIds)
    for k = 1:numel(categories)
        row = string(T.motif_candidate_id) == candidateIds(c) & ...
            T.category_index == categories(k);
        assert(nnz(row) == 1, ...
            'make_run09_candidate_qc_figures:BadCompositionCell', ...
            'Every candidate-category composition cell must occur once.');
        Z(c, k) = T.enrichment_ratio(row);
    end
end
end

function [sizes, entropy] = i_candidate_composition_summary(T, candidateIds)
sizes = zeros(numel(candidateIds), 1);
entropy = zeros(numel(candidateIds), 1);
for c = 1:numel(candidateIds)
    rows = string(T.motif_candidate_id) == candidateIds(c);
    sizes(c) = sum(T.node_count(rows));
    entropy(c) = T.normalized_mixing_entropy(find(rows, 1));
end
end

function [Z, selectedNames] = i_profile_matrix( ...
        Profile, candidateIds, maximumFeatures)
allNames = unique(string(Profile.feature_name), 'stable');
Full = zeros(numel(candidateIds), numel(allNames));
for c = 1:numel(candidateIds)
    for f = 1:numel(allNames)
        row = string(Profile.motif_candidate_id) == candidateIds(c) & ...
            string(Profile.feature_name) == allNames(f);
        assert(nnz(row) == 1, ...
            'make_run09_candidate_qc_figures:BadProfileCell', ...
            'Every candidate-feature profile cell must occur once.');
        Full(c, f) = Profile.standardized_mean_difference(row);
    end
end
[~, order] = sort(var(Full, 0, 1), 'descend');
order = order(1:min(maximumFeatures, numel(order)));
Z = Full(:, order);
selectedNames = allNames(order);
end

function i_ticks(categories, candidates)
ax = gca;
ax.XTick = 1:numel(categories);
ax.XTickLabel = compose('%.0f', categories);
ax.YTick = 1:numel(candidates);
ax.YTickLabel = erase(candidates, "MC_M000750_");
ax.TickLabelInterpreter = 'none';
end

function i_profile_ticks(features, candidates)
ax = gca;
ax.XTick = 1:numel(features);
ax.XTickLabel = features;
ax.XTickLabelRotation = 45;
ax.YTick = 1:numel(candidates);
ax.YTickLabel = erase(candidates, "MC_M000750_");
ax.TickLabelInterpreter = 'none';
end
