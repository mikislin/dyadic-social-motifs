function FigureManifest = make_run09_candidate_paper_figures( ...
        outRoot, Postfit, TopologyAudits, params)
%MAKE_RUN09_CANDIDATE_PAPER_FIGURES Create paper-facing run_09 figures.
%
% These figures summarize frozen graph results and post-fit characterization.
% Figure appearance is never consumed by membership, eligibility, exemplar
% selection, annotation, or any later experimental-label analysis.

outRoot = string(outRoot);
figureRoot = fullfile(outRoot, params.paper_figure_directory);
if ~isfolder(figureRoot)
    mkdir(figureRoot);
end

Resolution = i_read(outRoot, params.resolution_sensitivity_file);
Hierarchy = i_read(outRoot, params.hierarchy_resolution_audit_file);
Topology = TopologyAudits.candidate_topology;
Objective = TopologyAudits.objective_concordance;
records = struct('id', {}, 'file', {}, 'role', {}, 'inputs', {}, ...
    'format', {}, 'paperRole', {});

fig = i_partition_validation_figure( ...
    Resolution, Hierarchy, Topology, Objective, params);
records = [records; i_export_pair(fig, figureRoot, ... %#ok<AGROW>
    "run09_paper_partition_validation", ...
    ['Multiresolution retention, consensus stability, candidate topology, ' ...
     'failed nesting, and objective-challenger audit'], ...
    strjoin([params.resolution_sensitivity_file, ...
    params.hierarchy_resolution_audit_file, ...
    params.candidate_topology_audit_file, ...
    params.objective_concordance_audit_file], ';'), ...
    "primary_methodological_figure", params)];

fig = i_candidate_organization_figure(Postfit, Topology);
records = [records; i_export_pair(fig, figureRoot, ... %#ok<AGROW>
    "run09_paper_candidate_organization", ...
    ['Graph-sample prevalence, multiscale organization, session mixing, ' ...
     'and provisional post-freeze annotation'], ...
    strjoin([params.candidate_topology_audit_file, ...
    params.scale_composition_audit_file, ...
    params.session_composition_audit_file, params.annotation_file], ';'), ...
    "primary_candidate_summary_figure", params)];

fig = i_postfit_profile_figure(Postfit);
records = [records; i_export_pair(fig, figureRoot, ... %#ok<AGROW>
    "run09_paper_condition_blind_profiles", ...
    ['Condition-blind behavioral and event characterization after the ' ...
     'candidate freeze; not independent behavioral validation'], ...
    strjoin([params.behavioral_profile_file, ...
    params.event_profile_file, params.annotation_file], ';'), ...
    "supplementary_postfit_characterization", params)];

n = numel(records);
figureId = strings(n, 1);
fileName = strings(n, 1);
figurePath = strings(n, 1);
scientificRole = strings(n, 1);
inputArtifacts = strings(n, 1);
format = strings(n, 1);
recommendedPaperRole = strings(n, 1);
fileBytes = zeros(n, 1);
sha256 = strings(n, 1);
for i = 1:n
    figureId(i) = records(i).id;
    fileName(i) = records(i).file;
    figurePath(i) = fullfile(figureRoot, records(i).file);
    scientificRole(i) = records(i).role;
    inputArtifacts(i) = records(i).inputs;
    format(i) = records(i).format;
    recommendedPaperRole(i) = records(i).paperRole;
    info = dir(figurePath(i));
    fileBytes(i) = info.bytes;
    sha256(i) = compute_file_sha256(figurePath(i));
end
FigureManifest = table(figureId, fileName, figurePath, ...
    scientificRole, inputArtifacts, format, recommendedPaperRole, ...
    fileBytes, sha256, true(n, 1), false(n, 1), false(n, 1), ...
    false(n, 1), repmat(string(Postfit.membership_sha256), n, 1), ...
    repmat("none", n, 1), ...
    'VariableNames', {'figure_id', 'file_name', 'figure_path', ...
    'scientific_role', 'input_artifacts', 'format', ...
    'recommended_paper_role', 'file_bytes', 'sha256', ...
    'generated_after_candidate_freeze', 'used_for_membership', ...
    'used_for_interpretation_eligibility', 'used_for_annotation', ...
    'frozen_membership_sha256', 'experimental_labels_used'});
end

function fig = i_partition_validation_figure( ...
        Resolution, Hierarchy, Topology, Objective, params)
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [50 50 1800 1100]);
layout = tiledlayout(fig, 2, 3, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
retained = logical(Resolution.retained_by_predeclared_rules);
x = Resolution.density_multiplier;

ax = nexttile(layout);
semilogx(ax, x, Resolution.representative_candidate_count, ...
    '-o', 'Color', [0.25 0.25 0.25], 'MarkerFaceColor', [0.7 0.7 0.7]);
hold(ax, 'on');
scatter(ax, x(retained), ...
    Resolution.representative_candidate_count(retained), 90, ...
    [0.00 0.45 0.70], 'filled', 'MarkerEdgeColor', 'k');
hold(ax, 'off');
xlabel(ax, 'CPM density multiplier');
ylabel(ax, 'Consensus candidate count');
title(ax, '(A) Resolution-dependent fragmentation');
grid(ax, 'on');

ax = nexttile(layout);
semilogx(ax, x, Resolution.minimum_consensus_split_ari, ...
    '-o', 'Color', [0.00 0.45 0.70], 'DisplayName', 'Minimum split ARI');
hold(ax, 'on');
semilogx(ax, x, Resolution.mean_consensus_split_node_stability, ...
    '-s', 'Color', [0.80 0.40 0.00], ...
    'DisplayName', 'Mean node stability');
yline(ax, params.consensus_retention_min_split_ari, '--', ...
    'Color', [0.35 0.35 0.35], 'DisplayName', 'Retention threshold');
scatter(ax, x(retained), ...
    Resolution.minimum_consensus_split_ari(retained), 80, ...
    [0.00 0.45 0.70], 'filled', 'MarkerEdgeColor', 'k', ...
    'HandleVisibility', 'off');
hold(ax, 'off');
xlabel(ax, 'CPM density multiplier');
ylabel(ax, 'Stability');
ylim(ax, [0 1]);
title(ax, '(B) Split-consensus reproducibility');
legend(ax, 'Location', 'southwest');
grid(ax, 'on');

ax = nexttile(layout);
eligible = logical(Topology.eligible_for_behavioral_interpretation);
scatter(ax, Topology.node_count(~eligible), ...
    Topology.weighted_conductance(~eligible), 36, ...
    [0.65 0.65 0.65], 'filled', 'DisplayName', 'Residual/unstable');
hold(ax, 'on');
scatter(ax, Topology.node_count(eligible), ...
    Topology.weighted_conductance(eligible), 66, ...
    [0.00 0.45 0.70], 'filled', 'MarkerEdgeColor', 'k', ...
    'DisplayName', 'Interpretation-eligible');
i_label_candidates(ax, Topology, eligible, ...
    Topology.node_count, Topology.weighted_conductance);
hold(ax, 'off');
set(ax, 'XScale', 'log');
xlabel(ax, 'Graph node count');
ylabel(ax, 'Weighted conductance');
title(ax, '(C) Candidate topology');
legend(ax, 'Location', 'southwest');
grid(ax, 'on');

ax = nexttile(layout);
sizeScale = 20 + 70 * sqrt(Topology.node_count ./ ...
    max(Topology.node_count));
scatter(ax, Topology.mean_node_stability(~eligible), ...
    Topology.boundary_node_fraction(~eligible), ...
    sizeScale(~eligible), [0.65 0.65 0.65], 'filled');
hold(ax, 'on');
scatter(ax, Topology.mean_node_stability(eligible), ...
    Topology.boundary_node_fraction(eligible), ...
    sizeScale(eligible), [0.00 0.45 0.70], 'filled', ...
    'MarkerEdgeColor', 'k');
xline(ax, params.interpretation_min_mean_node_stability, '--', ...
    'Color', [0.35 0.35 0.35]);
yline(ax, params.interpretation_max_boundary_node_fraction, '--', ...
    'Color', [0.35 0.35 0.35]);
i_label_candidates(ax, Topology, ...
    eligible | Topology.node_count > params.interpretation_small_candidate_max_nodes, ...
    Topology.mean_node_stability, Topology.boundary_node_fraction);
hold(ax, 'off');
xlabel(ax, 'Mean node stability');
ylabel(ax, 'Boundary-node fraction');
xlim(ax, [0 1]);
ylim(ax, [0 1]);
title(ax, '(D) Candidate uncertainty');
grid(ax, 'on');

ax = nexttile(layout);
pairKey = string(Hierarchy.coarse_resolution_id) + "->" + ...
    string(Hierarchy.fine_resolution_id);
[~, firstPairRow] = unique(pairKey, 'stable');
pairRows = Hierarchy(firstPairRow, :);
pairIndex = (1:height(pairRows))';
plot(ax, pairIndex, ...
    pairRows.pair_weighted_child_purity, '-o', ...
    'Color', [0.00 0.45 0.70], 'DisplayName', 'Weighted child purity');
hold(ax, 'on');
plot(ax, pairIndex, ...
    pairRows.pair_linked_node_fraction, '-s', ...
    'Color', [0.80 0.40 0.00], 'DisplayName', 'Linked-node fraction');
yline(ax, params.hierarchy_min_pair_weighted_child_purity, '--', ...
    'Color', [0.35 0.35 0.35], 'DisplayName', 'Required threshold');
hold(ax, 'off');
xlabel(ax, 'Adjacent resolution pair');
ylabel(ax, 'Nesting support');
ylim(ax, [0 1]);
title(ax, '(E) Hierarchy audit (no multilevel hierarchy retained)');
legend(ax, 'Location', 'southwest');
grid(ax, 'on');

ax = nexttile(layout);
isReference = Objective.analysis_role == "frozen_primary_reference";
isMedoid = logical(Objective.is_within_resolution_ari_medoid);
challenger = Objective.analysis_role == ...
    "audit_only_modularity_challenger" & isMedoid;
scatter(ax, Objective.candidate_count(challenger), ...
    Objective.node_weighted_mean_conductance(challenger), 70, ...
    [0.80 0.40 0.00], 's', 'filled', ...
    'DisplayName', 'Modularity challenger');
hold(ax, 'on');
scatter(ax, Objective.candidate_count(isReference), ...
    Objective.node_weighted_mean_conductance(isReference), 90, ...
    [0.00 0.45 0.70], 'o', 'filled', 'MarkerEdgeColor', 'k', ...
    'DisplayName', 'Frozen CPM');
rows = find(challenger | isReference);
for i = rows'
    if isReference(i)
        label = "CPM";
    else
        label = "\gamma=" + compose('%.1g', Objective.resolution(i));
    end
    text(ax, Objective.candidate_count(i), ...
        Objective.node_weighted_mean_conductance(i), " " + label, ...
        'FontSize', 9, 'Interpreter', 'tex');
end
hold(ax, 'off');
xlabel(ax, 'Candidate count');
ylabel(ax, 'Node-weighted mean conductance');
title(ax, '(F) Objective challenger');
legend(ax, 'Location', 'best');
grid(ax, 'on');

sgtitle(layout, ...
    'Condition-blind motif-candidate discovery and validation', ...
    'FontWeight', 'bold');
end

function fig = i_candidate_organization_figure(Postfit, Topology)
candidateIds = Postfit.candidate_ids(Postfit.eligible);
shortIds = i_short_candidate_ids(candidateIds);
nCandidates = numel(candidateIds);
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [50 50 1800 1050]);
layout = tiledlayout(fig, 2, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

[present, topologyLoc] = ismember(candidateIds, ...
    string(Topology.motif_candidate_id));
assert(all(present), ...
    'make_run09_candidate_paper_figures:MissingTopologyCandidate', ...
    'Every eligible paper candidate requires a topology row.');
eligibleTopology = Topology(topologyLoc, :);

ax = nexttile(layout);
barh(ax, 100 * eligibleTopology.node_prevalence, ...
    'FaceColor', [0.00 0.45 0.70], 'EdgeColor', 'none');
set(ax, 'YTick', 1:nCandidates, 'YTickLabel', shortIds, ...
    'YDir', 'reverse');
xlabel(ax, 'Graph-sample prevalence (%)');
ylabel(ax, 'Frozen candidate');
title(ax, '(A) Candidate prevalence');
grid(ax, 'on');

ax = nexttile(layout);
Scale = Postfit.scale_composition;
scaleRows = string(Scale.motif_candidate_id) == candidateIds(1);
scaleSeconds = i_numeric(Scale.category_label(scaleRows));
scaleIndex = Scale.category_index(scaleRows);
[scaleIndex, scaleOrder] = sort(scaleIndex);
scaleSeconds = scaleSeconds(scaleOrder);
Z = zeros(nCandidates, numel(scaleIndex));
for iCandidate = 1:nCandidates
    rows = string(Scale.motif_candidate_id) == candidateIds(iCandidate);
    candidateScale = Scale(rows, :);
    [present, loc] = ismember(scaleIndex, candidateScale.category_index);
    assert(all(present), ...
        'make_run09_candidate_paper_figures:MissingScaleCell', ...
        'Every eligible candidate requires every scale cell.');
    Z(iCandidate, :) = log2(max(candidateScale.enrichment_ratio(loc), 1e-3));
end
imagesc(ax, Z);
i_apply_diverging_colormap(ax, Z);
set(ax, 'XTick', 1:numel(scaleSeconds), ...
    'XTickLabel', compose('%.2g', scaleSeconds), ...
    'YTick', 1:nCandidates, 'YTickLabel', shortIds, ...
    'TickLabelInterpreter', 'none');
xlabel(ax, 'Chunk duration (s)');
ylabel(ax, 'Frozen candidate');
title(ax, '(B) Post-fit scale enrichment, log_2(observed/expected)');
colorbar(ax);

ax = nexttile(layout);
Session = Postfit.session_composition;
entropy = zeros(nCandidates, 1);
maximumShare = zeros(nCandidates, 1);
effectiveSessions = zeros(nCandidates, 1);
for iCandidate = 1:nCandidates
    rows = string(Session.motif_candidate_id) == candidateIds(iCandidate);
    entropy(iCandidate) = Session.normalized_mixing_entropy(find(rows, 1));
    effectiveSessions(iCandidate) = ...
        Session.effective_category_count(find(rows, 1));
    maximumShare(iCandidate) = max(Session.candidate_fraction(rows));
end
plot(ax, entropy, 1:nCandidates, 'o', ...
    'Color', [0.00 0.45 0.70], 'MarkerFaceColor', [0.00 0.45 0.70], ...
    'DisplayName', 'Normalized entropy');
hold(ax, 'on');
plot(ax, maximumShare, 1:nCandidates, 's', ...
    'Color', [0.80 0.40 0.00], 'MarkerFaceColor', [0.80 0.40 0.00], ...
    'DisplayName', 'Largest-session fraction');
hold(ax, 'off');
xlabel(ax, 'Session-mixing metric');
xlim(ax, [0 1]);
set(ax, 'YTick', 1:nCandidates, 'YTickLabel', shortIds, ...
    'YDir', 'reverse');
ylabel(ax, 'Frozen candidate');
title(ax, '(C) Post-fit session mixing');
legend(ax, 'Location', 'southwest');
grid(ax, 'on');
for iCandidate = 1:nCandidates
    if effectiveSessions(iCandidate) < 50
        text(ax, maximumShare(iCandidate), iCandidate, ...
            sprintf('  n_{eff}=%.1f', effectiveSessions(iCandidate)), ...
            'FontSize', 9);
    end
end

ax = nexttile(layout);
Annotation = Postfit.annotation;
[present, annotationLoc] = ismember(candidateIds, ...
    string(Annotation.motif_candidate_id));
assert(all(present), ...
    'make_run09_candidate_paper_figures:MissingAnnotationCandidate', ...
    'Every eligible paper candidate requires an annotation row.');
scores = Annotation.social_evidence_score(annotationLoc);
barh(ax, scores, 'FaceColor', [0.55 0.55 0.55], ...
    'EdgeColor', 'none');
hold(ax, 'on');
xline(ax, Annotation.social_score_threshold(annotationLoc(1)), '--', ...
    'Color', [0.00 0.45 0.70]);
xline(ax, Annotation.non_social_score_threshold(annotationLoc(1)), '--', ...
    'Color', [0.80 0.40 0.00]);
social = Annotation.provisional_annotation(annotationLoc) == ...
    "provisional_social";
scatter(ax, scores(social), find(social), 60, ...
    [0.00 0.45 0.70], 'filled', 'MarkerEdgeColor', 'k');
hold(ax, 'off');
set(ax, 'YTick', 1:nCandidates, 'YTickLabel', shortIds, ...
    'YDir', 'reverse');
xlabel(ax, 'Predeclared social-evidence score');
ylabel(ax, 'Frozen candidate');
title(ax, '(D) Provisional post-freeze annotation');
grid(ax, 'on');

sgtitle(layout, ...
    'Organization of interpretation-eligible motif candidates', ...
    'FontWeight', 'bold');
end

function fig = i_postfit_profile_figure(Postfit)
candidateIds = Postfit.candidate_ids(Postfit.eligible);
shortIds = i_short_candidate_ids(candidateIds);
behaviorFeatures = [ ...
    "centroid_dist", "head2head_dist", "radial_speed_12", ...
    "approach_speed_1", "approach_speed_2", "mutual_facing", ...
    "in_contact", "head_close", "close_pair", ...
    "mutual_approach", "withdrawal", "asym_investigate"];
behaviorLabels = [ ...
    "Centroid distance", "Head distance", "Radial speed", ...
    "Approach mouse 1", "Approach mouse 2", "Mutual facing", ...
    "Contact", "Head close", "Close pair", ...
    "Mutual approach", "Withdrawal", "Role asymmetry"];
eventFeatures = [ ...
    "contact_dwell_fraction", "mutual_approach_dwell_fraction", ...
    "withdrawal_dwell_fraction", "centroid_distance_mean_mm", ...
    "centroid_distance_min_mm", "radial_speed_mean_mm_s", ...
    "mutual_facing_mean", "role_approach_imbalance"];
eventLabels = [ ...
    "Contact dwell", "Mutual approach dwell", "Withdrawal dwell", ...
    "Mean distance", "Minimum distance", "Radial speed", ...
    "Mutual facing", "Role imbalance"];

fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [50 50 1800 900]);
layout = tiledlayout(fig, 1, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(layout);
Z = i_profile_matrix(Postfit.behavioral_profile, ...
    candidateIds, behaviorFeatures);
imagesc(ax, Z);
i_apply_diverging_colormap(ax, Z);
i_profile_axes(ax, behaviorLabels, shortIds);
title(ax, '(A) Behavioral profile');
colorbar(ax);

ax = nexttile(layout);
Z = i_profile_matrix(Postfit.event_profile, ...
    candidateIds, eventFeatures);
imagesc(ax, Z);
i_apply_diverging_colormap(ax, Z);
i_profile_axes(ax, eventLabels, shortIds);
title(ax, '(B) Event-summary profile');
colorbar(ax);

sgtitle(layout, ...
    ['Post-freeze condition-blind characterization ' ...
     '(not independent behavioral validation)'], ...
    'FontWeight', 'bold');
end

function Z = i_profile_matrix(Profile, candidateIds, featureNames)
Z = zeros(numel(candidateIds), numel(featureNames));
for iCandidate = 1:numel(candidateIds)
    for iFeature = 1:numel(featureNames)
        row = string(Profile.motif_candidate_id) == ...
            candidateIds(iCandidate) & ...
            string(Profile.feature_name) == featureNames(iFeature);
        assert(nnz(row) == 1, ...
            'make_run09_candidate_paper_figures:MissingProfileCell', ...
            'Each selected candidate-profile cell must occur exactly once.');
        Z(iCandidate, iFeature) = ...
            Profile.standardized_mean_difference(row);
    end
end
end

function i_profile_axes(ax, featureLabels, candidateLabels)
set(ax, 'XTick', 1:numel(featureLabels), ...
    'XTickLabel', featureLabels, 'XTickLabelRotation', 40, ...
    'YTick', 1:numel(candidateLabels), ...
    'YTickLabel', candidateLabels, 'TickLabelInterpreter', 'none');
xlabel(ax, 'Condition-blind feature');
ylabel(ax, 'Frozen candidate');
end

function i_label_candidates(ax, Topology, mask, x, y)
rows = find(mask);
labels = i_short_candidate_ids( ...
    string(Topology.motif_candidate_id(rows)));
if isempty(rows)
    return;
end

xLimits = xlim(ax);
yLimits = ylim(ax);
ySpan = max(diff(yLimits), eps);
[~, order] = sort(y(rows));
orderedRows = rows(order);
orderedLabels = labels(order);
labelY = y(orderedRows);
minimumSeparation = 0.025 * ySpan;
for i = 2:numel(labelY)
    labelY(i) = max(labelY(i), labelY(i - 1) + minimumSeparation);
end
upperLimit = yLimits(2) - 0.015 * ySpan;
if labelY(end) > upperLimit
    labelY = labelY - (labelY(end) - upperLimit);
end
lowerLimit = yLimits(1) + 0.015 * ySpan;
for i = numel(labelY) - 1:-1:1
    labelY(i) = min(labelY(i), labelY(i + 1) - minimumSeparation);
end
if labelY(1) < lowerLimit
    labelY = labelY + (lowerLimit - labelY(1));
end

if string(ax.XScale) == "log"
    labelX = x(orderedRows) * 1.04;
else
    labelX = x(orderedRows) + 0.008 * diff(xLimits);
end
for i = 1:numel(orderedRows)
    plot(ax, [x(orderedRows(i)), labelX(i)], ...
        [y(orderedRows(i)), labelY(i)], '-', ...
        'Color', [0.55 0.55 0.55], 'LineWidth', 0.5, ...
        'HandleVisibility', 'off');
    text(ax, labelX(i), labelY(i), " " + orderedLabels(i), ...
        'FontSize', 8, 'Interpreter', 'none', ...
        'VerticalAlignment', 'middle');
end
end

function values = i_numeric(x)
if isnumeric(x)
    values = double(x);
else
    values = str2double(string(x));
end
end

function labels = i_short_candidate_ids(candidateIds)
labels = regexprep(string(candidateIds), '^MC_M[0-9]+_', '');
end

function i_apply_diverging_colormap(ax, values)
maxAbs = max(abs(values), [], 'all', 'omitnan');
if ~isfinite(maxAbs) || maxAbs == 0
    maxAbs = 1;
end
colormap(ax, i_diverging_map(257));
clim(ax, [-maxAbs maxAbs]);
end

function map = i_diverging_map(n)
blue = [0.00 0.45 0.70];
middle = [0.97 0.97 0.97];
orange = [0.80 0.40 0.00];
nLeft = ceil(n / 2);
nRight = n - nLeft + 1;
left = [linspace(blue(1), middle(1), nLeft)', ...
    linspace(blue(2), middle(2), nLeft)', ...
    linspace(blue(3), middle(3), nLeft)'];
right = [linspace(middle(1), orange(1), nRight)', ...
    linspace(middle(2), orange(2), nRight)', ...
    linspace(middle(3), orange(3), nRight)'];
map = [left; right(2:end, :)];
end

function records = i_export_pair( ...
        fig, root, id, role, inputs, paperRole, params)
rasterName = id + "." + params.paper_figure_raster_format;
vectorName = id + "." + params.paper_figure_vector_format;
exportgraphics(fig, fullfile(root, rasterName), ...
    'Resolution', params.paper_figure_raster_dpi);
exportgraphics(fig, fullfile(root, vectorName), ...
    'ContentType', 'vector');
close(fig);
records = [ ...
    struct('id', id + "_raster", 'file', rasterName, 'role', role, ...
    'inputs', inputs, 'format', params.paper_figure_raster_format, ...
    'paperRole', paperRole)
    struct('id', id + "_vector", 'file', vectorName, 'role', role, ...
    'inputs', inputs, 'format', params.paper_figure_vector_format, ...
    'paperRole', paperRole)];
end

function T = i_read(root, fileName)
T = readtable(fullfile(root, fileName), ...
    'Delimiter', ',', 'TextType', 'string');
end
