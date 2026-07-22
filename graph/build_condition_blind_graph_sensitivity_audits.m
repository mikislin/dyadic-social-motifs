function Audit = build_condition_blind_graph_sensitivity_audits(X, nodeManifest, primaryGraph, params, outRoot)
%BUILD_CONDITION_BLIND_GRAPH_SENSITIVITY_AUDITS Run condition-blind graph audits.
%
% Alternative graphs are sensitivity artifacts only. Numeric run_07 scores are
% the sole distance inputs. Anchor stage and session identity can restrict the
% sensitivity node/candidate universe but never modify the primary graph.

outRoot = string(outRoot);
X = double(X);
assert(all(isfinite(X), 'all'), ...
    'build_condition_blind_graph_sensitivity_audits:NonFiniteInput', ...
    'Sensitivity input scores must be finite.');
assert(istable(nodeManifest) && height(nodeManifest) == size(X, 1), ...
    'build_condition_blind_graph_sensitivity_audits:BadNodeManifest', ...
    'nodeManifest must align one-to-one with X.');

paths = local_paths(outRoot);
n = height(nodeManifest);
globalIds = double(nodeManifest.graph_node_id);
if ismember('anchor_stage', nodeManifest.Properties.VariableNames)
    stage = string(nodeManifest.anchor_stage);
else
    stage = repmat("base_time_even", n, 1);
    nodeManifest.anchor_stage = stage;
end
baseMask = stage == "base_time_even";
rareMask = stage == "rare_strata_enriched";
if ~any(baseMask)
    baseMask = ~rareMask;
end

rng(params.audit_random_seed, 'twister');
variantNames = ["combined"; "base_only"; "rare_only"; "stage_balanced"];
variantMasks = {true(n, 1); baseMask; rareMask; local_stage_balanced_mask(nodeManifest, stage)};
StageSummary = table();
StageByScale = table();
StageSelection = table();
StageNodeAudit = table();
baseGraph = struct();
baseMeta = table();

primaryReferenceKeys = local_edge_keys(primaryGraph.Edges, globalIds, n);
for v = 1:numel(variantNames)
    variant = variantNames(v);
    mask = logical(variantMasks{v});
    if nnz(mask) <= params.k_neighbors
        StageSummary = [StageSummary; local_unavailable_variant_row(variant, nnz(mask), params)]; %#ok<AGROW>
        continue
    end
    Meta = nodeManifest(mask, :);
    ids = globalIds(mask);
    if variant == "combined"
        G = primaryGraph;
    else
        G = build_condition_blind_knn_graph(X(mask, :), Meta, params);
    end
    [oneSummary, oneScale, oneNode] = local_graph_metrics( ...
        G, Meta, ids, variant, primaryReferenceKeys, n);
    StageSummary = [StageSummary; oneSummary]; %#ok<AGROW>
    StageByScale = [StageByScale; oneScale]; %#ok<AGROW>
    StageNodeAudit = [StageNodeAudit; oneNode]; %#ok<AGROW>
    StageSelection = [StageSelection; local_selection_rows( ...
        variant, mask, nodeManifest, stage, params.audit_random_seed)]; %#ok<AGROW>
    if variant == "base_only"
        baseGraph = G;
        baseMeta = Meta;
        local_write_base_graph(paths, G, Meta, ids, oneSummary, oneNode);
    end
end

writetable(StageSummary, paths.stageSummary);
writetable(StageByScale, paths.stageByScale);
writetable(StageSelection, paths.stageSelection);
writetable(StageNodeAudit, paths.stageNode);

if logical(params.audit_session_excluded_enabled)
    SessionGraph = build_condition_blind_session_excluded_knn_audit(X, nodeManifest, params);
    [SessionSummary, SessionNode] = local_session_excluded_metrics( ...
        primaryGraph, SessionGraph, nodeManifest);
    SessionEdges = local_edges_with_global_ids(SessionGraph.Edges, globalIds, ...
        "session_excluded_audit_only");
    writetable(SessionSummary, paths.sessionSummary);
    writetable(SessionNode, paths.sessionNode);
    if logical(params.audit_persist_session_excluded_edges)
        writetable(SessionEdges, paths.sessionEdges);
    end
else
    SessionSummary = local_disabled_summary("session_excluded", params);
    SessionNode = table();
    writetable(SessionSummary, paths.sessionSummary);
end

[Resampling, Panel, ResamplingSelection] = local_resampling_audit(X, nodeManifest, params);
writetable(Resampling, paths.resampling);
writetable(Panel, paths.resamplingPanel);
writetable(ResamplingSelection, paths.resamplingSelection);

Audit = struct();
Audit.stageSummary = StageSummary;
Audit.stageByScale = StageByScale;
Audit.stageSelection = StageSelection;
Audit.stageNodeAudit = StageNodeAudit;
Audit.sessionSummary = SessionSummary;
Audit.sessionNodeAudit = SessionNode;
Audit.resampling = Resampling;
Audit.resamplingPanel = Panel;
Audit.resamplingSelection = ResamplingSelection;
Audit.paths = paths;
Audit.base_graph_available = ~isempty(fieldnames(baseGraph));
Audit.base_graph_n_nodes = height(baseMeta);
Audit.labels_used_for_distance = "none";
Audit.arena_used_for_distance = false;
Audit.condition_used_for_distance = false;
end

function paths = local_paths(outRoot)
paths = struct();
paths.stageSummary = fullfile(outRoot, 'graph_anchor_stage_sensitivity_audit.csv');
paths.stageByScale = fullfile(outRoot, 'graph_anchor_stage_sensitivity_by_scale.csv');
paths.stageSelection = fullfile(outRoot, 'graph_anchor_stage_sensitivity_selection_manifest.csv');
paths.stageNode = fullfile(outRoot, 'graph_anchor_stage_sensitivity_node_audit.csv');
paths.baseNode = fullfile(outRoot, 'graph_primary_base_only_node_manifest.csv');
paths.baseEdges = fullfile(outRoot, 'graph_primary_base_only_edge_list.csv');
paths.baseTopology = fullfile(outRoot, 'graph_primary_base_only_topology_summary.csv');
paths.baseDegree = fullfile(outRoot, 'graph_primary_base_only_degree_audit.csv');
paths.sessionSummary = fullfile(outRoot, 'graph_session_excluded_sensitivity_audit.csv');
paths.sessionNode = fullfile(outRoot, 'graph_session_excluded_node_audit.csv');
paths.sessionEdges = fullfile(outRoot, 'graph_session_excluded_edge_list_audit.csv');
paths.resampling = fullfile(outRoot, 'graph_neighborhood_resampling_audit.csv');
paths.resamplingPanel = fullfile(outRoot, 'graph_neighborhood_resampling_panel_manifest.csv');
paths.resamplingSelection = fullfile(outRoot, 'graph_neighborhood_resampling_selection_manifest.csv');
end

function mask = local_stage_balanced_mask(Meta, stage)
n = height(Meta);
mask = false(n, 1);
scales = unique(double(Meta.scale_index), 'stable')';
for s = scales
    base = find(double(Meta.scale_index) == s & stage == "base_time_even");
    rare = find(double(Meta.scale_index) == s & stage == "rare_strata_enriched");
    take = min(numel(base), numel(rare));
    if take == 0
        continue
    end
    mask(base(local_random_take(numel(base), take))) = true;
    mask(rare(local_random_take(numel(rare), take))) = true;
end
end

function idx = local_random_take(n, take)
if take >= n
    idx = (1:n)';
else
    order = randperm(n, take);
    idx = sort(order(:));
end
end

function [Summary, ByScale, NodeAudit] = local_graph_metrics( ...
    G, Meta, globalIds, variant, referenceKeys, nGlobal)
Edges = G.Edges;
n = height(Meta);
source = double(Edges.source_node_id);
target = double(Edges.target_node_id);
outDegree = accumarray(source, 1, [n 1], @sum, 0);
inDegree = accumarray(target, 1, [n 1], @sum, 0);
radius = accumarray(source, double(Edges.neighbor_distance), [n 1], @max, NaN);
meanDistance = accumarray(source, double(Edges.neighbor_distance), [n 1], @mean, NaN);
mutualCount = accumarray(source, double(Edges.is_mutual_neighbor), [n 1], @sum, 0);
[sUndir, tUndir] = local_undirected_pairs(Edges);
H = graph(sUndir, tUndir, [], n);
component = conncomp(H)';
undirectedDegree = degree(H);
sameScale = double(Meta.scale_index(source)) == double(Meta.scale_index(target));
sameSession = string(Meta.session_index(source)) == string(Meta.session_index(target));
variantKeys = local_edge_keys(Edges, globalIds, nGlobal);
retainedEdge = ismember(variantKeys, referenceKeys);
retentionByNode = accumarray(source, double(retainedEdge), [n 1], @mean, NaN);

Summary = table();
Summary.graph_variant = string(variant);
Summary.graph_role = "condition_blind_anchor_stage_sensitivity_not_primary_graph";
Summary.status = "completed";
Summary.n_nodes = n;
Summary.n_base_nodes = nnz(string(Meta.anchor_stage) == "base_time_even");
Summary.n_rare_enriched_nodes = nnz(string(Meta.anchor_stage) == "rare_strata_enriched");
Summary.n_graph_dimensions = G.n_dims;
Summary.k_neighbors = G.k;
Summary.n_directed_edges = height(Edges);
Summary.n_undirected_edges = numel(sUndir);
Summary.n_components = numel(unique(component));
Summary.largest_component_fraction = max(accumarray(component, 1)) ./ n;
Summary.median_in_degree = median(inDegree, 'omitnan');
Summary.median_undirected_degree = median(undirectedDegree, 'omitnan');
Summary.median_knn_radius = median(radius, 'omitnan');
Summary.mean_same_scale_neighbor_fraction = mean(sameScale, 'omitnan');
Summary.mean_same_session_neighbor_fraction = mean(sameSession, 'omitnan');
Summary.directed_edge_retention_to_combined_graph = mean(retainedEdge, 'omitnan');
Summary.mean_node_neighbor_retention_to_combined_graph = mean(retentionByNode, 'omitnan');
Summary.node_subset_rule = local_variant_rule(variant);
Summary.rare_stratum_labels_used_for_edge_construction = false;
Summary.labels_used_for_distance = "none";
Summary.arena_used_for_distance = false;
Summary.condition_used_for_distance = false;

NodeAudit = table();
NodeAudit.graph_variant = repmat(string(variant), n, 1);
NodeAudit.sensitivity_node_id = (1:n)';
NodeAudit.graph_node_id = globalIds(:);
NodeAudit.embedding_row_id = Meta.embedding_row_id;
NodeAudit.scale_index = double(Meta.scale_index);
NodeAudit.chunk_sec = double(Meta.chunk_sec);
NodeAudit.anchor_stage = string(Meta.anchor_stage);
NodeAudit.knn_in_degree = inDegree;
NodeAudit.knn_out_degree = outDegree;
NodeAudit.undirected_degree = undirectedDegree;
NodeAudit.mutual_neighbor_fraction = mutualCount ./ max(outDegree, 1);
NodeAudit.mean_neighbor_distance = meanDistance;
NodeAudit.knn_radius = radius;
NodeAudit.component_id = component;
NodeAudit.neighbor_retention_to_combined_graph = retentionByNode;
NodeAudit.audit_role = repmat("condition_blind_anchor_stage_sensitivity", n, 1);
NodeAudit.labels_used_for_distance = repmat("none", n, 1);
NodeAudit.arena_used_for_distance = false(n, 1);
NodeAudit.condition_used_for_distance = false(n, 1);

ByScale = table();
scales = unique(double(Meta.scale_index), 'stable')';
for s = scales
    nodeMask = double(Meta.scale_index) == s;
    edgeMask = nodeMask(source);
    one = table();
    one.graph_variant = string(variant);
    one.scale_index = s;
    one.chunk_sec = double(Meta.chunk_sec(find(nodeMask, 1)));
    one.n_nodes = nnz(nodeMask);
    one.n_base_nodes = nnz(nodeMask & string(Meta.anchor_stage) == "base_time_even");
    one.n_rare_enriched_nodes = nnz(nodeMask & string(Meta.anchor_stage) == "rare_strata_enriched");
    one.mean_same_scale_neighbor_fraction = mean(sameScale(edgeMask), 'omitnan');
    one.mean_same_session_neighbor_fraction = mean(sameSession(edgeMask), 'omitnan');
    one.median_knn_radius = median(radius(nodeMask), 'omitnan');
    one.median_undirected_degree = median(undirectedDegree(nodeMask), 'omitnan');
    one.mean_neighbor_retention_to_combined_graph = mean(retentionByNode(nodeMask), 'omitnan');
    one.audit_role = "condition_blind_anchor_stage_sensitivity_by_scale";
    one.labels_used_for_distance = "none";
    one.arena_used_for_distance = false;
    one.condition_used_for_distance = false;
    ByScale = [ByScale; one]; %#ok<AGROW>
end
end

function row = local_unavailable_variant_row(variant, nNodes, params)
row = table();
row.graph_variant = string(variant);
row.graph_role = "condition_blind_anchor_stage_sensitivity_not_primary_graph";
row.status = "not_applicable_insufficient_stage_nodes";
row.n_nodes = nNodes;
row.n_base_nodes = NaN;
row.n_rare_enriched_nodes = NaN;
row.n_graph_dimensions = params.graph_n_pcs;
row.k_neighbors = params.k_neighbors;
row.n_directed_edges = NaN;
row.n_undirected_edges = NaN;
row.n_components = NaN;
row.largest_component_fraction = NaN;
row.median_in_degree = NaN;
row.median_undirected_degree = NaN;
row.median_knn_radius = NaN;
row.mean_same_scale_neighbor_fraction = NaN;
row.mean_same_session_neighbor_fraction = NaN;
row.directed_edge_retention_to_combined_graph = NaN;
row.mean_node_neighbor_retention_to_combined_graph = NaN;
row.node_subset_rule = local_variant_rule(variant);
row.rare_stratum_labels_used_for_edge_construction = false;
row.labels_used_for_distance = "none";
row.arena_used_for_distance = false;
row.condition_used_for_distance = false;
end

function rule = local_variant_rule(variant)
switch string(variant)
    case "combined"
        rule = "all_run07_rows_primary_graph_reuse";
    case "base_only"
        rule = "anchor_stage_equals_base_time_even";
    case "rare_only"
        rule = "anchor_stage_equals_rare_strata_enriched";
    case "stage_balanced"
        rule = "within_scale_equal_base_and_rare_counts_fixed_seed";
    otherwise
        rule = "unknown";
end
end

function Selection = local_selection_rows(variant, mask, Meta, stage, seed)
idx = find(mask);
Selection = table();
Selection.graph_variant = repmat(string(variant), numel(idx), 1);
Selection.graph_node_id = double(Meta.graph_node_id(idx));
Selection.embedding_row_id = Meta.embedding_row_id(idx);
Selection.scale_index = double(Meta.scale_index(idx));
Selection.chunk_sec = double(Meta.chunk_sec(idx));
Selection.anchor_stage = stage(idx);
Selection.within_scale_stage_selection_probability = ones(numel(idx), 1);
if variant == "stage_balanced"
    scales = unique(double(Meta.scale_index), 'stable')';
    for s = scales
        for st = ["base_time_even", "rare_strata_enriched"]
            pool = double(Meta.scale_index) == s & stage == st;
            selected = mask & pool;
            p = nnz(selected) ./ max(nnz(pool), 1);
            row = double(Meta.scale_index(idx)) == s & stage(idx) == st;
            Selection.within_scale_stage_selection_probability(row) = p;
        end
    end
end
Selection.selection_seed = repmat(seed, numel(idx), 1);
Selection.node_subset_rule = repmat(local_variant_rule(variant), numel(idx), 1);
Selection.rare_stratum_labels_used_for_selection = false(numel(idx), 1);
Selection.labels_used_for_distance = repmat("none", numel(idx), 1);
Selection.condition_used_for_selection = false(numel(idx), 1);
end

function local_write_base_graph(paths, G, Meta, ids, Summary, NodeAudit)
BaseNode = Meta;
BaseNode.primary_base_graph_node_id = (1:height(BaseNode))';
BaseNode.combined_graph_node_id = ids(:);
BaseNode = movevars(BaseNode, {'primary_base_graph_node_id','combined_graph_node_id'}, 'Before', 1);
BaseNode.graph_role = repmat("corrected_primary_base_only_same_run07_representation", height(BaseNode), 1);
BaseNode.labels_used_for_distance = repmat("none", height(BaseNode), 1);
BaseNode.arena_used_for_distance = false(height(BaseNode), 1);
BaseNode.condition_used_for_distance = false(height(BaseNode), 1);
BaseEdges = local_edges_with_global_ids(G.Edges, ids, "corrected_primary_base_only_sensitivity");
writetable(BaseNode, paths.baseNode);
writetable(BaseEdges, paths.baseEdges);
writetable(Summary, paths.baseTopology);
writetable(NodeAudit, paths.baseDegree);
end

function E = local_edges_with_global_ids(Edges, ids, role)
E = Edges;
E.source_combined_graph_node_id = ids(double(Edges.source_node_id));
E.target_combined_graph_node_id = ids(double(Edges.target_node_id));
E.sensitivity_graph_role = repmat(string(role), height(E), 1);
E.rare_stratum_labels_used_for_edge_construction = false(height(E), 1);
E.labels_used_for_distance = repmat("none", height(E), 1);
E.arena_used_for_distance = false(height(E), 1);
E.condition_used_for_distance = false(height(E), 1);
end

function [Summary, NodeAudit] = local_session_excluded_metrics(Primary, Session, Meta)
n = height(Meta);
k = Primary.k;
primaryTargets = local_target_matrix(Primary.Edges, n, k);
sessionTargets = local_target_matrix(Session.Edges, n, Session.k);
overlap = zeros(n, 1);
jaccard = zeros(n, 1);
for i = 1:n
    overlap(i) = numel(intersect(primaryTargets(i, :), sessionTargets(i, :)));
    jaccard(i) = overlap(i) ./ max(2 * k - overlap(i), 1);
end
primaryRadius = max(reshape(double(Primary.Edges.neighbor_distance), k, n), [], 1)';
sessionRadius = max(reshape(double(Session.Edges.neighbor_distance), Session.k, n), [], 1)';
source = double(Session.Edges.source_node_id);
target = double(Session.Edges.target_node_id);
sameScale = double(Meta.scale_index(source)) == double(Meta.scale_index(target));
sameSession = string(Meta.session_index(source)) == string(Meta.session_index(target));

NodeAudit = table();
NodeAudit.graph_node_id = double(Meta.graph_node_id);
NodeAudit.embedding_row_id = Meta.embedding_row_id;
NodeAudit.scale_index = double(Meta.scale_index);
NodeAudit.chunk_sec = double(Meta.chunk_sec);
NodeAudit.n_primary_neighbors_retained = overlap;
NodeAudit.primary_neighbor_retention_fraction = overlap ./ k;
NodeAudit.neighbor_jaccard_to_primary = jaccard;
NodeAudit.primary_knn_radius = primaryRadius;
NodeAudit.session_excluded_knn_radius = sessionRadius;
NodeAudit.session_excluded_over_primary_radius = sessionRadius ./ max(primaryRadius, eps);
NodeAudit.session_identity_used_for_exclusion = true(n, 1);
NodeAudit.audit_only_not_primary_graph = true(n, 1);
NodeAudit.labels_used_for_distance = repmat("none", n, 1);
NodeAudit.arena_used_for_distance = false(n, 1);
NodeAudit.condition_used_for_distance = false(n, 1);

Summary = table();
for s = [NaN unique(double(Meta.scale_index), 'stable')']
    if isnan(s)
        nodeMask = true(n, 1);
        scope = "all_scales";
        sec = NaN;
    else
        nodeMask = double(Meta.scale_index) == s;
        scope = "source_scale";
        sec = double(Meta.chunk_sec(find(nodeMask, 1)));
    end
    edgeMask = nodeMask(source);
    one = table();
    one.source_scope = scope;
    one.scale_index = s;
    one.chunk_sec = sec;
    one.n_nodes = nnz(nodeMask);
    one.k_neighbors = Session.k;
    one.mean_primary_neighbor_retention_fraction = mean(overlap(nodeMask) ./ k, 'omitnan');
    one.median_neighbor_jaccard_to_primary = median(jaccard(nodeMask), 'omitnan');
    one.p10_neighbor_jaccard_to_primary = quantile(jaccard(nodeMask), 0.10);
    one.mean_session_excluded_over_primary_radius = mean( ...
        sessionRadius(nodeMask) ./ max(primaryRadius(nodeMask), eps), 'omitnan');
    one.mean_same_scale_neighbor_fraction = mean(sameScale(edgeMask), 'omitnan');
    one.mean_same_session_neighbor_fraction = mean(sameSession(edgeMask), 'omitnan');
    one.session_identity_used_for_exclusion = true;
    one.audit_only_not_primary_graph = true;
    one.labels_used_for_distance = "none";
    one.arena_used_for_distance = false;
    one.condition_used_for_distance = false;
    Summary = [Summary; one]; %#ok<AGROW>
end
end

function [Audit, Panel, Selection] = local_resampling_audit(X, Meta, params)
n = height(Meta);
panelMax = min(n, floor(params.audit_resample_panel_max_nodes));
scales = unique(double(Meta.scale_index), 'stable')';
perScale = max(1, floor(panelMax ./ numel(scales)));
panelMask = false(n, 1);
selectionProb = zeros(n, 1);
for s = scales
    pool = find(double(Meta.scale_index) == s);
    take = min(numel(pool), perScale);
    chosen = pool(local_random_take(numel(pool), take));
    panelMask(chosen) = true;
    selectionProb(chosen) = take ./ numel(pool);
end
remaining = panelMax - nnz(panelMask);
if remaining > 0
    pool = find(~panelMask);
    chosen = pool(local_random_take(numel(pool), min(remaining, numel(pool))));
    panelMask(chosen) = true;
    selectionProb(chosen) = min(remaining, numel(pool)) ./ max(numel(pool), 1);
end

panelRows = find(panelMask);
PanelMeta = Meta(panelMask, :);
Xpanel = X(panelMask, :);
Panel = table();
Panel.panel_node_id = (1:numel(panelRows))';
Panel.graph_node_id = double(Meta.graph_node_id(panelRows));
Panel.embedding_row_id = Meta.embedding_row_id(panelRows);
Panel.scale_index = double(Meta.scale_index(panelRows));
Panel.chunk_sec = double(Meta.chunk_sec(panelRows));
Panel.session_index = Meta.session_index(panelRows);
Panel.panel_selection_probability = selectionProb(panelRows);
Panel.panel_selection_rule = repmat("fixed_seed_approximately_equal_count_per_scale", numel(panelRows), 1);
Panel.selection_seed = repmat(params.audit_random_seed, numel(panelRows), 1);
Panel.labels_used_for_panel_selection = repmat("none", numel(panelRows), 1);
Panel.arena_used_for_panel_selection = false(numel(panelRows), 1);
Panel.condition_used_for_panel_selection = false(numel(panelRows), 1);

Ref = build_condition_blind_knn_graph(Xpanel, PanelMeta, params);
Audit = local_resample_row("reference", 0, params.audit_random_seed, Ref, Ref, ...
    PanelMeta, size(Xpanel, 2), numel(unique(string(PanelMeta.session_index))));
Selection = table();

for r = 1:params.audit_resample_replicates
    seed = params.audit_random_seed + r;
    rng(seed, 'twister');
    nDims = max(2, round(size(Xpanel, 2) .* params.audit_resample_dimension_fraction));
    dims = sort(randperm(size(Xpanel, 2), nDims));
    G = build_condition_blind_knn_graph(Xpanel(:, dims), PanelMeta, params);
    Audit = [Audit; local_resample_row("dimension_subsample", r, seed, G, Ref, ...
        PanelMeta, nDims, numel(unique(string(PanelMeta.session_index))))]; %#ok<AGROW>
    Selection = [Selection; local_resample_selection_rows( ...
        "dimension_subsample", r, seed, "graph_pc_index", dims(:))]; %#ok<AGROW>

    rng(seed + 100000, 'twister');
    sessions = unique(string(PanelMeta.session_index), 'stable');
    nSessions = max(2, round(numel(sessions) .* params.audit_resample_session_fraction));
    selectedSessions = sessions(sort(randperm(numel(sessions), nSessions)));
    keep = ismember(string(PanelMeta.session_index), selectedSessions);
    SessionMeta = PanelMeta(keep, :);
    Gs = build_condition_blind_knn_graph(Xpanel(keep, :), SessionMeta, params);
    Audit = [Audit; local_resample_row_subset("session_subsample", r, seed + 100000, ...
        Gs, Ref, PanelMeta, keep, size(Xpanel, 2), nSessions)]; %#ok<AGROW>
    Selection = [Selection; local_resample_selection_rows( ...
        "session_subsample", r, seed + 100000, "session_index", selectedSessions(:))]; %#ok<AGROW>
end
end

function row = local_resample_row(type, replicate, seed, G, Ref, Meta, nDims, nSessions)
n = height(Meta);
k = min(G.k, Ref.k);
targets = local_target_matrix(G.Edges, n, G.k);
refTargets = local_target_matrix(Ref.Edges, n, Ref.k);
[overlapFraction, jaccard] = local_neighbor_overlap(targets(:, 1:k), refTargets(:, 1:k));
row = local_resample_summary(type, replicate, seed, G, Meta, nDims, nSessions, ...
    overlapFraction, jaccard);
end

function row = local_resample_row_subset(type, replicate, seed, G, Ref, FullMeta, keep, nDims, nSessions)
fullIndex = find(keep);
nFull = height(FullMeta);
n = numel(fullIndex);
targetsLocal = local_target_matrix(G.Edges, n, G.k);
targetsFull = fullIndex(targetsLocal);
refTargets = local_target_matrix(Ref.Edges, nFull, Ref.k);
refTargets = refTargets(fullIndex, :);
[overlapFraction, jaccard] = local_neighbor_overlap(targetsFull, refTargets);
row = local_resample_summary(type, replicate, seed, G, FullMeta(keep, :), ...
    nDims, nSessions, overlapFraction, jaccard);
end

function [overlapFraction, jaccard] = local_neighbor_overlap(targets, reference)
n = size(targets, 1);
overlapFraction = zeros(n, 1);
jaccard = zeros(n, 1);
for i = 1:n
    inter = numel(intersect(targets(i, :), reference(i, :)));
    overlapFraction(i) = inter ./ max(size(reference, 2), 1);
    jaccard(i) = inter ./ max(size(targets, 2) + size(reference, 2) - inter, 1);
end
end

function row = local_resample_summary(type, replicate, seed, G, Meta, nDims, nSessions, overlap, jaccard)
source = double(G.Edges.source_node_id);
target = double(G.Edges.target_node_id);
sameScale = double(Meta.scale_index(source)) == double(Meta.scale_index(target));
[s, t] = local_undirected_pairs(G.Edges);
H = graph(s, t, [], height(Meta));
component = conncomp(H);
row = table();
row.resampling_type = string(type);
row.replicate = replicate;
row.random_seed = seed;
row.n_nodes = height(Meta);
row.n_graph_dimensions = nDims;
row.n_sessions = nSessions;
row.k_neighbors = G.k;
row.mean_reference_neighbor_retention_fraction = mean(overlap, 'omitnan');
row.median_reference_neighbor_retention_fraction = median(overlap, 'omitnan');
row.p10_reference_neighbor_retention_fraction = quantile(overlap, 0.10);
row.mean_neighbor_jaccard_to_panel_reference = mean(jaccard, 'omitnan');
row.median_neighbor_jaccard_to_panel_reference = median(jaccard, 'omitnan');
row.p10_neighbor_jaccard_to_panel_reference = quantile(jaccard, 0.10);
row.mean_same_scale_neighbor_fraction = mean(sameScale, 'omitnan');
row.n_components = numel(unique(component));
row.largest_component_fraction = max(accumarray(component(:), 1)) ./ height(Meta);
row.audit_panel_rule = "fixed_seed_scale_balanced_condition_blind_node_panel";
row.audit_role = "neighborhood_resampling_sensitivity_not_primary_graph";
row.labels_used_for_distance = "none";
row.arena_used_for_distance = false;
row.condition_used_for_distance = false;
end

function T = local_resample_selection_rows(type, replicate, seed, unit, values)
values = string(values(:));
T = table();
T.resampling_type = repmat(string(type), numel(values), 1);
T.replicate = repmat(replicate, numel(values), 1);
T.random_seed = repmat(seed, numel(values), 1);
T.selection_unit = repmat(string(unit), numel(values), 1);
T.selected_value = values;
T.labels_used_for_distance = repmat("none", numel(values), 1);
T.arena_used_for_distance = false(numel(values), 1);
T.condition_used_for_distance = false(numel(values), 1);
end

function T = local_disabled_summary(name, params)
T = table();
T.audit_name = string(name);
T.status = "disabled_by_config";
T.random_seed = params.audit_random_seed;
T.labels_used_for_distance = "none";
T.arena_used_for_distance = false;
T.condition_used_for_distance = false;
end

function targetMatrix = local_target_matrix(Edges, n, k)
assert(height(Edges) == n * k, ...
    'build_condition_blind_graph_sensitivity_audits:BadEdgeShape', ...
    'Expected exactly n*k directed edges ordered by source and rank.');
targetMatrix = reshape(double(Edges.target_node_id), k, n)';
end

function keys = local_edge_keys(Edges, ids, nGlobal)
sourceGlobal = uint64(ids(double(Edges.source_node_id)) - 1);
targetGlobal = uint64(ids(double(Edges.target_node_id)));
keys = sourceGlobal .* uint64(nGlobal) + targetGlobal;
end

function [sUndir, tUndir] = local_undirected_pairs(Edges)
s = min(double(Edges.source_node_id), double(Edges.target_node_id));
t = max(double(Edges.source_node_id), double(Edges.target_node_id));
pairs = unique([s t], 'rows');
sUndir = pairs(:, 1);
tUndir = pairs(:, 2);
end
