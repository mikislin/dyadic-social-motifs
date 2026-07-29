function Audit = build_condition_blind_consensus_neighborhood( ...
    X, nodeManifest, primaryGraph, eventNodeAudit, params, outRoot)
%BUILD_CONDITION_BLIND_CONSENSUS_NEIGHBORHOOD Stabilize run-08 neighborhoods.
%
% This run-08 extension repeatedly builds exact numeric-score kNN graphs on
% within-scale anchor-stage-balanced node sets and predefined PC subsets. It
% aggregates only edge recurrence, co-inclusion, rank, and distance. No
% condition, cohort, arena, event, or rare-stratum identity enters a distance
% or consensus edge weight. Event channels are joined after the frozen graph
% is selected for coverage auditing only.

X = double(X);
outRoot = string(outRoot);
assert(all(isfinite(X), 'all') && height(nodeManifest) == size(X, 1), ...
    'build_condition_blind_consensus_neighborhood:BadInput', ...
    'X must be finite and align one-to-one with nodeManifest.');
assert(size(X, 2) == params.graph_n_pcs, ...
    'build_condition_blind_consensus_neighborhood:UnexpectedDimensionCount', ...
    'Consensus must use the same locked %d-PC graph representation.', params.graph_n_pcs);

paths = local_paths(outRoot);
n = height(nodeManifest);
d = size(X, 2);
R = floor(params.consensus_replicates);
kValues = unique(floor(params.consensus_k_values(:)'), 'stable');
thresholds = unique(double(params.consensus_support_thresholds(:)'), 'stable');
primaryK = params.k_neighbors;
primaryKIndex = find(kValues == primaryK, 1);
assert(~isempty(primaryKIndex), ...
    'build_condition_blind_consensus_neighborhood:MissingPrimaryK', ...
    'consensus_k_values must include k_neighbors.');

if ismember('anchor_stage', nodeManifest.Properties.VariableNames)
    stage = string(nodeManifest.anchor_stage);
else
    stage = repmat("base_time_even", n, 1);
end
scaleIndex = double(nodeManifest.scale_index);
sessionIndex = string(nodeManifest.session_index);

eligible = false(n, R);
dimensionRows = table();
selectionRows = cell(R, 1);
replicateRows = table();
keyCells = cell(R, 1);
anyCells = cell(R, 1);
mutualCells = cell(R, 1);
loToHiCells = cell(R, 1);
hiToLoCells = cell(R, 1);
rankCells = cell(R, 1);
distanceCells = cell(R, 1);
designCells = cell(R, 1);

fprintf('Building %d condition-blind consensus graph replicates...\n', R);
for r = 1:R
    seed = params.audit_random_seed + 200000 + r;
    [dims, designName, designCode, dimRow] = local_dimension_design(d, r, seed, params);
    stageOrderSeed = params.audit_random_seed + 300000;
    [keep, selectionRow] = local_stage_balanced_selection( ...
        nodeManifest, stage, scaleIndex, r, seed, stageOrderSeed, ...
        designCode, designName);
    eligible(:, r) = keep;
    dimensionRows = [dimensionRows; dimRow]; %#ok<AGROW>
    selectionRows{r} = selectionRow;

    buildParams = params;
    buildParams.k_neighbors = min(params.consensus_candidate_k, nnz(keep) - 1);
    ticReplicate = tic;
    G = build_condition_blind_knn_graph(X(keep, dims), nodeManifest(keep, :), buildParams);
    elapsed = toc(ticReplicate);
    ids = double(nodeManifest.graph_node_id(keep));
    P = local_replicate_pair_records(G.Edges, ids, n, kValues, primaryK);
    keyCells{r} = P.key;
    anyCells{r} = P.any_by_k;
    mutualCells{r} = P.mutual_by_k;
    loToHiCells{r} = P.lo_to_hi_primary;
    hiToLoCells{r} = P.hi_to_lo_primary;
    rankCells{r} = P.best_rank;
    distanceCells{r} = P.best_distance;
    designCells{r} = repmat(uint8(designCode), numel(P.key), 1);

    source = double(G.Edges.source_node_id);
    target = double(G.Edges.target_node_id);
    one = table();
    one.replicate_id = r;
    one.random_seed = seed;
    one.dimension_design = string(designName);
    one.dimension_design_code = designCode;
    one.n_graph_dimensions = numel(dims);
    one.selected_graph_pc_indices = strjoin(string(dims), ";");
    one.n_nodes = nnz(keep);
    one.n_base_nodes = nnz(keep & stage ~= "rare_strata_enriched");
    one.n_rare_enriched_nodes = nnz(keep & stage == "rare_strata_enriched");
    one.candidate_k = G.k;
    one.n_directed_candidate_edges = height(G.Edges);
    one.n_undirected_candidate_pairs = numel(P.key);
    one.mean_same_scale_candidate_fraction = mean( ...
        scaleIndex(ids(source)) == scaleIndex(ids(target)), 'omitnan');
    one.mean_same_session_candidate_fraction = mean( ...
        sessionIndex(ids(source)) == sessionIndex(ids(target)), 'omitnan');
    one.elapsed_seconds = elapsed;
    one.node_selection_rule = ...
        "all_base_plus_within_scale_dimension_specific_partitioned_rare_stage";
    one.dimension_selection_rule = "four_design_cycle_inside_locked_graph_pc_set";
    one.labels_used_for_distance = "none";
    one.arena_used_for_distance = false;
    one.condition_used_for_distance = false;
    one.event_channels_used_for_distance = false;
    one.rare_stratum_labels_used_for_distance = false;
    replicateRows = [replicateRows; one]; %#ok<AGROW>
    fprintf('  replicate %d/%d: %s, %d dims, %d nodes, %.1f s\n', ...
        r, R, designName, numel(dims), nnz(keep), elapsed);
end

stageSelection = vertcat(selectionRows{:});
writetable(dimensionRows, paths.dimensionManifest);
writetable(stageSelection, paths.stageManifest);
writetable(replicateRows, paths.replicateManifest);

fprintf('Aggregating sparse co-inclusion-aware edge support...\n');
allKey = vertcat(keyCells{:});
allAny = vertcat(anyCells{:});
allMutual = vertcat(mutualCells{:});
allLoToHi = vertcat(loToHiCells{:});
allHiToLo = vertcat(hiToLoCells{:});
allRank = vertcat(rankCells{:});
allDistance = vertcat(distanceCells{:});
allDesign = vertcat(designCells{:});
clear keyCells anyCells mutualCells loToHiCells hiToLoCells rankCells distanceCells designCells

[allKey, order] = sort(allKey);
allAny = allAny(order, :);
allMutual = allMutual(order, :);
allLoToHi = allLoToHi(order);
allHiToLo = allHiToLo(order);
allRank = allRank(order);
allDistance = allDistance(order);
allDesign = allDesign(order);
group = cumsum([1; double(diff(allKey) ~= 0)]);
nUnion = group(end);
unionKey = allKey([true; diff(allKey) ~= 0]);
[lo, hi] = local_decode_pair_keys(unionKey, n);
coInclusion = local_co_inclusion_counts(eligible, lo, hi);

recurrence = zeros(nUnion, numel(kValues));
mutualRecurrence = zeros(nUnion, numel(kValues));
for j = 1:numel(kValues)
    recurrence(:, j) = accumarray(group, double(allAny(:, j)), [nUnion 1], @sum, 0);
    mutualRecurrence(:, j) = accumarray(group, double(allMutual(:, j)), [nUnion 1], @sum, 0);
end
loToHiRecurrence = accumarray(group, double(allLoToHi), [nUnion 1], @sum, 0);
hiToLoRecurrence = accumarray(group, double(allHiToLo), [nUnion 1], @sum, 0);
medianRank = accumarray(group, double(allRank), [nUnion 1], @median, NaN);
medianDistance = accumarray(group, double(allDistance), [nUnion 1], @median, NaN);
minDistance = accumarray(group, double(allDistance), [nUnion 1], @min, NaN);
designRecurrence = zeros(nUnion, 4);
for code = 1:4
    designRecurrence(:, code) = accumarray(group, ...
        double(allAny(:, primaryKIndex) & allDesign == code), [nUnion 1], @sum, 0);
end
clear allKey allAny allMutual allLoToHi allHiToLo allRank allDistance allDesign group order

EdgeSupport = table();
EdgeSupport.source_node_id = lo;
EdgeSupport.target_node_id = hi;
EdgeSupport.co_inclusion_replicates = coInclusion;
EdgeSupport.eligible_for_support = coInclusion >= params.consensus_min_co_inclusion_replicates;
for j = 1:numel(kValues)
    k = kValues(j);
    EdgeSupport.(sprintf('neighbor_recurrence_k%d', k)) = recurrence(:, j);
    EdgeSupport.(sprintf('mutual_recurrence_k%d', k)) = mutualRecurrence(:, j);
    EdgeSupport.(sprintf('conditional_neighbor_support_k%d', k)) = ...
        recurrence(:, j) ./ max(coInclusion, 1);
    EdgeSupport.(sprintf('conditional_mutual_support_k%d', k)) = ...
        mutualRecurrence(:, j) ./ max(coInclusion, 1);
end
EdgeSupport.source_to_target_recurrence_primary_k = loToHiRecurrence;
EdgeSupport.target_to_source_recurrence_primary_k = hiToLoRecurrence;
EdgeSupport.source_to_target_conditional_support_primary_k = loToHiRecurrence ./ max(coInclusion, 1);
EdgeSupport.target_to_source_conditional_support_primary_k = hiToLoRecurrence ./ max(coInclusion, 1);
EdgeSupport.full_dimension_recurrence_primary_k = designRecurrence(:, 1);
EdgeSupport.medium_prefix_recurrence_primary_k = designRecurrence(:, 2);
EdgeSupport.short_prefix_recurrence_primary_k = designRecurrence(:, 3);
EdgeSupport.core_tail_recurrence_primary_k = designRecurrence(:, 4);
EdgeSupport.median_best_neighbor_rank = medianRank;
EdgeSupport.median_min_neighbor_distance = medianDistance;
EdgeSupport.minimum_neighbor_distance = minDistance;
EdgeSupport.support_denominator_rule = repmat( ...
    "replicates_where_both_endpoint_nodes_were_stage_eligible", nUnion, 1);
EdgeSupport.edge_support_role = repmat( ...
    "condition_blind_numeric_neighborhood_recurrence", nUnion, 1);
EdgeSupport.labels_used_for_edge_support = repmat("none", nUnion, 1);
EdgeSupport.arena_used_for_edge_support = false(nUnion, 1);
EdgeSupport.condition_used_for_edge_support = false(nUnion, 1);
EdgeSupport.event_channels_used_for_edge_support = false(nUnion, 1);
EdgeSupport.rare_stratum_labels_used_for_edge_support = false(nUnion, 1);
if logical(params.consensus_write_edge_support)
    writetable(EdgeSupport, paths.edgeSupport);
end

primarySource = double(primaryGraph.Edges.source_node_id);
primaryTarget = double(primaryGraph.Edges.target_node_id);
primarySameScale = mean(scaleIndex(primarySource) == scaleIndex(primaryTarget), 'omitnan');
expectedSameSession = local_same_category_null(sessionIndex);

Sensitivity = table();
metricCache = cell(numel(kValues), numel(thresholds));
for j = 1:numel(kValues)
    support = recurrence(:, j) ./ max(coInclusion, 1);
    for t = 1:numel(thresholds)
        mask = coInclusion >= params.consensus_min_co_inclusion_replicates & ...
            support >= thresholds(t);
        M = local_consensus_metrics(lo, hi, support, mask, n, stage, scaleIndex, ...
            sessionIndex, expectedSameSession, primarySameScale, params, kValues(j), thresholds(t));
        metricCache{j, t} = M;
        Sensitivity = [Sensitivity; M.summary]; %#ok<AGROW>
    end
end

primaryRows = Sensitivity.k_neighbors == primaryK;
passing = primaryRows & Sensitivity.selection_topology_gates_pass;
nPassing = nnz(passing);
if any(passing)
    candidateRows = find(passing);
    [~, pick] = max(Sensitivity.support_threshold(candidateRows));
    selectedRow = candidateRows(pick);
else
    candidateRows = find(primaryRows);
    [~, pick] = min(Sensitivity.support_threshold(candidateRows));
    selectedRow = candidateRows(pick);
end
handoffReady = Sensitivity.all_handoff_gates_pass(selectedRow) && ...
    nPassing >= params.consensus_min_passing_thresholds;
Sensitivity.selected_for_frozen_handoff = false(height(Sensitivity), 1);
Sensitivity.selected_for_frozen_handoff(selectedRow) = true;
Sensitivity.n_primary_k_passing_thresholds = repmat(nPassing, height(Sensitivity), 1);
Sensitivity.n_primary_k_topology_passing_thresholds = repmat(nPassing, height(Sensitivity), 1);
Sensitivity.minimum_passing_thresholds_required = repmat( ...
    params.consensus_min_passing_thresholds, height(Sensitivity), 1);
Sensitivity.handoff_ready = repmat(handoffReady, height(Sensitivity), 1);
writetable(Sensitivity, paths.thresholdSensitivity);

selectedK = Sensitivity.k_neighbors(selectedRow);
selectedThreshold = Sensitivity.support_threshold(selectedRow);
j = find(kValues == selectedK, 1);
t = find(abs(thresholds - selectedThreshold) < 1e-12, 1);
Selected = metricCache{j, t};
selectedSupport = recurrence(:, j) ./ max(coInclusion, 1);
selectedMask = coInclusion >= params.consensus_min_co_inclusion_replicates & ...
    selectedSupport >= selectedThreshold;

ConsensusEdges = EdgeSupport(selectedMask, :);
ConsensusEdges.consensus_edge_weight = selectedSupport(selectedMask);
ConsensusEdges.selected_k_neighbors = repmat(selectedK, height(ConsensusEdges), 1);
ConsensusEdges.selected_support_threshold = repmat(selectedThreshold, height(ConsensusEdges), 1);
ConsensusEdges.consensus_edge_rule = repmat( ...
    "undirected_pair_conditional_any_direction_recurrence", height(ConsensusEdges), 1);
ConsensusEdges.labels_used_for_consensus_edge = repmat("none", height(ConsensusEdges), 1);
ConsensusEdges.arena_used_for_consensus_edge = false(height(ConsensusEdges), 1);
ConsensusEdges.condition_used_for_consensus_edge = false(height(ConsensusEdges), 1);
ConsensusEdges.event_channels_used_for_consensus_edge = false(height(ConsensusEdges), 1);
ConsensusEdges.rare_stratum_labels_used_for_consensus_edge = false(height(ConsensusEdges), 1);
writetable(ConsensusEdges, paths.consensusEdges);

NodeStability = local_node_stability(nodeManifest, stage, eligible, Selected, ...
    selectedK, selectedThreshold, params, handoffReady);
ConsensusNodes = nodeManifest;
ConsensusNodes.consensus_degree = NodeStability.consensus_degree;
ConsensusNodes.consensus_weighted_degree = NodeStability.consensus_weighted_degree;
ConsensusNodes.consensus_component_id = NodeStability.consensus_component_id;
ConsensusNodes.consensus_in_largest_component = NodeStability.consensus_in_largest_component;
ConsensusNodes.consensus_stable_for_run09 = NodeStability.consensus_stable_for_run09;
ConsensusNodes.unstable_for_run09 = NodeStability.unstable_for_run09;
ConsensusNodes.consensus_eligible_replicates = NodeStability.consensus_eligible_replicates;
ConsensusNodes.consensus_selected_k = repmat(selectedK, n, 1);
ConsensusNodes.consensus_selected_support_threshold = repmat(selectedThreshold, n, 1);
ConsensusNodes.consensus_handoff_ready = repmat(handoffReady, n, 1);
ConsensusNodes.labels_used_for_consensus = repmat("none", n, 1);
ConsensusNodes.arena_used_for_consensus = false(n, 1);
ConsensusNodes.condition_used_for_consensus = false(n, 1);
ConsensusNodes.event_channels_used_for_consensus = false(n, 1);
writetable(ConsensusNodes, paths.consensusNodes);
writetable(NodeStability, paths.nodeStability);

stableForRun09 = logical(NodeStability.consensus_stable_for_run09);
Run09Nodes = NodeStability(stableForRun09, {'graph_node_id','embedding_row_id', ...
    'scale_index','chunk_sec','consensus_eligible_replicates','consensus_degree', ...
    'consensus_stable_induced_degree','consensus_weighted_degree', ...
    'selected_k_neighbors','selected_support_threshold'});
Run09Nodes.run09_node_input_role = repmat( ...
    "frozen_condition_blind_stable_consensus_node", height(Run09Nodes), 1);
Run09Nodes.labels_available_to_run09_graph_input = repmat("none", height(Run09Nodes), 1);
Run09Nodes.arena_available_to_run09_graph_input = false(height(Run09Nodes), 1);
Run09Nodes.condition_available_to_run09_graph_input = false(height(Run09Nodes), 1);
Run09Nodes.anchor_stage_available_to_run09_graph_input = false(height(Run09Nodes), 1);
Run09Nodes.event_channels_available_to_run09_graph_input = false(height(Run09Nodes), 1);
Run09Nodes.rare_strata_available_to_run09_graph_input = false(height(Run09Nodes), 1);
writetable(Run09Nodes, paths.run09Nodes);

stableEdge = stableForRun09(ConsensusEdges.source_node_id) & ...
    stableForRun09(ConsensusEdges.target_node_id);
Run09Edges = ConsensusEdges(stableEdge, {'source_node_id','target_node_id', ...
    'co_inclusion_replicates','median_best_neighbor_rank', ...
    'median_min_neighbor_distance','consensus_edge_weight', ...
    'selected_k_neighbors','selected_support_threshold'});
selectedValues = ConsensusEdges.(sprintf('neighbor_recurrence_k%d', selectedK));
Run09Edges.neighbor_recurrence_selected_k = selectedValues(stableEdge);
selectedValues = ConsensusEdges.(sprintf('conditional_neighbor_support_k%d', selectedK));
Run09Edges.conditional_neighbor_support_selected_k = selectedValues(stableEdge);
selectedValues = ConsensusEdges.(sprintf('mutual_recurrence_k%d', selectedK));
Run09Edges.mutual_recurrence_selected_k = selectedValues(stableEdge);
selectedValues = ConsensusEdges.(sprintf('conditional_mutual_support_k%d', selectedK));
Run09Edges.conditional_mutual_support_selected_k = selectedValues(stableEdge);
Run09Edges.run09_edge_input_role = repmat( ...
    "frozen_condition_blind_stable_consensus_edge", height(Run09Edges), 1);
Run09Edges.labels_available_to_run09_graph_input = repmat("none", height(Run09Edges), 1);
Run09Edges.arena_available_to_run09_graph_input = false(height(Run09Edges), 1);
Run09Edges.condition_available_to_run09_graph_input = false(height(Run09Edges), 1);
Run09Edges.anchor_stage_available_to_run09_graph_input = false(height(Run09Edges), 1);
Run09Edges.event_channels_available_to_run09_graph_input = false(height(Run09Edges), 1);
Run09Edges.rare_strata_available_to_run09_graph_input = false(height(Run09Edges), 1);
writetable(Run09Edges, paths.run09Edges);

Topology = Sensitivity(selectedRow, :);
Topology.consensus_graph_role = "frozen_run08_condition_blind_handoff_substrate";
Topology.handoff_status = string(local_handoff_status(handoffReady));
Topology.motifs_defined_in_run08 = false;
Topology.communities_defined_in_run08 = false;
Topology.n_run09_stable_nodes = height(Run09Nodes);
Topology.n_run09_stable_edges = height(Run09Edges);
if any(stableForRun09)
    Topology.minimum_run09_stable_induced_degree = min( ...
        NodeStability.consensus_stable_induced_degree(stableForRun09));
else
    Topology.minimum_run09_stable_induced_degree = NaN;
end
Topology.run09_stable_subgraph_rule = ...
    "largest_iterative_minimum_degree_core_of_selected_consensus_graph";
Topology.labels_used_for_consensus = "none";
Topology.arena_used_for_consensus = false;
Topology.condition_used_for_consensus = false;
writetable(Topology, paths.topology);

ScaleMixing = local_scale_mixing(Selected.lo, Selected.hi, scaleIndex, ...
    double(nodeManifest.chunk_sec));
writetable(ScaleMixing, paths.scaleMixing);
SessionSensitivity = local_session_sensitivity(Selected.lo, Selected.hi, ...
    Selected.weight, n, sessionIndex, params);
writetable(SessionSensitivity, paths.sessionSensitivity);
EventCoverage = local_event_coverage(eventNodeAudit, NodeStability, nodeManifest);
writetable(EventCoverage, paths.eventCoverage);

FreezeConfig = local_freeze_config(params, n, d, selectedK, selectedThreshold, ...
    nPassing, handoffReady);
writetable(FreezeConfig, paths.freezeConfig);
figureManifest = make_run08_consensus_qc_figures(outRoot, params);
writetable(figureManifest, paths.figureManifest);

artifactPaths = [paths.run09Edges; paths.run09Nodes; paths.freezeConfig; paths.topology; ...
    paths.consensusEdges; paths.consensusNodes; paths.nodeStability; ...
    paths.thresholdSensitivity; paths.edgeSupport; ...
    paths.dimensionManifest; paths.stageManifest; paths.replicateManifest; ...
    paths.scaleMixing; paths.sessionSensitivity; paths.eventCoverage; ...
    paths.figureManifest; fullfile(outRoot, 'graph_edge_list.csv'); ...
    fullfile(outRoot, 'graph_node_manifest.csv'); ...
    fullfile(outRoot, 'graph_event_node_audit.csv'); ...
    fullfile(outRoot, 'graph_score_preprocess_audit.csv'); ...
    fullfile(outRoot, 'graph_input_manifest.csv'); string(params.config_path); ...
    string(mfilename('fullpath')) + ".m"; ...
    fullfile(fileparts(mfilename('fullpath')), 'build_condition_blind_knn_graph.m'); ...
    fullfile(fileparts(mfilename('fullpath')), 'make_run08_consensus_qc_figures.m'); ...
    fullfile(fileparts(mfilename('fullpath')), 'load_multiscale_graph_config.m'); ...
    fullfile(fileparts(mfilename('fullpath')), 'refresh_run08_consensus_handoff_manifest.m'); ...
    fullfile(fileparts(mfilename('fullpath')), 'compute_file_sha256.m')];
artifactRows = [height(Run09Edges); height(Run09Nodes); height(FreezeConfig); height(Topology); ...
    height(ConsensusEdges); height(ConsensusNodes); height(NodeStability); ...
    height(Sensitivity); height(EdgeSupport); ...
    height(dimensionRows); height(stageSelection); height(replicateRows); ...
    height(ScaleMixing); height(SessionSensitivity); height(EventCoverage); ...
    height(figureManifest); height(primaryGraph.Edges); n; height(eventNodeAudit); d; 1; ...
    height(params.config_table); NaN; NaN; NaN; NaN; NaN; NaN];
required = false(numel(artifactPaths), 1);
required(1:4) = true;
Handoff = local_handoff_manifest(artifactPaths, artifactRows, required, ...
    outRoot, handoffReady, selectedK, selectedThreshold);
writetable(Handoff, paths.handoffManifest);

Audit = struct();
Audit.thresholdSensitivity = Sensitivity;
Audit.topologySummary = Topology;
Audit.nodeStability = NodeStability;
Audit.scaleMixing = ScaleMixing;
Audit.sessionSensitivity = SessionSensitivity;
Audit.eventCoverage = EventCoverage;
Audit.replicateManifest = replicateRows;
Audit.handoffManifest = Handoff;
Audit.handoffReady = handoffReady;
Audit.selectedK = selectedK;
Audit.selectedThreshold = selectedThreshold;
Audit.nConsensusEdges = height(ConsensusEdges);
Audit.paths = paths;
Audit.labels_used_for_consensus = "none";
Audit.arena_used_for_consensus = false;
Audit.condition_used_for_consensus = false;

fprintf('Consensus selected k=%d, support>=%.3f, edges=%d, handoff=%s\n', ...
    selectedK, selectedThreshold, height(ConsensusEdges), local_handoff_status(handoffReady));
end

function paths = local_paths(outRoot)
paths = struct();
paths.dimensionManifest = fullfile(outRoot, 'graph_dimension_resample_manifest.csv');
paths.stageManifest = fullfile(outRoot, 'graph_stage_balanced_resample_manifest.csv');
paths.replicateManifest = fullfile(outRoot, 'graph_replicate_manifest.csv');
paths.edgeSupport = fullfile(outRoot, 'graph_resampled_edge_support_audit.csv');
paths.consensusEdges = fullfile(outRoot, 'graph_consensus_edge_list.csv');
paths.consensusNodes = fullfile(outRoot, 'graph_consensus_node_manifest.csv');
paths.nodeStability = fullfile(outRoot, 'graph_consensus_node_stability_audit.csv');
paths.thresholdSensitivity = fullfile(outRoot, 'graph_consensus_threshold_sensitivity.csv');
paths.topology = fullfile(outRoot, 'graph_consensus_topology_summary.csv');
paths.scaleMixing = fullfile(outRoot, 'graph_consensus_scale_mixing_audit.csv');
paths.sessionSensitivity = fullfile(outRoot, 'graph_consensus_session_sensitivity_audit.csv');
paths.eventCoverage = fullfile(outRoot, 'graph_consensus_rare_event_coverage_audit.csv');
paths.freezeConfig = fullfile(outRoot, 'graph_consensus_freeze_config.csv');
paths.figureManifest = fullfile(outRoot, 'graph_consensus_qc_figure_manifest.csv');
paths.handoffManifest = fullfile(outRoot, 'run08_to_run09_handoff_manifest.csv');
paths.run09Nodes = fullfile(outRoot, 'run08_to_run09_node_input.csv');
paths.run09Edges = fullfile(outRoot, 'run08_to_run09_edge_list.csv');
end

function [dims, name, code, Rows] = local_dimension_design(d, replicate, seed, params)
design = mod(replicate - 1, 4) + 1;
switch design
    case 1
        dims = 1:d;
        name = "full_locked_pc_set";
    case 2
        dims = 1:min(d, floor(params.consensus_prefix_medium_n_pcs));
        name = "medium_pc_prefix";
    case 3
        dims = 1:min(d, floor(params.consensus_prefix_short_n_pcs));
        name = "short_pc_prefix";
    otherwise
        core = min(floor(params.consensus_core_n_pcs), d - 1);
        nTake = min(d, max(core + 1, round(d .* params.consensus_dimension_fraction)));
        tail = (core + 1):d;
        rng(seed, 'twister');
        tailTake = sort(tail(randperm(numel(tail), nTake - core)));
        dims = [1:core tailTake];
        name = "core_plus_tail_draw";
end
dims = double(dims(:)');
Rows = table();
Rows.replicate_id = repmat(replicate, numel(dims), 1);
Rows.random_seed = repmat(seed, numel(dims), 1);
Rows.dimension_design = repmat(string(name), numel(dims), 1);
Rows.dimension_design_code = repmat(design, numel(dims), 1);
Rows.graph_pc_index = dims(:);
Rows.is_core_pc = dims(:) <= min(params.consensus_core_n_pcs, d);
Rows.selected_for_replicate = true(numel(dims), 1);
Rows.dimension_source = repmat("locked_run07_global_pc_graph_matrix", numel(dims), 1);
Rows.labels_used_for_dimension_selection = repmat("none", numel(dims), 1);
Rows.arena_used_for_dimension_selection = false(numel(dims), 1);
Rows.condition_used_for_dimension_selection = false(numel(dims), 1);
code = design;
end

function [keep, Rows] = local_stage_balanced_selection( ...
    Meta, stage, scales, replicate, seed, stageOrderSeed, designCode, designName)
n = height(Meta);
keep = stage ~= "rare_strata_enriched";
selectionProbability = ones(n, 1);
rareMask = stage == "rare_strata_enriched";
if any(rareMask)
    uniqueScales = unique(scales, 'stable')';
    for s = uniqueScales
        base = find(scales == s & ~rareMask);
        rare = find(scales == s & rareMask);
        if isempty(rare)
            continue
        end
        take = min(numel(base), numel(rare));
        if take == 0
            take = min(numel(rare), max(1, ceil(numel(rare) ./ 3)));
        end
        % Each dimension design has its own fixed within-scale partition.
        % Across the three occurrences of that design in the 12-replicate
        % production cycle, every current rare-stage node is exposed once.
        % Independent partitions across designs broaden rare-rare endpoint
        % co-inclusion without replacing exact exposure by expected exposure.
        rng(stageOrderSeed + 100000 .* designCode + round(1000 .* s), 'twister');
        order = rare(randperm(numel(rare)));
        designOccurrence = floor((replicate - 1) ./ 4);
        offset = mod(designOccurrence .* take, numel(rare));
        positions = mod(offset + (0:(take - 1)), numel(rare)) + 1;
        chosen = order(positions);
        keep(chosen) = true;
        selectionProbability(rare) = take ./ numel(rare);
    end
end
idx = find(keep);
Rows = table();
Rows.replicate_id = repmat(replicate, numel(idx), 1);
Rows.random_seed = repmat(seed, numel(idx), 1);
Rows.stage_order_seed = repmat(stageOrderSeed, numel(idx), 1);
Rows.dimension_design = repmat(string(designName), numel(idx), 1);
Rows.dimension_design_code = repmat(designCode, numel(idx), 1);
Rows.graph_node_id = double(Meta.graph_node_id(idx));
Rows.embedding_row_id = Meta.embedding_row_id(idx);
Rows.scale_index = scales(idx);
Rows.chunk_sec = double(Meta.chunk_sec(idx));
Rows.anchor_stage = stage(idx);
Rows.replicate_inclusion_probability = selectionProbability(idx);
Rows.selection_rule = repmat( ...
    "all_base_plus_within_scale_dimension_specific_partitioned_rare_stage", numel(idx), 1);
Rows.rare_stratum_labels_used_for_selection = false(numel(idx), 1);
Rows.labels_used_for_selection = repmat("none", numel(idx), 1);
Rows.arena_used_for_selection = false(numel(idx), 1);
Rows.condition_used_for_selection = false(numel(idx), 1);
end

function P = local_replicate_pair_records(Edges, ids, nGlobal, kValues, primaryK)
source = ids(double(Edges.source_node_id));
target = ids(double(Edges.target_node_id));
rank = double(Edges.neighbor_rank);
distance = double(Edges.neighbor_distance);
lo = min(source, target);
hi = max(source, target);
key = uint64(lo - 1) .* uint64(nGlobal) + uint64(hi);
[uniqueKey, ~, group] = unique(key, 'sorted');
nPair = numel(uniqueKey);
loToHi = source == lo;
anyByK = false(nPair, numel(kValues));
mutualByK = false(nPair, numel(kValues));
for j = 1:numel(kValues)
    within = rank <= kValues(j);
    directionA = accumarray(group, double(within & loToHi), [nPair 1], @max, 0) > 0;
    directionB = accumarray(group, double(within & ~loToHi), [nPair 1], @max, 0) > 0;
    anyByK(:, j) = directionA | directionB;
    mutualByK(:, j) = directionA & directionB;
end
withinPrimary = rank <= primaryK;
P = struct();
P.key = uniqueKey;
P.any_by_k = anyByK;
P.mutual_by_k = mutualByK;
P.lo_to_hi_primary = accumarray(group, double(withinPrimary & loToHi), ...
    [nPair 1], @max, 0) > 0;
P.hi_to_lo_primary = accumarray(group, double(withinPrimary & ~loToHi), ...
    [nPair 1], @max, 0) > 0;
P.best_rank = uint16(accumarray(group, rank, [nPair 1], @min, 0));
P.best_distance = single(accumarray(group, distance, [nPair 1], @min, NaN));
end

function [lo, hi] = local_decode_pair_keys(key, n)
keyDouble = double(key);
lo = floor((keyDouble - 1) ./ n) + 1;
hi = keyDouble - (lo - 1) .* n;
lo = double(lo(:));
hi = double(hi(:));
end

function counts = local_co_inclusion_counts(eligible, lo, hi)
counts = zeros(numel(lo), 1);
block = 250000;
for first = 1:block:numel(lo)
    last = min(numel(lo), first + block - 1);
    counts(first:last) = sum(eligible(lo(first:last), :) & eligible(hi(first:last), :), 2);
end
end

function M = local_consensus_metrics(loAll, hiAll, supportAll, mask, n, stage, ...
    scales, sessions, expectedSameSession, primarySameScale, params, k, threshold)
lo = loAll(mask);
hi = hiAll(mask);
weight = supportAll(mask);
if isempty(lo)
    H = graph([], [], [], n);
else
    H = graph(lo, hi, weight, n);
end
deg = degree(H);
weightedDegree = full(sum(adjacency(H, 'weighted'), 2));
component = conncomp(H)';
componentCounts = accumarray(component, 1);
[largestCount, largestId] = max(componentCounts);
sameScale = scales(lo) == scales(hi);
sameSession = sessions(lo) == sessions(hi);
[stable, stableInducedDegree] = local_stable_core(H, params.consensus_stable_min_degree);
rare = stage == "rare_strata_enriched";
if any(rare)
    rareStableFraction = mean(deg(rare) >= params.consensus_rare_min_degree);
else
    rareStableFraction = 1;
end
if isempty(lo)
    sameScaleFraction = 0;
    sameSessionFraction = 0;
else
    sameScaleFraction = mean(sameScale);
    sameSessionFraction = mean(sameSession);
end
sessionRatio = sameSessionFraction ./ max(expectedSameSession, eps);
gateStable = mean(stable) >= params.consensus_gate_min_stable_node_fraction;
gateComponent = largestCount ./ n >= params.consensus_gate_min_largest_component_fraction;
gateRare = rareStableFraction >= params.consensus_gate_min_rare_stable_fraction;
gateSession = sessionRatio <= params.consensus_gate_max_same_session_null_ratio;
gateScale = sameScaleFraction <= primarySameScale + params.consensus_gate_max_same_scale_fraction_delta;

S = table();
S.k_neighbors = k;
S.support_threshold = threshold;
S.minimum_co_inclusion_replicates = params.consensus_min_co_inclusion_replicates;
S.n_nodes = n;
S.n_edges = numel(lo);
S.n_components = numel(componentCounts);
S.largest_component_fraction = largestCount ./ n;
S.mean_degree = mean(deg);
S.median_degree = median(deg);
S.p10_degree = quantile(deg, 0.10);
S.stable_min_degree = params.consensus_stable_min_degree;
S.stable_node_fraction = mean(stable);
S.rare_stage_min_degree = params.consensus_rare_min_degree;
S.rare_stage_stable_fraction = rareStableFraction;
S.same_scale_edge_fraction = sameScaleFraction;
S.primary_graph_same_scale_neighbor_fraction = primarySameScale;
S.same_scale_fraction_delta_vs_primary = sameScaleFraction - primarySameScale;
S.same_session_edge_fraction = sameSessionFraction;
S.random_mixing_same_session_fraction = expectedSameSession;
S.same_session_observed_over_null = sessionRatio;
S.gate_stable_node_fraction_pass = gateStable;
S.gate_largest_component_pass = gateComponent;
S.gate_rare_stage_coverage_pass = gateRare;
S.gate_session_composition_pass = gateSession;
S.gate_scale_composition_pass = gateScale;
S.selection_topology_gates_pass = gateStable & gateComponent & gateRare;
S.postselection_composition_audits_pass = gateSession & gateScale;
S.all_handoff_gates_pass = S.selection_topology_gates_pass & ...
    S.postselection_composition_audits_pass;
S.session_or_scale_provenance_used_for_threshold_selection = false;
S.edge_rule = "conditional_any_direction_recurrence_given_endpoint_co_inclusion";
S.labels_used_for_consensus = "none";
S.arena_used_for_consensus = false;
S.condition_used_for_consensus = false;
S.event_channels_used_for_consensus = false;
S.rare_stratum_labels_used_for_consensus = false;

M = struct();
M.summary = S;
M.lo = lo;
M.hi = hi;
M.weight = weight;
M.degree = deg;
M.weightedDegree = weightedDegree;
M.component = component;
M.largestComponentId = largestId;
M.stable = stable;
M.stableInducedDegree = stableInducedDegree;
end

function T = local_node_stability(Meta, stage, eligible, Selected, k, threshold, params, ready)
n = height(Meta);
stable = Selected.stable;
T = table();
T.graph_node_id = double(Meta.graph_node_id);
T.embedding_row_id = Meta.embedding_row_id;
T.scale_index = double(Meta.scale_index);
T.chunk_sec = double(Meta.chunk_sec);
T.anchor_stage = stage;
T.consensus_eligible_replicates = sum(eligible, 2);
T.consensus_degree = Selected.degree;
T.consensus_stable_induced_degree = Selected.stableInducedDegree;
T.consensus_weighted_degree = Selected.weightedDegree;
T.consensus_component_id = Selected.component;
T.consensus_in_largest_component = Selected.component == Selected.largestComponentId;
T.consensus_stable_for_run09 = stable;
T.unstable_for_run09 = ~stable;
T.unstable_reason = repmat("none", n, 1);
T.unstable_reason(Selected.degree < params.consensus_stable_min_degree) = "degree_below_predefined_gate";
T.unstable_reason(Selected.degree >= params.consensus_stable_min_degree & ~stable) = ...
    "removed_by_iterative_minimum_degree_core_or_component_gate";
T.selected_k_neighbors = repmat(k, n, 1);
T.selected_support_threshold = repmat(threshold, n, 1);
T.handoff_ready = repmat(ready, n, 1);
T.stability_role = repmat("run08_graph_stability_not_cluster_assignment", n, 1);
T.labels_used_for_stability = repmat("none", n, 1);
T.arena_used_for_stability = false(n, 1);
T.condition_used_for_stability = false(n, 1);
T.event_channels_used_for_stability = false(n, 1);
T.rare_stratum_labels_used_for_stability = false(n, 1);
end

function T = local_scale_mixing(lo, hi, scales, chunkSec)
directedSource = [lo; hi];
directedTarget = [hi; lo];
uniqueScales = unique(scales, 'stable')';
T = table();
for a = uniqueScales
    sourceNodes = scales == a;
    sourceEdges = scales(directedSource) == a;
    sourceEdgeCount = nnz(sourceEdges);
    for b = uniqueScales
        count = nnz(sourceEdges & scales(directedTarget) == b);
        targetCount = nnz(scales == b);
        if a == b
            expected = max(targetCount - 1, 0) ./ max(numel(scales) - 1, 1);
        else
            expected = targetCount ./ max(numel(scales) - 1, 1);
        end
        one = table();
        one.source_scale_index = a;
        one.source_chunk_sec = median(chunkSec(sourceNodes), 'omitnan');
        one.target_scale_index = b;
        one.target_chunk_sec = median(chunkSec(scales == b), 'omitnan');
        one.n_source_nodes = nnz(sourceNodes);
        one.n_target_nodes = targetCount;
        one.n_source_directed_edges = sourceEdgeCount;
        one.observed_directed_edge_count = count;
        one.observed_neighbor_fraction = count ./ max(sourceEdgeCount, 1);
        one.random_mixing_neighbor_fraction = expected;
        one.observed_over_random_ratio = one.observed_neighbor_fraction ./ max(expected, eps);
        one.mixing_role = "selected_consensus_graph_complete_scale_matrix";
        one.labels_used_for_mixing = "none";
        one.arena_used_for_mixing = false;
        one.condition_used_for_mixing = false;
        T = [T; one]; %#ok<AGROW>
    end
end
end

function T = local_session_sensitivity(lo, hi, weight, n, sessions, params)
sameSession = sessions(lo) == sessions(hi);
variants = ["selected_consensus"; "same_session_edges_removed_audit"];
masks = {true(numel(lo), 1); ~sameSession};
T = table();
for i = 1:2
    mask = masks{i};
    H = graph(lo(mask), hi(mask), weight(mask), n);
    deg = degree(H);
    stableCore = local_stable_core(H, params.consensus_stable_min_degree);
    component = conncomp(H)';
    counts = accumarray(component, 1);
    one = table();
    one.graph_variant = variants(i);
    one.n_nodes = n;
    one.n_edges = nnz(mask);
    one.edge_retention_to_selected_consensus = nnz(mask) ./ max(numel(lo), 1);
    one.n_components = numel(counts);
    one.largest_component_fraction = max(counts) ./ n;
    one.stable_node_fraction = mean(stableCore);
    one.median_degree = median(deg);
    one.same_session_edge_fraction = mean(sameSession(mask), 'omitnan');
    one.session_identity_used_for_edge_filter = i == 2;
    one.audit_only_not_consensus_edge_input = i == 2;
    one.labels_used_for_distance = "none";
    one.arena_used_for_distance = false;
    one.condition_used_for_distance = false;
    T = [T; one]; %#ok<AGROW>
end
end

function T = local_event_coverage(Event, Stability, Meta)
T = table();
if isempty(Event) || height(Event) ~= height(Stability)
    one = table("unavailable_or_row_mismatch", "postfit_event_coverage_not_graph_input", ...
        'VariableNames', {'status','coverage_role'});
    T = one;
    return
end
[tf, loc] = ismember(double(Meta.graph_node_id), double(Event.graph_node_id));
assert(all(tf), 'build_condition_blind_consensus_neighborhood:MissingEventRows', ...
    'eventNodeAudit does not cover every consensus node.');
Event = Event(loc, :);
exclude = ["graph_node_id","embedding_row_id","scale_index","chunk_sec", ...
    "labels_used_for_event_node_audit","arena_used_for_event_node_audit", ...
    "condition_used_for_event_node_audit"];
eventNames = setdiff(string(Event.Properties.VariableNames), exclude, 'stable');
scales = unique(double(Meta.scale_index), 'stable')';
for e = 1:numel(eventNames)
    flag = logical(Event.(eventNames(e)));
    for s = [NaN scales]
        if isnan(s)
            scope = true(height(Meta), 1);
            scopeName = "all_scales";
            sec = NaN;
        else
            scope = double(Meta.scale_index) == s;
            scopeName = "source_scale";
            sec = double(Meta.chunk_sec(find(scope, 1)));
        end
        one = table();
        one.event_id = eventNames(e);
        one.source_scope = scopeName;
        one.scale_index = s;
        one.chunk_sec = sec;
        one.n_nodes = nnz(scope);
        one.n_event_nodes = nnz(scope & flag);
        one.event_fraction = nnz(scope & flag) ./ max(nnz(scope), 1);
        one.n_stable_nodes = nnz(scope & Stability.consensus_stable_for_run09);
        one.n_stable_event_nodes = nnz(scope & flag & Stability.consensus_stable_for_run09);
        one.event_node_stable_fraction = nnz(scope & flag & Stability.consensus_stable_for_run09) ./ ...
            max(nnz(scope & flag), 1);
        one.nonevent_node_stable_fraction = nnz(scope & ~flag & Stability.consensus_stable_for_run09) ./ ...
            max(nnz(scope & ~flag), 1);
        one.median_degree_event_nodes = median(Stability.consensus_degree(scope & flag), 'omitnan');
        one.median_degree_nonevent_nodes = median(Stability.consensus_degree(scope & ~flag), 'omitnan');
        one.coverage_role = "postfit_event_coverage_not_consensus_input";
        one.labels_used_for_consensus = "none";
        one.arena_used_for_consensus = false;
        one.condition_used_for_consensus = false;
        one.event_used_for_consensus_edge = false;
        T = [T; one]; %#ok<AGROW>
    end
end
end

function p = local_same_category_null(values)
[~, ~, g] = unique(string(values));
counts = accumarray(g, 1);
n = numel(values);
p = sum(counts .* (counts - 1)) ./ max(n .* (n - 1), 1);
end

function [stable, inducedDegree] = local_stable_core(H, minimumDegree)
n = numnodes(H);
stable = true(n, 1);
while any(stable)
    ids = find(stable);
    sub = subgraph(H, ids);
    remove = degree(sub) < minimumDegree;
    if ~any(remove)
        break
    end
    stable(ids(remove)) = false;
end
inducedDegree = zeros(n, 1);
if ~any(stable)
    return
end
ids = find(stable);
sub = subgraph(H, ids);
component = conncomp(sub)';
counts = accumarray(component, 1);
[~, largest] = max(counts);
stable(ids(component ~= largest)) = false;
ids = find(stable);
sub = subgraph(H, ids);
inducedDegree(ids) = degree(sub);
end

function T = local_freeze_config(params, n, d, k, threshold, nPassing, ready)
keys = ["consensus_algorithm_version";"n_nodes";"locked_graph_dimensions"; ...
    "replicates";"candidate_k";"consensus_k_values";"selected_k"; ...
    "support_thresholds_audited"; ...
    "selected_support_threshold";"minimum_co_inclusion_replicates"; ...
    "dimension_design_cycle";"core_n_pcs";"core_tail_dimension_fraction"; ...
    "short_prefix_n_pcs";"medium_prefix_n_pcs"; ...
    "stage_selection_rule";"random_seed_base"; ...
    "stable_min_degree";"stable_node_rule";"rare_stage_min_degree"; ...
    "gate_min_stable_node_fraction";"gate_min_largest_component_fraction"; ...
    "gate_min_rare_stage_fraction";"gate_max_same_session_null_ratio"; ...
    "gate_max_same_scale_fraction_delta";"minimum_passing_thresholds"; ...
    "passing_primary_k_topology_thresholds"; ...
    "handoff_ready";"labels_used_for_consensus";"motifs_defined_in_run08"];
values = ["run08_consensus_v1"; string(n); string(d); string(params.consensus_replicates); ...
    string(params.consensus_candidate_k); strjoin(string(params.consensus_k_values), ";"); string(k); ...
    strjoin(string(params.consensus_support_thresholds), ";"); string(threshold); ...
    string(params.consensus_min_co_inclusion_replicates); ...
    "full|medium_prefix|short_prefix|core_plus_tail"; ...
    string(params.consensus_core_n_pcs); string(params.consensus_dimension_fraction); ...
    string(params.consensus_prefix_short_n_pcs); string(params.consensus_prefix_medium_n_pcs); ...
    "all_base_plus_within_scale_dimension_specific_partitioned_rare_stage"; ...
    string(params.audit_random_seed); string(params.consensus_stable_min_degree); ...
    "largest_iterative_minimum_degree_core"; string(params.consensus_rare_min_degree); ...
    string(params.consensus_gate_min_stable_node_fraction); ...
    string(params.consensus_gate_min_largest_component_fraction); ...
    string(params.consensus_gate_min_rare_stable_fraction); ...
    string(params.consensus_gate_max_same_session_null_ratio); ...
    string(params.consensus_gate_max_same_scale_fraction_delta); ...
    string(params.consensus_min_passing_thresholds); string(nPassing); string(ready); ...
    "none"; "false"];
T = table(keys, values, 'VariableNames', {'parameter','effective_value'});
T.freeze_role = repmat("run08_to_run09_reproducibility_contract", height(T), 1);
end

function T = local_handoff_manifest(paths, rows, required, outRoot, ready, k, threshold)
n = numel(paths);
T = table();
T.artifact_id = strings(n, 1);
T.artifact_path = string(paths(:));
T.sha256 = strings(n, 1);
T.file_bytes = zeros(n, 1);
T.n_rows = rows(:);
T.required_for_run09 = logical(required(:));
T.artifact_status = strings(n, 1);
for i = 1:n
    [~, stem, ext] = fileparts(paths(i));
    T.artifact_id(i) = string(stem) + string(ext);
    info = dir(paths(i));
    if isempty(info)
        T.artifact_status(i) = "missing";
        continue
    end
    T.file_bytes(i) = info.bytes;
    T.sha256(i) = compute_file_sha256(paths(i));
    T.artifact_status(i) = "present_hashed_byte_exact";
end
edgeHash = T.sha256(T.artifact_id == "run08_to_run09_edge_list.csv");
configHash = T.sha256(T.artifact_id == "graph_consensus_freeze_config.csv");
freezeId = "run08_consensus_" + extractBefore(edgeHash, 13) + "_" + extractBefore(configHash, 9);
T.consensus_freeze_id = repmat(freezeId, n, 1);
T.output_root = repmat(string(outRoot), n, 1);
T.selected_k_neighbors = repmat(k, n, 1);
T.selected_support_threshold = repmat(threshold, n, 1);
T.handoff_ready = repmat(ready, n, 1);
T.handoff_status = repmat(string(local_handoff_status(ready)), n, 1);
T.run09_may_modify_edges = false(n, 1);
T.run09_may_cluster_only_stable_nodes = true(n, 1);
T.motifs_defined_in_run08 = false(n, 1);
T.labels_used_for_consensus = repmat("none", n, 1);
T.arena_used_for_consensus = false(n, 1);
T.condition_used_for_consensus = false(n, 1);
T.hash_algorithm = repmat("SHA-256", n, 1);
T.hash_scope = repmat("exact_file_bytes", n, 1);
T.hash_implementation = repmat( ...
    "matlab_fread_uint8_to_java_message_digest", n, 1);
T.manifest_schema_version = repmat( ...
    "run08_handoff_manifest_v2_byte_exact", n, 1);
T.artifact_role = repmat("run08_consensus_audit", n, 1);
T.permitted_run09_use = repmat("audit_only_not_graph_input", n, 1);
inputMask = ismember(T.artifact_id, ["run08_to_run09_edge_list.csv", ...
    "run08_to_run09_node_input.csv"]);
T.artifact_role(inputMask) = "frozen_run09_graph_input";
T.permitted_run09_use(T.artifact_id == "run08_to_run09_edge_list.csv") = ...
    "immutable_weighted_edge_input";
T.permitted_run09_use(T.artifact_id == "run08_to_run09_node_input.csv") = ...
    "immutable_stable_node_identity_scale_provenance_only";
configMask = T.artifact_id == "graph_consensus_freeze_config.csv";
T.artifact_role(configMask) = "frozen_run09_graph_contract";
T.permitted_run09_use(configMask) = "read_only_reproducibility_contract";
provenanceMask = ismember(T.artifact_id, ["graph_edge_list.csv", ...
    "graph_node_manifest.csv","graph_event_node_audit.csv", ...
    "graph_score_preprocess_audit.csv","graph_input_manifest.csv"]);
T.artifact_role(provenanceMask) = "locked_run08_input_provenance";
T.permitted_run09_use(provenanceMask) = "audit_only_not_graph_input";
sourceMask = endsWith(T.artifact_id, ".m") | ...
    T.artifact_id == "multiscale_graph_config.csv";
T.artifact_role(sourceMask) = "source_code_or_config_provenance";
T.permitted_run09_use(sourceMask) = "audit_only_exact_implementation_trace";
end

function text = local_handoff_status(ready)
if ready
    text = "ready_frozen_consensus_graph";
else
    text = "not_ready_predefined_graph_stability_gates_failed";
end
end
