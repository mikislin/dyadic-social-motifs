function outputs = run_condition_blind_motif_graph_build(repoRoot, opts)
%RUN_CONDITION_BLIND_MOTIF_GRAPH_BUILD Build run-08 graph/topology artifacts.
%
% Condition-blind contract
%   - Graph edges are built only from numeric run_07 global embedding scores.
%   - Condition, cohort, group, drug, genotype, outcome, arena, and event
%     labels do not enter score preprocessing, kNN graph construction, or
%     topology sensitivity.
%   - Arena and event channels are joined only after graph construction for
%     audit/coverage diagnostics.

if nargin < 1 || strlength(string(repoRoot)) == 0
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
end
if nargin < 2 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'configPath') || strlength(string(opts.configPath)) == 0
    opts.configPath = fullfile(repoRoot, 'config', 'multiscale_graph_config.csv');
end

repoRoot = string(repoRoot);
cd(repoRoot);
addpath(genpath(repoRoot));

params = load_multiscale_graph_config(opts.configPath);
outRoot = resolve_repo_path(repoRoot, params.output_dir);
figDir = fullfile(outRoot, 'figures');
logDir = fullfile(outRoot, 'logs');
local_ensure_dir(outRoot);
local_ensure_dir(figDir);
local_ensure_dir(logDir);

paths = local_output_paths(outRoot, params);
diary(fullfile(logDir, 'run_08_build_condition_blind_motif_graph_latest.log'));
cleanup = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('run_08_build_condition_blind_motif_graph\n');
fprintf('Repo root: %s\n', repoRoot);
fprintf('Output root: %s\n', outRoot);
fprintf('Run mode: %s\n', params.run_mode);
fprintf('Anchor manifest mode: %s\n', params.anchor_manifest_mode);

embeddingRoot = resolve_repo_path(repoRoot, params.embedding_input_dir);
fprintf('Run_07 embedding input root: %s\n', embeddingRoot);
local_assert_required_embedding_files(embeddingRoot, params);

writetable(params.config_table, paths.parameterAudit);

[scoreTable, rowManifest, inputScaleAudit, featureDictionary] = local_load_inputs(embeddingRoot, params);
local_assert_no_label_value_columns(scoreTable, "run_07 global score table");
[scoreTable, rowManifest] = local_align_score_and_rows(scoreTable, rowManifest);

[Xraw, scoreNames] = local_graph_score_matrix(scoreTable, params);
[X, scorePreprocessAudit] = local_preprocess_graph_scores(Xraw, scoreNames, params);
nodeManifest = local_node_manifest(scoreTable, scoreNames);

inputManifest = local_input_manifest(embeddingRoot, params, scoreTable, rowManifest, ...
    inputScaleAudit, scoreNames, Xraw, X);
writetable(inputManifest, paths.inputManifest);
writetable(scorePreprocessAudit, paths.scorePreprocessAudit);
writetable(nodeManifest, paths.nodeManifest);

maxK = max([params.k_neighbors, params.k_sensitivity_values(:)']);
buildParams = params;
buildParams.k_neighbors = maxK;
GraphMax = build_condition_blind_knn_graph(X, nodeManifest, buildParams);
Graph = local_subset_graph(GraphMax, params.k_neighbors);

[topologySummary, degreeAudit, componentAudit, componentId] = ...
    local_topology_audits(Graph, nodeManifest, params);
[arenaLabel, arenaSource] = local_posthoc_arena_labels(rowManifest, nodeManifest);
neighborAudit = local_neighbor_composition_audit(Graph, nodeManifest, arenaLabel);
[scaleMixingAudit, neighborNullAudit] = ...
    compute_graph_scale_null_audits(Graph, nodeManifest, arenaLabel);
kSensitivity = local_k_sensitivity_audit(GraphMax, params, nodeManifest);

[eventCoverageAudit, eventNodeAudit, eventSource] = local_event_coverage_audit( ...
    rowManifest, nodeManifest, Graph, inputScaleAudit, params);

writetable(Graph.Edges, paths.edgeList);
writetable(topologySummary, paths.topologySummary);
writetable(degreeAudit, paths.degreeAudit);
writetable(componentAudit, paths.componentAudit);
writetable(neighborAudit, paths.neighborAudit);
writetable(scaleMixingAudit, paths.scaleMixingAudit);
writetable(neighborNullAudit, paths.neighborNullAudit);
writetable(kSensitivity, paths.kSensitivityAudit);
writetable(eventCoverageAudit, paths.eventCoverageAudit);
writetable(eventNodeAudit, paths.eventNodeAudit);

eventChunkRoot = string(fileparts(char(eventSource)));
[rareDefinition, rareMembership, rareSeedManifest] = define_condition_blind_rare_strata(outRoot, ...
    'embeddingRoot', embeddingRoot, ...
    'chunkRoot', eventChunkRoot, ...
    'eventSummaryFile', local_active_event_summary_file(params), ...
    'nodeManifest', nodeManifest, ...
    'degreeAudit', degreeAudit, ...
    'rowManifest', rowManifest, ...
    'lowDensityQuantile', params.rare_low_density_quantile, ...
    'peripheryInDegreeQuantile', params.rare_periphery_indegree_quantile, ...
    'highRadialQuantile', params.high_quantile_threshold, ...
    'undercoveredCellQuantile', params.rare_undercovered_cell_quantile, ...
    'writeOutputs', false);
rareDefinition.definition_stage = repmat("postfit_current_graph_re_evaluation", height(rareDefinition), 1);
rareDefinition.active_for_anchor_selection = false(height(rareDefinition), 1);
rareMembership.definition_stage = repmat("postfit_current_graph_re_evaluation", height(rareMembership), 1);
rareMembership.active_for_anchor_selection = false(height(rareMembership), 1);
rareSeedManifest.definition_stage = repmat("postfit_current_graph_re_evaluation", height(rareSeedManifest), 1);
rareSeedManifest.active_for_anchor_selection = false(height(rareSeedManifest), 1);
rareCompositionAudit = local_rare_composition_audit(rareDefinition, rareMembership, nodeManifest, params);
rareNeighborAudit = local_rare_neighbor_retention_audit(rareMembership, nodeManifest, Graph, params);
baselineCoverageAudit = local_baseline_vs_enriched_coverage_audit(repoRoot, eventCoverageAudit, params);
eventPrevalenceFoldAudit = compute_graph_event_prevalence_fold_audit(baselineCoverageAudit);
writetable(rareCompositionAudit, paths.rareCompositionAudit);
writetable(rareNeighborAudit, paths.rareNeighborAudit);
writetable(baselineCoverageAudit, paths.baselineCoverageAudit);
writetable(eventPrevalenceFoldAudit, paths.eventPrevalenceFoldAudit);
rareDefinitionStageAudit = persist_run08_rare_definition_stages( ...
    repoRoot, outRoot, params, rareDefinition, rareMembership, rareSeedManifest);

fprintf('Building anchor-stage, session-excluded, and resampling graph audits...\n');
sensitivityAudit = build_condition_blind_graph_sensitivity_audits( ...
    X, nodeManifest, Graph, params, outRoot);
fprintf('Building global PCA and visualization-only UMAP audits...\n');
visualizationAudit = build_run08_embedding_visualization_audits( ...
    embeddingRoot, outRoot, X, nodeManifest, Graph, params);

GraphModel = struct();
GraphModel.params = params;
GraphModel.score_names = scoreNames;
GraphModel.X = X;
GraphModel.nodeManifest = nodeManifest;
GraphModel.Edges = Graph.Edges;
GraphModel.component_id = componentId;
GraphModel.topologySummary = topologySummary;
GraphModel.degreeAudit = degreeAudit;
GraphModel.componentAudit = componentAudit;
GraphModel.neighborAudit = neighborAudit;
GraphModel.scaleMixingAudit = scaleMixingAudit;
GraphModel.neighborNullAudit = neighborNullAudit;
GraphModel.kSensitivity = kSensitivity;
GraphModel.eventCoverageAudit = eventCoverageAudit;
GraphModel.eventNodeAudit = eventNodeAudit;
GraphModel.rareDefinition = rareDefinition;
GraphModel.rareMembership = rareMembership;
GraphModel.rareSeedManifest = rareSeedManifest;
GraphModel.rareCompositionAudit = rareCompositionAudit;
GraphModel.rareNeighborAudit = rareNeighborAudit;
GraphModel.baselineCoverageAudit = baselineCoverageAudit;
GraphModel.eventPrevalenceFoldAudit = eventPrevalenceFoldAudit;
GraphModel.rareDefinitionProvenanceAudit = rareDefinitionStageAudit.provenance;
GraphModel.rareDefinitionComparisonAudit = rareDefinitionStageAudit.comparison;
GraphModel.anchorStageSensitivityAudit = sensitivityAudit.stageSummary;
GraphModel.sessionExcludedSensitivityAudit = sensitivityAudit.sessionSummary;
GraphModel.neighborhoodResamplingAudit = sensitivityAudit.resampling;
GraphModel.globalPcaVarianceAudit = visualizationAudit.globalVariance;
GraphModel.umapStatusAudit = visualizationAudit.umapStatus;
GraphModel.labels_used_for_graph = "none";
GraphModel.arena_used_for_graph = false;
GraphModel.condition_used_for_graph = false;
save(paths.graphModelMat, 'GraphModel', '-v7.3');

figureManifest = make_run08_graph_qc_figures(outRoot, params);
writetable(figureManifest, paths.figureManifest);

fprintf('Graph nodes: %d\n', height(nodeManifest));
fprintf('Primary k: %d | directed edges: %d\n', Graph.k, height(Graph.Edges));
fprintf('Arena audit source: %s\n', arenaSource);
fprintf('Event audit source: %s\n', eventSource);
fprintf('Wrote parameter audit: %s\n', paths.parameterAudit);
fprintf('Wrote input manifest: %s\n', paths.inputManifest);
fprintf('Wrote node manifest: %s\n', paths.nodeManifest);
fprintf('Wrote edge list: %s\n', paths.edgeList);
fprintf('Wrote topology summary: %s\n', paths.topologySummary);
fprintf('Wrote rare-event coverage audit: %s\n', paths.eventCoverageAudit);
fprintf('Wrote complete scale-mixing audit: %s\n', paths.scaleMixingAudit);
fprintf('Wrote null-normalized neighbor audit: %s\n', paths.neighborNullAudit);
fprintf('Wrote rare-strata composition audit: %s\n', paths.rareCompositionAudit);
fprintf('Wrote rare-strata neighbor-retention audit: %s\n', paths.rareNeighborAudit);
fprintf('Wrote baseline-versus-enriched coverage audit: %s\n', paths.baselineCoverageAudit);
fprintf('Wrote event prevalence-fold audit: %s\n', paths.eventPrevalenceFoldAudit);
fprintf('Wrote figure manifest: %s\n', paths.figureManifest);
fprintf('Wrote ignored MAT artifact: %s\n', paths.graphModelMat);

outputs = struct();
outputs.output_root = string(outRoot);
outputs.embedding_input_root = string(embeddingRoot);
outputs.n_nodes = height(nodeManifest);
outputs.k_neighbors = Graph.k;
outputs.n_directed_edges = height(Graph.Edges);
outputs.topology_summary_path = string(paths.topologySummary);
outputs.edge_list_path = string(paths.edgeList);
outputs.event_coverage_audit_path = string(paths.eventCoverageAudit);
outputs.rare_composition_audit_path = string(paths.rareCompositionAudit);
outputs.rare_neighbor_retention_audit_path = string(paths.rareNeighborAudit);
outputs.baseline_coverage_audit_path = string(paths.baselineCoverageAudit);
outputs.figure_manifest_path = string(paths.figureManifest);
outputs.graph_model_mat_path = string(paths.graphModelMat);
end

function paths = local_output_paths(outRoot, params)
paths = struct();
paths.parameterAudit = fullfile(outRoot, 'graph_parameter_audit.csv');
paths.inputManifest = fullfile(outRoot, 'graph_input_manifest.csv');
paths.scorePreprocessAudit = fullfile(outRoot, 'graph_score_preprocess_audit.csv');
paths.nodeManifest = fullfile(outRoot, 'graph_node_manifest.csv');
paths.edgeList = fullfile(outRoot, 'graph_edge_list.csv');
paths.topologySummary = fullfile(outRoot, 'graph_topology_summary.csv');
paths.degreeAudit = fullfile(outRoot, 'graph_degree_audit.csv');
paths.componentAudit = fullfile(outRoot, 'graph_component_audit.csv');
paths.neighborAudit = fullfile(outRoot, 'graph_neighbor_composition_audit.csv');
paths.scaleMixingAudit = fullfile(outRoot, 'graph_scale_mixing_matrix_audit.csv');
paths.neighborNullAudit = fullfile(outRoot, 'graph_neighbor_null_normalized_audit.csv');
paths.kSensitivityAudit = fullfile(outRoot, 'graph_k_sensitivity_audit.csv');
paths.eventCoverageAudit = fullfile(outRoot, 'graph_rare_event_coverage_audit.csv');
paths.eventNodeAudit = fullfile(outRoot, 'graph_event_node_audit.csv');
paths.rareCompositionAudit = fullfile(outRoot, 'graph_rare_strata_composition_audit.csv');
paths.rareNeighborAudit = fullfile(outRoot, 'graph_rare_strata_neighbor_retention_audit.csv');
paths.baselineCoverageAudit = fullfile(outRoot, 'graph_baseline_vs_enriched_coverage_audit.csv');
paths.eventPrevalenceFoldAudit = fullfile(outRoot, 'graph_event_prevalence_fold_audit.csv');
paths.figureManifest = fullfile(outRoot, 'graph_qc_figure_manifest.csv');
paths.graphModelMat = fullfile(outRoot, char(params.graph_model_mat_file));
end

function local_assert_required_embedding_files(embeddingRoot, params)
required = [params.embedding_scores_file, params.embedding_row_manifest_file, ...
    params.embedding_input_scale_audit_file, params.embedding_feature_dictionary_file];
for i = 1:numel(required)
    p = fullfile(embeddingRoot, char(required(i)));
    if ~isfile(p)
        if params.run_mode == "full"
            error('run_condition_blind_motif_graph_build:MissingRun07Input', ...
                ['Missing required full run_07 input: %s\n' ...
                'Full run_08 reads embedding_input_dir=%s.\n' ...
                'Run full run_07 first, or switch run_08 to smoke mode with ' ...
                'RUN08_GRAPH_RUN_MODE=smoke and smoke input/output roots.'], ...
                p, string(params.embedding_input_dir));
        else
            error('run_condition_blind_motif_graph_build:MissingRun07Input', ...
                ['Missing required smoke run_07 input: %s\n' ...
                'Smoke run_08 reads embedding_input_dir=%s. Run smoke run_07 first.'], ...
                p, string(params.embedding_input_dir));
        end
    end
end
end

function [scoreTable, rowManifest, inputScaleAudit, featureDictionary] = local_load_inputs(embeddingRoot, params)
scoreTable = local_read_csv(fullfile(embeddingRoot, char(params.embedding_scores_file)));
rowManifest = local_read_csv(fullfile(embeddingRoot, char(params.embedding_row_manifest_file)));
inputScaleAudit = local_read_csv(fullfile(embeddingRoot, char(params.embedding_input_scale_audit_file)));
featureDictionary = local_read_csv(fullfile(embeddingRoot, char(params.embedding_feature_dictionary_file)));
end

function [scoreTable, rowManifest] = local_align_score_and_rows(scoreTable, rowManifest)
assert(ismember('embedding_row_id', scoreTable.Properties.VariableNames) && ...
    ismember('embedding_row_id', rowManifest.Properties.VariableNames), ...
    'run_condition_blind_motif_graph_build:MissingEmbeddingRowId', ...
    'Run_07 score and row tables must contain embedding_row_id.');
[tf, loc] = ismember(scoreTable.embedding_row_id, rowManifest.embedding_row_id);
assert(all(tf), 'run_condition_blind_motif_graph_build:MissingRowManifestRows', ...
    'Some run_07 score rows could not be matched to the row manifest.');
rowManifest = rowManifest(loc, :);
end

function local_assert_no_label_value_columns(T, context)
names = string(T.Properties.VariableNames);
lowerNames = lower(names);
blocked = ["condition", "cohort", "group", "drug", "genotype", "outcome"];
bad = false(size(names));
for i = 1:numel(blocked)
    bad = bad | contains(lowerNames, blocked(i));
end
auditFlag = startsWith(lowerNames, "labels_used") | contains(lowerNames, "_used_for_");
bad(auditFlag) = false;
assert(~any(bad), 'run_condition_blind_motif_graph_build:LabelValueColumn', ...
    '%s contains label-like value columns: %s', context, strjoin(names(bad), ', '));
end

function [X, scoreNames] = local_graph_score_matrix(scoreTable, params)
allNames = string(scoreTable.Properties.VariableNames);
scoreNames = allNames(startsWith(allNames, "global_pc"));
assert(numel(scoreNames) >= params.graph_n_pcs, ...
    'run_condition_blind_motif_graph_build:InsufficientGlobalPCs', ...
    'Requested %d graph PCs, but only %d global_pc columns are present.', ...
    params.graph_n_pcs, numel(scoreNames));
scoreNames = scoreNames(1:params.graph_n_pcs);
X = zeros(height(scoreTable), numel(scoreNames));
for j = 1:numel(scoreNames)
    X(:, j) = double(scoreTable.(scoreNames(j)));
end
assert(all(isfinite(X), 'all'), ...
    'run_condition_blind_motif_graph_build:NonFiniteScores', ...
    'Selected global PC scores contain non-finite values.');
end

function [Xout, audit] = local_preprocess_graph_scores(X, scoreNames, params)
Xout = double(X);
nDim = size(Xout, 2);
center = zeros(nDim, 1);
scale = ones(nDim, 1);
inputStd = zeros(nDim, 1);
outputStd = zeros(nDim, 1);
if ~logical(params.standardize_graph_scores)
    audit = table((1:nDim)', scoreNames(:), center, scale, inputStd, std(Xout, 0, 1)', ...
        'VariableNames', {'graph_pc_index','score_name','robust_center','robust_scale', ...
        'input_std','output_std'});
    audit.standardized_for_graph = false(height(audit), 1);
    audit.labels_used_for_score_preprocessing = repmat("none", height(audit), 1);
    audit.arena_used_for_score_preprocessing = false(height(audit), 1);
    audit.condition_used_for_score_preprocessing = false(height(audit), 1);
    return
end
for j = 1:nDim
    x = Xout(:, j);
    inputStd(j) = std(x, 0, 'omitnan');
    med = median(x, 'omitnan');
    sc = iqr(x);
    if ~(isfinite(sc) && sc > 1e-10)
        sc = std(x, 0, 'omitnan');
    end
    if ~(isfinite(sc) && sc > 1e-10)
        Xout(:, j) = 0;
        center(j) = med;
        scale(j) = 1;
        outputStd(j) = 0;
        continue
    end
    z = (x - med) ./ sc;
    z = min(max(z, -params.graph_score_winsor_abs), params.graph_score_winsor_abs);
    Xout(:, j) = z;
    center(j) = med;
    scale(j) = sc;
    outputStd(j) = std(z, 0, 'omitnan');
end
audit = table((1:nDim)', scoreNames(:), center, scale, inputStd, outputStd, ...
    'VariableNames', {'graph_pc_index','score_name','robust_center','robust_scale', ...
    'input_std','output_std'});
audit.standardized_for_graph = true(height(audit), 1);
audit.score_winsor_abs = repmat(params.graph_score_winsor_abs, height(audit), 1);
audit.labels_used_for_score_preprocessing = repmat("none", height(audit), 1);
audit.arena_used_for_score_preprocessing = false(height(audit), 1);
audit.condition_used_for_score_preprocessing = false(height(audit), 1);
end

function nodeManifest = local_node_manifest(scoreTable, scoreNames)
candidate = ["embedding_row_id", "scale_index", "primary_scale_rank", "chunk_sec", ...
    "chunk_frames", "session_index", "raw_index", "anchor_frame", "anchor_time_s", ...
    "start_time_s", "stop_time_s", "valid_frac", "feature_finite_frac", ...
    "hierarchical_role", "session_id", "qc_set", "primary_anchor_global_id", ...
    "expanded_anchor_global_id", "anchor_manifest_mode", "anchor_stage", ...
    "selection_phase", "quota_requested_stratum_id", ...
    "rare_stratum_id", "final_assigned_rare_stratum_id", ...
    "rare_stratum_rule", "rare_stratum_score", ...
    "rare_strata_membership_ids", "rare_strata_membership_count", ...
    "selection_composite_score", "fill_composite_score", "fill_reason", ...
    "base_inclusion_probability", "rare_inclusion_probability", ...
    "final_inclusion_probability", "pca_training_weight", ...
    "graph_training_weight", "audit_inverse_probability_weight", ...
    "audit_weight_interpretation", ...
    "run07_embedding_anchor_id", "run07_matrix_role", "labels_used_for_embedding_matrix", ...
    "arena_used_for_embedding_matrix", "condition_used_for_embedding_matrix"];
keep = candidate(ismember(candidate, string(scoreTable.Properties.VariableNames)));
nodeManifest = scoreTable(:, keep);
nodeManifest.graph_node_id = (1:height(nodeManifest))';
nodeManifest = movevars(nodeManifest, 'graph_node_id', 'Before', 1);
if ismember("global_pc01", scoreNames)
    nodeManifest.graph_plot_pc1 = double(scoreTable.global_pc01);
end
if ismember("global_pc02", scoreNames)
    nodeManifest.graph_plot_pc2 = double(scoreTable.global_pc02);
end
if ismember("global_pc03", scoreNames)
    nodeManifest.graph_plot_pc3 = double(scoreTable.global_pc03);
end
nodeManifest.labels_used_for_graph_node = repmat("none", height(nodeManifest), 1);
nodeManifest.arena_used_for_graph_node = false(height(nodeManifest), 1);
nodeManifest.condition_used_for_graph_node = false(height(nodeManifest), 1);
end

function inputManifest = local_input_manifest(embeddingRoot, params, scoreTable, rowManifest, ...
    inputScaleAudit, scoreNames, Xraw, X)
inputManifest = table();
inputManifest.embedding_input_root = string(embeddingRoot);
inputManifest.embedding_scores_csv = string(fullfile(embeddingRoot, char(params.embedding_scores_file)));
inputManifest.embedding_row_manifest_csv = string(fullfile(embeddingRoot, char(params.embedding_row_manifest_file)));
inputManifest.n_score_rows = height(scoreTable);
inputManifest.n_row_manifest_rows = height(rowManifest);
inputManifest.n_graph_score_columns_available = nnz(startsWith(string(scoreTable.Properties.VariableNames), "global_pc"));
inputManifest.n_graph_score_columns_used = numel(scoreNames);
inputManifest.graph_score_columns_used = strjoin(scoreNames, ";");
inputManifest.raw_score_finite_fraction = mean(isfinite(Xraw), 'all');
inputManifest.graph_input_finite_fraction = mean(isfinite(X), 'all');
inputManifest.run07_input_scale_count = height(inputScaleAudit);
if ismember('used_reviewed_snapshot_fallback', inputScaleAudit.Properties.VariableNames)
    inputManifest.run07_used_reviewed_snapshot_fallback = any(logical(inputScaleAudit.used_reviewed_snapshot_fallback));
else
    inputManifest.run07_used_reviewed_snapshot_fallback = false;
end
inputManifest.graph_rule = string(params.graph_rule);
inputManifest.labels_used_for_graph_input = "none";
inputManifest.arena_used_for_graph_input = false;
inputManifest.condition_used_for_graph_input = false;
end

function Graph = local_subset_graph(GraphMax, k)
k = min(floor(k), GraphMax.k);
Graph = GraphMax;
Graph.Edges = GraphMax.Edges(GraphMax.Edges.neighbor_rank <= k, :);
Graph.k = k;
Graph.distance_scale = median(Graph.Edges.neighbor_distance(Graph.Edges.neighbor_distance > 0), 'omitnan');
if ~(isfinite(Graph.distance_scale) && Graph.distance_scale > 0)
    Graph.distance_scale = 1;
end
Graph.Edges.edge_weight = exp(-Graph.Edges.neighbor_distance ./ Graph.distance_scale);
end

function [summary, degreeAudit, componentAudit, componentId] = local_topology_audits(Graph, nodeManifest, params)
n = Graph.n_nodes;
Edges = Graph.Edges;
inDegree = accumarray(Edges.target_node_id, 1, [n 1], @sum, 0);
outDegree = accumarray(Edges.source_node_id, 1, [n 1], @sum, 0);
mutualCount = accumarray(Edges.source_node_id, double(Edges.is_mutual_neighbor), [n 1], @sum, 0);
meanDistance = accumarray(Edges.source_node_id, Edges.neighbor_distance, [n 1], @mean, NaN);
radiusK = accumarray(Edges.source_node_id, Edges.neighbor_distance, [n 1], @max, NaN);

[sUndir, tUndir] = local_undirected_edge_pairs(Edges);
G = graph(sUndir, tUndir, [], n);
componentId = conncomp(G)';
undirectedDegree = degree(G);

degreeAudit = nodeManifest(:, {'graph_node_id', 'embedding_row_id', 'scale_index', 'chunk_sec'});
degreeAudit.knn_in_degree = inDegree;
degreeAudit.knn_out_degree = outDegree;
degreeAudit.undirected_degree = undirectedDegree;
degreeAudit.mutual_neighbor_fraction = mutualCount ./ max(outDegree, 1);
degreeAudit.mean_neighbor_distance = meanDistance;
degreeAudit.knn_radius = radiusK;
degreeAudit.component_id = componentId;
degreeAudit.labels_used_for_degree_audit = repmat("none", n, 1);
degreeAudit.arena_used_for_degree_audit = false(n, 1);
degreeAudit.condition_used_for_degree_audit = false(n, 1);

u = unique(componentId, 'stable');
componentAudit = table();
for i = 1:numel(u)
    idx = componentId == u(i);
    one = table();
    one.component_id = u(i);
    one.n_nodes = nnz(idx);
    one.node_fraction = nnz(idx) ./ n;
    one.n_scales = numel(unique(nodeManifest.scale_index(idx)));
    one.n_sessions = numel(unique(nodeManifest.session_index(idx)));
    one.mean_undirected_degree = mean(undirectedDegree(idx), 'omitnan');
    one.median_knn_radius = median(radiusK(idx), 'omitnan');
    one.labels_used_for_component_audit = "none";
    one.arena_used_for_component_audit = false;
    one.condition_used_for_component_audit = false;
    componentAudit = [componentAudit; one]; %#ok<AGROW>
end
componentAudit = sortrows(componentAudit, 'n_nodes', 'descend');

summary = table();
summary.n_nodes = n;
summary.n_graph_dims = Graph.n_dims;
summary.k_neighbors = Graph.k;
summary.n_directed_edges = height(Edges);
summary.n_undirected_edges = numel(sUndir);
summary.n_components = height(componentAudit);
summary.largest_component_fraction = max(componentAudit.node_fraction);
summary.median_in_degree = median(inDegree, 'omitnan');
summary.median_undirected_degree = median(undirectedDegree, 'omitnan');
summary.median_mutual_neighbor_fraction = median(degreeAudit.mutual_neighbor_fraction, 'omitnan');
summary.median_knn_radius = median(radiusK, 'omitnan');
summary.graph_score_winsor_abs = params.graph_score_winsor_abs;
summary.labels_used_for_graph = "none";
summary.arena_used_for_graph = false;
summary.condition_used_for_graph = false;
summary.motifs_defined_in_run08 = false;
end

function [arenaLabel, sourceText] = local_posthoc_arena_labels(rowManifest, nodeManifest)
arenaLabel = strings(height(nodeManifest), 1);
sourceText = "unavailable";
if ismember('arena_label', rowManifest.Properties.VariableNames)
    [tf, loc] = ismember(nodeManifest.embedding_row_id, rowManifest.embedding_row_id);
    if all(tf)
        arenaLabel = string(rowManifest.arena_label(loc));
        sourceText = "embedding_row_manifest.csv#arena_label";
    end
end
end

function audit = local_neighbor_composition_audit(Graph, nodeManifest, arenaLabel)
Edges = Graph.Edges;
sourceScale = nodeManifest.scale_index(Edges.source_node_id);
targetScale = nodeManifest.scale_index(Edges.target_node_id);
sourceSession = nodeManifest.session_index(Edges.source_node_id);
targetSession = nodeManifest.session_index(Edges.target_node_id);
sameScale = sourceScale == targetScale;
sameSession = sourceSession == targetSession;
sameArena = false(height(Edges), 1);
if any(arenaLabel ~= "")
    sameArena = arenaLabel(Edges.source_node_id) == arenaLabel(Edges.target_node_id);
end

scales = unique(nodeManifest.scale_index, 'stable')';
audit = table();
for s = scales
    idxNode = nodeManifest.scale_index == s;
    idxEdge = idxNode(Edges.source_node_id);
    one = table();
    one.scale_index = double(s);
    one.chunk_sec = double(nodeManifest.chunk_sec(find(idxNode, 1)));
    one.n_nodes = nnz(idxNode);
    one.n_edges_from_scale = nnz(idxEdge);
    one.mean_same_scale_neighbor_fraction = mean(sameScale(idxEdge), 'omitnan');
    one.mean_same_session_neighbor_fraction = mean(sameSession(idxEdge), 'omitnan');
    one.mean_same_arena_neighbor_fraction_posthoc = mean(sameArena(idxEdge), 'omitnan');
    one.mean_neighbor_distance = mean(Edges.neighbor_distance(idxEdge), 'omitnan');
    one.labels_used_for_neighbor_audit = "none";
    one.arena_used_for_graph = false;
    one.arena_used_for_posthoc_neighbor_audit = any(arenaLabel ~= "");
    one.condition_used_for_neighbor_audit = false;
    audit = [audit; one]; %#ok<AGROW>
end
end

function audit = local_k_sensitivity_audit(GraphMax, params, nodeManifest)
primary = local_subset_graph(GraphMax, params.k_neighbors);
primaryKeys = local_undirected_edge_keys(primary.Edges, GraphMax.n_nodes);
kValues = unique([params.k_sensitivity_values(:); params.k_neighbors], 'stable')';
audit = table();
for k = kValues
    Gk = local_subset_graph(GraphMax, k);
    keys = local_undirected_edge_keys(Gk.Edges, GraphMax.n_nodes);
    inter = numel(intersect(primaryKeys, keys));
    uni = numel(union(primaryKeys, keys));
    [sUndir, tUndir] = local_undirected_edge_pairs(Gk.Edges);
    H = graph(sUndir, tUndir, [], GraphMax.n_nodes);
    bins = conncomp(H);
    one = table();
    one.k_neighbors = double(k);
    one.n_directed_edges = height(Gk.Edges);
    one.n_undirected_edges = numel(keys);
    one.edge_jaccard_to_primary_k = inter ./ max(uni, 1);
    one.n_components = numel(unique(bins));
    one.largest_component_fraction = max(accumarray(bins(:), 1)) ./ GraphMax.n_nodes;
    one.median_knn_radius = median(accumarray(Gk.Edges.source_node_id, ...
        Gk.Edges.neighbor_distance, [GraphMax.n_nodes 1], @max, NaN), 'omitnan');
    one.n_scales = numel(unique(nodeManifest.scale_index));
    one.labels_used_for_k_sensitivity = "none";
    one.arena_used_for_k_sensitivity = false;
    one.condition_used_for_k_sensitivity = false;
    audit = [audit; one]; %#ok<AGROW>
end
end

function [coverageAudit, eventNodeAudit, eventSource] = local_event_coverage_audit( ...
    rowManifest, nodeManifest, Graph, inputScaleAudit, params)
eventSource = "unavailable";
coverageAudit = table();
eventNodeAudit = nodeManifest(:, {'graph_node_id', 'embedding_row_id', 'scale_index', 'chunk_sec'});
eventNodeAudit.labels_used_for_event_node_audit = repmat("none", height(eventNodeAudit), 1);
eventNodeAudit.arena_used_for_event_node_audit = false(height(eventNodeAudit), 1);
eventNodeAudit.condition_used_for_event_node_audit = false(height(eventNodeAudit), 1);

if ~ismember('run06_input_root', inputScaleAudit.Properties.VariableNames)
    return
end
inputRoot = string(inputScaleAudit.run06_input_root(find(strlength(string(inputScaleAudit.run06_input_root)) > 0, 1)));
if strlength(inputRoot) == 0
    return
end
eventPath = fullfile(inputRoot, char(local_active_event_summary_file(params)));
if ~isfile(eventPath)
    return
end
E = local_read_csv(eventPath);
eventSource = string(eventPath);
local_assert_event_summary_label_free(E);

if ismember('expanded_anchor_global_id', rowManifest.Properties.VariableNames) && ...
        ismember('expanded_anchor_global_id', E.Properties.VariableNames) && ...
        all(isfinite(double(rowManifest.expanded_anchor_global_id)))
    joinKey = 'expanded_anchor_global_id';
elseif ismember('primary_anchor_global_id', rowManifest.Properties.VariableNames) && ...
        ismember('primary_anchor_global_id', E.Properties.VariableNames)
    joinKey = 'primary_anchor_global_id';
else
    error('run_condition_blind_motif_graph_build:MissingEventJoinKey', ...
        'Event coverage audit has no shared condition-blind anchor identifier.');
end
[tf, loc] = ismember(rowManifest.(joinKey), E.(joinKey));
assert(all(tf), 'run_condition_blind_motif_graph_build:MissingEventRows', ...
    'Some run_07 row-manifest rows could not be matched to run_06 event summary.');
Enode = E(loc, :);

eventDefs = local_event_definitions();
for e = 1:height(eventDefs)
    eventNodeAudit.(eventDefs.event_id(e)) = false(height(nodeManifest), 1);
end

Edges = Graph.Edges;
scales = unique(nodeManifest.scale_index, 'stable')';
for e = 1:height(eventDefs)
    nodeFlagAll = false(height(nodeManifest), 1);
    for s = scales
        fullScale = E.scale_index == s;
        nodeScale = nodeManifest.scale_index == s;
        [fullFlag, threshold] = local_event_flag(E, fullScale, eventDefs(e, :), params);
        [nodeFlag, ~] = local_event_flag(Enode, nodeScale, eventDefs(e, :), params, threshold);
        nodeFlagAll(nodeScale) = nodeFlag(nodeScale);

        validFull = fullScale & E.event_valid_fraction >= params.min_event_valid_fraction;
        validNode = nodeScale & Enode.event_valid_fraction >= params.min_event_valid_fraction;
        fullEvent = validFull & fullFlag;
        nodeEvent = validNode & nodeFlag;
        neighborFrac = local_neighbor_event_fraction(Edges, nodeFlag);

        one = table();
        one.event_id = eventDefs.event_id(e);
        one.event_rule = eventDefs.event_rule(e);
        one.source_column = eventDefs.source_column(e);
        one.scale_index = double(s);
        one.chunk_sec = double(nodeManifest.chunk_sec(find(nodeScale, 1)));
        one.threshold_value = threshold;
        one.n_available_primary_rows = nnz(validFull);
        one.n_available_event_rows = nnz(fullEvent);
        one.available_event_fraction = nnz(fullEvent) ./ max(nnz(validFull), 1);
        one.n_graph_nodes = nnz(validNode);
        one.n_graph_event_nodes = nnz(nodeEvent);
        one.graph_event_fraction = nnz(nodeEvent) ./ max(nnz(validNode), 1);
        one.selected_event_coverage_fraction = nnz(nodeEvent) ./ max(nnz(fullEvent), 1);
        one.mean_neighbor_event_fraction_all_nodes = mean(neighborFrac(validNode), 'omitnan');
        one.mean_neighbor_event_fraction_event_nodes = mean(neighborFrac(nodeEvent), 'omitnan');
        one.coverage_role = "postfit_audit_not_graph_input";
        one.labels_used_for_event_coverage = "none";
        one.arena_used_for_event_coverage = false;
        one.condition_used_for_event_coverage = false;
        coverageAudit = [coverageAudit; one]; %#ok<AGROW>
    end
    eventNodeAudit.(eventDefs.event_id(e)) = nodeFlagAll;
end
end

function fileName = local_active_event_summary_file(params)
if params.anchor_manifest_mode == "rare_enriched"
    fileName = string(params.expanded_event_summary_file);
else
    fileName = string(params.primary_event_summary_file);
end
end

function Audit = local_rare_composition_audit(Definition, Membership, nodeManifest, params)
Audit = Definition;
Audit.anchor_manifest_mode = repmat(string(params.anchor_manifest_mode), height(Audit), 1);
Audit.n_graph_defined_member_rows = Audit.n_member_nodes;
Audit.n_input_assigned_rows = zeros(height(Audit), 1);
Audit.n_assigned_base_rows = zeros(height(Audit), 1);
Audit.n_assigned_rare_enriched_rows = zeros(height(Audit), 1);
if ismember('rare_stratum_id', nodeManifest.Properties.VariableNames)
    assigned = string(nodeManifest.rare_stratum_id);
else
    assigned = repmat("none", height(nodeManifest), 1);
end
if ismember('anchor_stage', nodeManifest.Properties.VariableNames)
    stage = string(nodeManifest.anchor_stage);
else
    stage = repmat("base_time_even", height(nodeManifest), 1);
end
for i = 1:height(Audit)
    idx = double(nodeManifest.scale_index) == double(Audit.scale_index(i)) & ...
        assigned == string(Audit.rare_stratum_id(i));
    Audit.n_input_assigned_rows(i) = nnz(idx);
    Audit.n_assigned_base_rows(i) = nnz(idx & stage == "base_time_even");
    Audit.n_assigned_rare_enriched_rows(i) = nnz(idx & stage == "rare_strata_enriched");
end
Audit.n_membership_rows_total = repmat(height(Membership), height(Audit), 1);
Audit.labels_used_for_composition_audit = repmat("none", height(Audit), 1);
Audit.arena_used_for_composition_audit = false(height(Audit), 1);
Audit.condition_used_for_composition_audit = false(height(Audit), 1);
end

function Audit = local_rare_neighbor_retention_audit(Membership, nodeManifest, Graph, params)
Audit = table();
if isempty(Membership)
    return
end
n = height(nodeManifest);
keys = unique(Membership(:, {'scale_index','rare_stratum_id'}), 'rows', 'stable');
Edges = Graph.Edges;
for i = 1:height(keys)
    memberRows = double(Membership.scale_index) == double(keys.scale_index(i)) & ...
        string(Membership.rare_stratum_id) == string(keys.rare_stratum_id(i));
    member = false(n, 1);
    member(double(Membership.graph_node_id(memberRows))) = true;
    sourceEdge = member(double(Edges.source_node_id));
    targetRetained = member(double(Edges.target_node_id(sourceEdge)));
    scaleMask = double(nodeManifest.scale_index) == double(keys.scale_index(i));
    one = table();
    one.scale_index = double(keys.scale_index(i));
    one.chunk_sec = double(nodeManifest.chunk_sec(find(scaleMask, 1)));
    one.rare_stratum_id = string(keys.rare_stratum_id(i));
    one.anchor_manifest_mode = string(params.anchor_manifest_mode);
    one.n_member_nodes = nnz(member);
    one.n_outgoing_neighbor_edges = nnz(sourceEdge);
    one.n_retained_neighbor_edges = nnz(targetRetained);
    one.mean_neighbor_retention_fraction = mean(double(targetRetained), 'omitnan');
    one.background_fraction_within_scale = nnz(member) ./ max(nnz(scaleMask), 1);
    one.retention_over_background = one.mean_neighbor_retention_fraction ./ ...
        max(one.background_fraction_within_scale, eps);
    one.neighbor_input_rule = "postfit_graph_membership_join_not_edge_input";
    one.labels_used_for_neighbor_retention = "none";
    one.arena_used_for_neighbor_retention = false;
    one.condition_used_for_neighbor_retention = false;
    Audit = [Audit; one]; %#ok<AGROW>
end
end

function Audit = local_baseline_vs_enriched_coverage_audit(repoRoot, current, params)
Audit = current;
Audit.anchor_manifest_mode = repmat(string(params.anchor_manifest_mode), height(Audit), 1);
Audit.baseline_graph_root = repmat(resolve_repo_path(repoRoot, params.baseline_graph_input_dir), height(Audit), 1);
Audit.baseline_n_graph_nodes = nan(height(Audit), 1);
Audit.baseline_n_graph_event_nodes = nan(height(Audit), 1);
Audit.baseline_graph_event_fraction = nan(height(Audit), 1);
Audit.baseline_selected_event_coverage_fraction = nan(height(Audit), 1);
baselinePath = fullfile(resolve_repo_path(repoRoot, params.baseline_graph_input_dir), ...
    'graph_rare_event_coverage_audit.csv');
if isfile(baselinePath)
    Base = local_read_csv(baselinePath);
    keyCurrent = string(Audit.event_id) + "_" + string(round(double(Audit.scale_index)));
    keyBase = string(Base.event_id) + "_" + string(round(double(Base.scale_index)));
    [tf, loc] = ismember(keyCurrent, keyBase);
    Audit.baseline_n_graph_nodes(tf) = double(Base.n_graph_nodes(loc(tf)));
    Audit.baseline_n_graph_event_nodes(tf) = double(Base.n_graph_event_nodes(loc(tf)));
    Audit.baseline_graph_event_fraction(tf) = double(Base.graph_event_fraction(loc(tf)));
    Audit.baseline_selected_event_coverage_fraction(tf) = ...
        double(Base.selected_event_coverage_fraction(loc(tf)));
end
Audit.delta_graph_event_nodes_vs_baseline = Audit.n_graph_event_nodes - Audit.baseline_n_graph_event_nodes;
Audit.delta_graph_event_fraction_vs_baseline = Audit.graph_event_fraction - Audit.baseline_graph_event_fraction;
Audit.delta_selected_event_coverage_vs_baseline = Audit.selected_event_coverage_fraction - ...
    Audit.baseline_selected_event_coverage_fraction;
Audit.comparison_role = repmat("postfit_baseline_vs_current_audit_not_graph_input", height(Audit), 1);
Audit.labels_used_for_baseline_comparison = repmat("none", height(Audit), 1);
Audit.arena_used_for_baseline_comparison = false(height(Audit), 1);
Audit.condition_used_for_baseline_comparison = false(height(Audit), 1);
end

function defs = local_event_definitions()
defs = table();
defs.event_id = ["contact_present"; "contact_transition"; "mutual_approach"; ...
    "withdrawal"; "approach_withdraw_transition"; "strong_role_asymmetry"; ...
    "high_radial_speed"; "large_heading_turn"];
defs.source_column = ["contact_dwell_fraction"; "contact_transition_count"; ...
    "mutual_approach_dwell_fraction"; "withdrawal_dwell_fraction"; ...
    "approach_withdraw_transition_count"; "role_asymmetry_bias_mean"; ...
    "radial_speed_mean_mm_s"; "heading_large_turn_count"];
defs.event_rule = ["value_gt_zero"; "count_gt_zero"; "value_gt_zero"; ...
    "value_gt_zero"; "count_gt_zero"; "abs_value_ge_0p75"; ...
    "within_scale_ge_high_quantile"; "count_gt_zero"];
end

function [flag, threshold] = local_event_flag(T, rowMask, def, params, forcedThreshold)
if nargin < 5
    forcedThreshold = NaN;
end
col = string(def.source_column(1));
flag = false(height(T), 1);
threshold = NaN;
if ~ismember(col, T.Properties.VariableNames)
    return
end
valid = rowMask & isfinite(double(T.(col))) & T.event_valid_fraction >= params.min_event_valid_fraction;
x = double(T.(col));
switch string(def.event_rule(1))
    case {"value_gt_zero", "count_gt_zero"}
        threshold = 0;
        flag = valid & x > 0;
    case "abs_value_ge_0p75"
        threshold = 0.75;
        flag = valid & abs(x) >= threshold;
    case "within_scale_ge_high_quantile"
        if isfinite(forcedThreshold)
            threshold = forcedThreshold;
        elseif nnz(valid) >= 5
            threshold = quantile(x(valid), params.high_quantile_threshold);
        else
            threshold = Inf;
        end
        flag = valid & x >= threshold;
    otherwise
        error('run_condition_blind_motif_graph_build:UnknownEventRule', ...
            'Unknown event rule: %s', def.event_rule(1));
end
end

function frac = local_neighbor_event_fraction(Edges, eventFlag)
n = numel(eventFlag);
targetFlag = double(eventFlag(Edges.target_node_id));
frac = accumarray(Edges.source_node_id, targetFlag, [n 1], @mean, NaN);
end

function local_assert_event_summary_label_free(E)
if ismember('labels_used_for_event_summary', E.Properties.VariableNames)
    assert(all(string(E.labels_used_for_event_summary) == "none"), ...
        'run_condition_blind_motif_graph_build:EventSummaryLabelLeak', ...
        'Event summary must report labels_used_for_event_summary=none.');
end
if ismember('arena_used_for_event_summary', E.Properties.VariableNames)
    assert(~any(logical(E.arena_used_for_event_summary)), ...
        'run_condition_blind_motif_graph_build:EventSummaryArenaLeak', ...
        'Event summary must not use arena labels.');
end
if ismember('condition_used_for_event_summary', E.Properties.VariableNames)
    assert(~any(logical(E.condition_used_for_event_summary)), ...
        'run_condition_blind_motif_graph_build:EventSummaryConditionLeak', ...
        'Event summary must not use condition labels.');
end
end

function [sUndir, tUndir] = local_undirected_edge_pairs(Edges)
s = Edges.source_node_id;
t = Edges.target_node_id;
a = min(s, t);
b = max(s, t);
pairs = unique([a b], 'rows');
sUndir = pairs(:, 1);
tUndir = pairs(:, 2);
end

function keys = local_undirected_edge_keys(Edges, n)
[s, t] = local_undirected_edge_pairs(Edges);
keys = uint64(s) + uint64(n + 1) .* uint64(t);
end

function local_ensure_dir(pathText)
if ~exist(pathText, 'dir')
    mkdir(pathText);
end
end

function T = local_read_csv(pathText)
opts = detectImportOptions(pathText, 'FileType', 'text', ...
    'Delimiter', ',', 'TextType', 'string');
T = readtable(pathText, opts);
end
