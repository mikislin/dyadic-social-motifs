function Audit = build_run08_embedding_visualization_audits( ...
    embeddingRoot, outRoot, X, nodeManifest, Graph, params)
%BUILD_RUN08_EMBEDDING_VISUALIZATION_AUDITS Persist PCA variance and UMAP audits.
%
% UMAP is a deterministic visualization-only transform of numeric graph input
% scores. It is not used to construct graph edges, define rare strata, or
% define motifs or motif families.

embeddingRoot = string(embeddingRoot);
outRoot = string(outRoot);
modelPath = fullfile(embeddingRoot, char(params.embedding_model_mat_file));
assert(isfile(modelPath), ...
    'build_run08_embedding_visualization_audits:MissingRun07Model', ...
    'Missing run_07 embedding model: %s', modelPath);
S = load(modelPath, 'Embedding', 'Audit');
assert(isfield(S, 'Embedding') && isfield(S.Embedding, 'global_explained'), ...
    'build_run08_embedding_visualization_audits:MissingGlobalVariance', ...
    'The run_07 model does not contain Embedding.global_explained.');

explained = double(S.Embedding.global_explained(:));
Variance = table();
Variance.global_pc_index = (1:numel(explained))';
Variance.score_name = "global_pc" + compose('%02d', Variance.global_pc_index);
Variance.explained_variance_percent = explained;
Variance.cumulative_explained_variance_percent = cumsum(explained);
Variance.selected_for_run07_export = Variance.global_pc_index <= numel(S.Embedding.global_score_names);
Variance.selected_for_run08_graph = Variance.global_pc_index <= params.graph_n_pcs;
Variance.run07_embedding_model_path = repmat(modelPath, height(Variance), 1);
Variance.global_matrix_mode = repmat(local_global_matrix_mode(S), height(Variance), 1);
Variance.variance_role = repmat("run07_global_pca_variance_persisted_by_run08_audit", height(Variance), 1);
Variance.labels_used_for_global_pca = repmat("none", height(Variance), 1);
Variance.arena_used_for_global_pca = false(height(Variance), 1);
Variance.condition_used_for_global_pca = false(height(Variance), 1);
writetable(Variance, fullfile(outRoot, 'graph_global_pca_cumulative_variance_audit.csv'));

if ~logical(params.umap_enabled)
    Umap = table();
    UmapStatus = table("umap", "disabled_by_config", params.audit_random_seed, ...
        'VariableNames', {'audit_name','status','random_seed'});
    UmapStatus.audit_role = "visualization_only_not_graph_or_motif_input";
    UmapStatus.labels_used_for_umap = "none";
    UmapStatus.arena_used_for_umap = false;
    UmapStatus.condition_used_for_umap = false;
    writetable(UmapStatus, fullfile(outRoot, 'graph_umap_status_audit.csv'));
    Audit = struct('globalVariance', Variance, 'umapCoordinates', Umap, ...
        'umapStatus', UmapStatus);
    return
end

n = size(X, 1);
assert(n == height(nodeManifest), ...
    'build_run08_embedding_visualization_audits:RowMismatch', ...
    'X and nodeManifest must have identical row counts.');
assert(exist('umap', 'file') == 2, ...
    'build_run08_embedding_visualization_audits:UmapUnavailable', ...
    'MATLAB umap is required for the requested visualization audit.');

if n > params.umap_max_nodes
    sampleMask = local_scale_balanced_sample(nodeManifest, params.umap_max_nodes, params.audit_random_seed);
    Xfit = X(sampleMask, :);
    Meta = nodeManifest(sampleMask, :);
    fitParams = params;
    Gfit = build_condition_blind_knn_graph(Xfit, Meta, fitParams);
    inputRule = "fixed_seed_scale_balanced_visualization_sample";
else
    Xfit = X;
    Meta = nodeManifest;
    Gfit = Graph;
    inputRule = "all_graph_nodes";
end

neighborIndices = local_neighbor_matrix(Gfit.Edges, height(Meta), Gfit.k);
rng(params.umap_random_seed, 'twister');
Y2 = umap(Xfit, NumDimensions=2, NumNeighbors=Gfit.k, ...
    NeighborIndices=neighborIndices, NumEpochs=params.umap_num_epochs, ...
    Initialization="pca", Reproducible="on", Standardize=false, ...
    Distance="euclidean", EmbeddingDensity=params.umap_embedding_density);
rng(params.umap_random_seed, 'twister');
Y3 = umap(Xfit, NumDimensions=3, NumNeighbors=Gfit.k, ...
    NeighborIndices=neighborIndices, NumEpochs=params.umap_num_epochs, ...
    Initialization="pca", Reproducible="on", Standardize=false, ...
    Distance="euclidean", EmbeddingDensity=params.umap_embedding_density);
assert(all(isfinite(Y2), 'all') && all(isfinite(Y3), 'all'), ...
    'build_run08_embedding_visualization_audits:NonFiniteUmap', ...
    'UMAP produced non-finite coordinates.');

Umap = table();
Umap.umap_row_id = (1:height(Meta))';
Umap.graph_node_id = double(Meta.graph_node_id);
Umap.embedding_row_id = Meta.embedding_row_id;
Umap.scale_index = double(Meta.scale_index);
Umap.chunk_sec = double(Meta.chunk_sec);
if ismember('anchor_stage', Meta.Properties.VariableNames)
    Umap.anchor_stage = string(Meta.anchor_stage);
else
    Umap.anchor_stage = repmat("base_time_even", height(Meta), 1);
end
Umap.umap2_x = Y2(:, 1);
Umap.umap2_y = Y2(:, 2);
Umap.umap3_x = Y3(:, 1);
Umap.umap3_y = Y3(:, 2);
Umap.umap3_z = Y3(:, 3);
Umap.input_graph_dimensions = repmat(size(Xfit, 2), height(Meta), 1);
Umap.umap_num_neighbors = repmat(Gfit.k, height(Meta), 1);
Umap.umap_num_epochs = repmat(params.umap_num_epochs, height(Meta), 1);
Umap.umap_random_seed = repmat(params.umap_random_seed, height(Meta), 1);
Umap.umap_input_rule = repmat(inputRule, height(Meta), 1);
Umap.umap_role = repmat("visualization_only_not_graph_or_motif_input", height(Meta), 1);
Umap.rare_strata_used_for_umap = false(height(Meta), 1);
Umap.labels_used_for_umap = repmat("none", height(Meta), 1);
Umap.arena_used_for_umap = false(height(Meta), 1);
Umap.condition_used_for_umap = false(height(Meta), 1);
writetable(Umap, fullfile(outRoot, 'graph_umap_embedding_audit.csv'));

UmapStatus = table();
UmapStatus.audit_name = "umap_2d_and_3d";
UmapStatus.status = "completed";
UmapStatus.n_source_graph_nodes = n;
UmapStatus.n_umap_rows = height(Umap);
UmapStatus.n_input_dimensions = size(Xfit, 2);
UmapStatus.num_neighbors = Gfit.k;
UmapStatus.num_epochs = params.umap_num_epochs;
UmapStatus.random_seed = params.umap_random_seed;
UmapStatus.input_rule = inputRule;
UmapStatus.neighbor_indices_source = "run08_numeric_pc_knn_graph";
UmapStatus.audit_role = "visualization_only_not_graph_or_motif_input";
UmapStatus.labels_used_for_umap = "none";
UmapStatus.arena_used_for_umap = false;
UmapStatus.condition_used_for_umap = false;
writetable(UmapStatus, fullfile(outRoot, 'graph_umap_status_audit.csv'));

Audit = struct();
Audit.globalVariance = Variance;
Audit.umapCoordinates = Umap;
Audit.umapStatus = UmapStatus;
end

function mode = local_global_matrix_mode(S)
mode = "unknown";
if isfield(S, 'Audit') && isfield(S.Audit, 'global') && ...
        istable(S.Audit.global) && ismember('global_matrix_mode', S.Audit.global.Properties.VariableNames)
    mode = string(S.Audit.global.global_matrix_mode(1));
end
end

function matrix = local_neighbor_matrix(Edges, n, k)
assert(height(Edges) == n * k, ...
    'build_run08_embedding_visualization_audits:BadNeighborShape', ...
    'UMAP neighbor source must contain exactly n*k ordered edges.');
matrix = reshape(double(Edges.target_node_id), k, n)';
end

function mask = local_scale_balanced_sample(Meta, maxNodes, seed)
rng(seed, 'twister');
n = height(Meta);
mask = false(n, 1);
scales = unique(double(Meta.scale_index), 'stable')';
perScale = max(1, floor(maxNodes ./ numel(scales)));
for s = scales
    pool = find(double(Meta.scale_index) == s);
    take = min(perScale, numel(pool));
    if take == numel(pool)
        chosen = pool;
    else
        chosen = pool(sort(randperm(numel(pool), take)));
    end
    mask(chosen) = true;
end
remaining = min(maxNodes - nnz(mask), nnz(~mask));
if remaining > 0
    pool = find(~mask);
    mask(pool(sort(randperm(numel(pool), remaining)))) = true;
end
end
