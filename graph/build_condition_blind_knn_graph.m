function Graph = build_condition_blind_knn_graph(X, rowMeta, params)
%BUILD_CONDITION_BLIND_KNN_GRAPH Exact directed kNN graph from numeric scores.
%
% X must already contain only condition-blind numeric embedding columns.
% rowMeta is carried for identifiers only and is not used in distance
% calculations.

X = double(X);
assert(ismatrix(X) && size(X, 1) >= 3 && size(X, 2) >= 2, ...
    'build_condition_blind_knn_graph:BadInputMatrix', ...
    'X must be an n-by-d matrix with at least 3 rows and 2 columns.');
assert(all(isfinite(X), 'all'), ...
    'build_condition_blind_knn_graph:NonFiniteInput', ...
    'Graph input matrix contains non-finite values.');
assert(istable(rowMeta) && height(rowMeta) == size(X, 1), ...
    'build_condition_blind_knn_graph:BadRowMeta', ...
    'rowMeta must be a table with one row per X row.');

n = size(X, 1);
k = min(floor(params.k_neighbors), n - 1);
blockSize = min(max(floor(params.knn_block_size), 16), n);

nEdges = n * k;
sourceNode = zeros(nEdges, 1);
targetNode = zeros(nEdges, 1);
neighborRank = zeros(nEdges, 1);
neighborDistance = zeros(nEdges, 1);

rowNorm = sum(X.^2, 2)';
cursor = 1;
for startRow = 1:blockSize:n
    stopRow = min(n, startRow + blockSize - 1);
    rows = startRow:stopRow;
    Xb = X(rows, :);
    d2 = sum(Xb.^2, 2) + rowNorm - 2 .* (Xb * X');
    d2 = max(d2, 0);
    for i = 1:numel(rows)
        d2(i, rows(i)) = Inf;
    end
    [dSorted, idxSorted] = sort(d2, 2, 'ascend');
    idxSorted = idxSorted(:, 1:k);
    dSorted = sqrt(dSorted(:, 1:k));
    m = numel(rows) * k;
    edgeRows = cursor:(cursor + m - 1);
    sourceNode(edgeRows) = repelem(rows(:), k);
    targetNode(edgeRows) = reshape(idxSorted', [], 1);
    neighborRank(edgeRows) = repmat((1:k)', numel(rows), 1);
    neighborDistance(edgeRows) = reshape(dSorted', [], 1);
    cursor = cursor + m;
end

A = sparse(sourceNode, targetNode, true, n, n);
isMutual = full(A(sub2ind([n n], targetNode, sourceNode)));
distanceScale = median(neighborDistance(isfinite(neighborDistance) & neighborDistance > 0), 'omitnan');
if ~(isfinite(distanceScale) && distanceScale > 0)
    distanceScale = 1;
end
edgeWeight = exp(-neighborDistance ./ distanceScale);

sourceEmbeddingRow = local_table_value(rowMeta, 'embedding_row_id', sourceNode);
targetEmbeddingRow = local_table_value(rowMeta, 'embedding_row_id', targetNode);

Edges = table();
Edges.source_node_id = sourceNode;
Edges.target_node_id = targetNode;
Edges.neighbor_rank = neighborRank;
Edges.neighbor_distance = neighborDistance;
Edges.edge_weight = edgeWeight;
Edges.is_mutual_neighbor = logical(isMutual);
Edges.source_embedding_row_id = sourceEmbeddingRow;
Edges.target_embedding_row_id = targetEmbeddingRow;
Edges.labels_used_for_graph = repmat("none", nEdges, 1);
Edges.arena_used_for_graph = false(nEdges, 1);
Edges.condition_used_for_graph = false(nEdges, 1);

Graph = struct();
Graph.X = X;
Graph.Edges = Edges;
Graph.k = k;
Graph.n_nodes = n;
Graph.n_dims = size(X, 2);
Graph.distance_scale = distanceScale;
Graph.labels_used_for_graph = "none";
Graph.arena_used_for_graph = false;
Graph.condition_used_for_graph = false;
end

function values = local_table_value(T, varName, idx)
if ismember(varName, T.Properties.VariableNames)
    values = T.(varName)(idx);
else
    values = idx(:);
end
end
