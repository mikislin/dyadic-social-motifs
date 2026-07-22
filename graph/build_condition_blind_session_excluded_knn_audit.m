function Graph = build_condition_blind_session_excluded_knn_audit(X, rowMeta, params)
%BUILD_CONDITION_BLIND_SESSION_EXCLUDED_KNN_AUDIT Audit-only exact kNN graph.
%
% Distances use only numeric condition-blind embedding scores. Session identity
% is used solely to remove same-session candidates after distances are formed.
% This graph must never replace the primary run_08 graph or define motifs.

X = double(X);
assert(ismatrix(X) && size(X, 1) >= 3 && size(X, 2) >= 2, ...
    'build_condition_blind_session_excluded_knn_audit:BadInputMatrix', ...
    'X must be an n-by-d matrix with at least 3 rows and 2 columns.');
assert(all(isfinite(X), 'all'), ...
    'build_condition_blind_session_excluded_knn_audit:NonFiniteInput', ...
    'Audit input matrix contains non-finite values.');
assert(istable(rowMeta) && height(rowMeta) == size(X, 1) && ...
    ismember('session_index', rowMeta.Properties.VariableNames), ...
    'build_condition_blind_session_excluded_knn_audit:BadRowMeta', ...
    'rowMeta must contain session_index and have one row per X row.');

n = size(X, 1);
session = findgroups(string(rowMeta.session_index));
sessionCounts = local_group_counts(session);
k = min(floor(params.k_neighbors), n - max(sessionCounts));
assert(k >= 2, ...
    'build_condition_blind_session_excluded_knn_audit:TooFewCrossSessionCandidates', ...
    'At least two cross-session candidates are required per node.');
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
        d2(i, session == session(rows(i))) = Inf;
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
distanceScale = median(neighborDistance(neighborDistance > 0), 'omitnan');
if ~(isfinite(distanceScale) && distanceScale > 0)
    distanceScale = 1;
end

Edges = table();
Edges.source_node_id = sourceNode;
Edges.target_node_id = targetNode;
Edges.neighbor_rank = neighborRank;
Edges.neighbor_distance = neighborDistance;
Edges.edge_weight = exp(-neighborDistance ./ distanceScale);
Edges.is_mutual_neighbor = logical(isMutual);
Edges.source_embedding_row_id = local_table_value(rowMeta, 'embedding_row_id', sourceNode);
Edges.target_embedding_row_id = local_table_value(rowMeta, 'embedding_row_id', targetNode);
Edges.session_identity_used_for_audit_exclusion = true(nEdges, 1);
Edges.audit_only_not_primary_graph = true(nEdges, 1);
Edges.labels_used_for_graph_distance = repmat("none", nEdges, 1);
Edges.arena_used_for_graph_distance = false(nEdges, 1);
Edges.condition_used_for_graph_distance = false(nEdges, 1);

Graph = struct();
Graph.X = X;
Graph.Edges = Edges;
Graph.k = k;
Graph.n_nodes = n;
Graph.n_dims = size(X, 2);
Graph.distance_scale = distanceScale;
Graph.graph_role = "session_excluded_condition_blind_sensitivity_audit_only";
Graph.session_identity_used_for_exclusion = true;
Graph.labels_used_for_graph_distance = "none";
Graph.arena_used_for_graph_distance = false;
Graph.condition_used_for_graph_distance = false;
end

function counts = local_group_counts(group)
counts = accumarray(double(group), 1);
end

function values = local_table_value(T, varName, idx)
if ismember(varName, T.Properties.VariableNames)
    values = T.(varName)(idx);
else
    values = idx(:);
end
end
