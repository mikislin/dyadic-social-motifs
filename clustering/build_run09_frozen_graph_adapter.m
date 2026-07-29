function Adapter = build_run09_frozen_graph_adapter(Validation)
%BUILD_RUN09_FROZEN_GRAPH_ADAPTER Create a strict contiguous weighted graph.
%
% Only graph_node_id, source_node_id, target_node_id, and
% consensus_edge_weight enter the adapted MATLAB graph. Scale provenance and
% all other run_08 audit columns are deliberately excluded.

assert(isstruct(Validation) && isfield(Validation, 'pass') && ...
    Validation.pass, ...
    'build_run09_frozen_graph_adapter:FreezeNotValidated', ...
    'A passing frozen-handoff validation is required before adaptation.');

Nodes = Validation.nodes(:, {'graph_node_id'});
Edges = Validation.edges(:, {'source_node_id', 'target_node_id', ...
    'consensus_edge_weight'});
nodeIds = double(Nodes.graph_node_id);
sourceIds = double(Edges.source_node_id);
targetIds = double(Edges.target_node_id);
weights = double(Edges.consensus_edge_weight);
[sourcePresent, sourceIndex] = ismember(sourceIds, nodeIds);
[targetPresent, targetIndex] = ismember(targetIds, nodeIds);
assert(all(sourcePresent) && all(targetPresent), ...
    'build_run09_frozen_graph_adapter:MissingEndpoint', ...
    'All endpoints must map to the validated frozen node set.');

n = numel(nodeIds);
sourceIndex = double(sourceIndex);
targetIndex = double(targetIndex);
matlabGraph = graph(sourceIndex, targetIndex, weights, n);

degree = accumarray([sourceIndex; targetIndex], 1, [n 1]);
weightedDegree = accumarray([sourceIndex; targetIndex], ...
    [weights; weights], [n 1]);
NodeIndex = table(nodeIds, (1:n)', degree, weightedDegree, ...
    'VariableNames', {'graph_node_id', 'matlab_node_index', ...
    'degree', 'weighted_degree'});

endNodes = double(matlabGraph.Edges.EndNodes);
roundSourceIds = nodeIds(endNodes(:, 1));
roundTargetIds = nodeIds(endNodes(:, 2));
roundPairs = sort([roundSourceIds roundTargetIds], 2);
originalPairs = sort([sourceIds targetIds], 2);
roundTrip = sortrows([roundPairs double(matlabGraph.Edges.Weight)], [1 2]);
original = sortrows([originalPairs weights], [1 2]);
roundTripPass = isequal(roundTrip, original);
weightPass = isequal(roundTrip(:, 3), original(:, 3));
edgeVariables = string(matlabGraph.Edges.Properties.VariableNames);
nodeVariables = string(matlabGraph.Nodes.Properties.VariableNames);
schemaPass = isequal(edgeVariables, ["EndNodes", "Weight"]) && ...
    isempty(nodeVariables);
degreePass = isequal(degree, Validation.degree);
weightedDegreePass = isequal(weightedDegree, Validation.weighted_degree);
weightedDensity = sum(weights) / (n * (n - 1) / 2);

Audit = table( ...
    ["n_nodes"; "n_edges"; "source_columns_used"; "target_columns_used"; ...
     "matlab_graph_edge_schema"; "matlab_graph_node_attributes"; ...
     "edge_roundtrip_exact"; "weight_roundtrip_exact"; ...
     "degree_recomputation_exact"; "weighted_degree_recomputation_exact"; ...
     "weighted_graph_density"], ...
    [string(n); string(numel(weights)); "graph_node_id"; ...
     "source_node_id;target_node_id;consensus_edge_weight"; ...
     strjoin(edgeVariables, ";"); ...
     i_empty_text(nodeVariables); string(roundTripPass); ...
     string(weightPass); string(degreePass); ...
     string(weightedDegreePass); compose('%.17g', weightedDensity)], ...
    [true; true; true; true; schemaPass; isempty(nodeVariables); ...
     roundTripPass; weightPass; degreePass; weightedDegreePass; ...
     isfinite(weightedDensity) && weightedDensity > 0], ...
    [ ...
     "Validated frozen node count."
     "Validated frozen edge count."
     "Only immutable node identity was mapped."
     "Only endpoints and consensus edge weight were passed."
     "No additional edge attributes reached the backend graph."
     "No node provenance or audit attributes reached the backend graph."
     "Node-ID edge pairs survive contiguous mapping and inverse mapping exactly."
     "Every IEEE-754 edge weight survives the adapter exactly."
     "Induced degree is recomputed from adapter endpoints."
     "Weighted degree is recomputed from adapter endpoints and weights."
     "CPM resolutions are derived from this condition-blind numeric density."
     ], ...
    'VariableNames', {'audit_item', 'observed_value', 'pass', 'details'});

assert(all(Audit.pass), ...
    'build_run09_frozen_graph_adapter:RoundTripMismatch', ...
    'The strict graph adapter failed an exact schema or round-trip gate.');

Adapter = struct();
Adapter.freeze_id = Validation.freeze_id;
Adapter.node_ids = nodeIds;
Adapter.source_index = sourceIndex;
Adapter.target_index = targetIndex;
Adapter.weights = weights;
Adapter.graph = matlabGraph;
Adapter.degree = degree;
Adapter.weighted_degree = weightedDegree;
Adapter.weighted_density = weightedDensity;
Adapter.total_edge_weight = sum(weights);
Adapter.node_index = NodeIndex;
Adapter.audit = Audit;
Adapter.backend_edge_columns = ...
    ["source_index", "target_index", "consensus_edge_weight"];
Adapter.backend_node_columns = "graph_node_id";
end

function value = i_empty_text(values)
if isempty(values)
    value = "none";
else
    value = strjoin(values, ";");
end
end
