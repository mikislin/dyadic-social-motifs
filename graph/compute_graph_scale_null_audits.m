function [Mixing, NullAudit] = compute_graph_scale_null_audits(Graph, nodeManifest, arenaLabel)
%COMPUTE_GRAPH_SCALE_NULL_AUDITS Persist scale mixing and null-normalized composition.
%
% Graph edges must already have been constructed from condition-blind numeric
% scores. Scale and session identifiers are used only after graph construction
% for audit summaries. Arena is optional and is post-fit audit-only.

assert(isstruct(Graph) && isfield(Graph, 'Edges'), ...
    'compute_graph_scale_null_audits:BadGraph', 'Graph.Edges is required.');
assert(istable(nodeManifest) && all(ismember( ...
    {'graph_node_id','scale_index','chunk_sec','session_index'}, ...
    nodeManifest.Properties.VariableNames)), ...
    'compute_graph_scale_null_audits:BadNodeManifest', ...
    'nodeManifest is missing required condition-blind provenance columns.');

n = height(nodeManifest);
if nargin < 3 || isempty(arenaLabel)
    arenaLabel = strings(n, 1);
else
    arenaLabel = string(arenaLabel(:));
end
assert(numel(arenaLabel) == n, ...
    'compute_graph_scale_null_audits:BadArenaLabels', ...
    'arenaLabel must have one value per graph node.');

Edges = Graph.Edges;
source = double(Edges.source_node_id);
target = double(Edges.target_node_id);
scale = double(nodeManifest.scale_index);
chunkSec = double(nodeManifest.chunk_sec);
session = string(nodeManifest.session_index);
scales = unique(scale, 'stable')';

Mixing = table();
for sourceScale = scales
    sourceNode = scale == sourceScale;
    sourceEdge = sourceNode(source);
    nSourceNodes = nnz(sourceNode);
    nSourceEdges = nnz(sourceEdge);
    for targetScale = scales
        targetNode = scale == targetScale;
        nTargetNodes = nnz(targetNode);
        nEdges = nnz(sourceEdge & scale(target) == targetScale);
        observed = nEdges ./ max(nSourceEdges, 1);
        if sourceScale == targetScale
            expected = max(nTargetNodes - 1, 0) ./ max(n - 1, 1);
        else
            expected = nTargetNodes ./ max(n - 1, 1);
        end
        one = table();
        one.source_scale_index = sourceScale;
        one.source_chunk_sec = chunkSec(find(sourceNode, 1));
        one.target_scale_index = targetScale;
        one.target_chunk_sec = chunkSec(find(targetNode, 1));
        one.n_source_nodes = nSourceNodes;
        one.n_target_nodes = nTargetNodes;
        one.n_directed_edges = nEdges;
        one.observed_source_edge_fraction = observed;
        one.random_mixing_expected_fraction = expected;
        one.observed_over_random_ratio = observed ./ max(expected, eps);
        one.is_same_scale_cell = sourceScale == targetScale;
        one.audit_role = "postfit_scale_mixing_not_graph_input";
        one.labels_used_for_scale_mixing = "none";
        one.arena_used_for_scale_mixing = false;
        one.condition_used_for_scale_mixing = false;
        Mixing = [Mixing; one]; %#ok<AGROW>
    end
end

NullAudit = table();
for sourceScale = [NaN scales]
    if isnan(sourceScale)
        nodeMask = true(n, 1);
        scope = "all_scales";
        scopeScale = NaN;
        scopeSec = NaN;
    else
        nodeMask = scale == sourceScale;
        scope = "source_scale";
        scopeScale = sourceScale;
        scopeSec = chunkSec(find(nodeMask, 1));
    end
    edgeMask = nodeMask(source);

    expectedSameScale = local_expected_same_group(scale, nodeMask);
    observedSameScale = mean(scale(source(edgeMask)) == scale(target(edgeMask)), 'omitnan');
    NullAudit = [NullAudit; local_null_row(scope, scopeScale, scopeSec, ...
        "same_scale", nnz(nodeMask), nnz(edgeMask), observedSameScale, ...
        expectedSameScale, "scale_index_postfit_audit", false)]; %#ok<AGROW>

    expectedSameSession = local_expected_same_group(session, nodeMask);
    observedSameSession = mean(session(source(edgeMask)) == session(target(edgeMask)), 'omitnan');
    NullAudit = [NullAudit; local_null_row(scope, scopeScale, scopeSec, ...
        "same_session", nnz(nodeMask), nnz(edgeMask), observedSameSession, ...
        expectedSameSession, "session_index_postfit_audit", false)]; %#ok<AGROW>

    validArenaNode = arenaLabel ~= "" & ~ismissing(arenaLabel);
    arenaSourceNode = nodeMask & validArenaNode;
    arenaEdge = edgeMask & validArenaNode(source) & validArenaNode(target);
    if any(arenaSourceNode) && any(arenaEdge)
        expectedSameArena = local_expected_same_group(arenaLabel(validArenaNode), ...
            nodeMask(validArenaNode));
        observedSameArena = mean(arenaLabel(source(arenaEdge)) == arenaLabel(target(arenaEdge)), 'omitnan');
    else
        expectedSameArena = NaN;
        observedSameArena = NaN;
    end
    NullAudit = [NullAudit; local_null_row(scope, scopeScale, scopeSec, ...
        "same_arena", nnz(arenaSourceNode), nnz(arenaEdge), observedSameArena, ...
        expectedSameArena, "arena_label_postfit_audit_only", true)]; %#ok<AGROW>
end
end

function expected = local_expected_same_group(group, sourceMask)
group = string(group(:));
sourceMask = logical(sourceMask(:));
n = numel(group);
if n <= 1 || ~any(sourceMask)
    expected = NaN;
    return
end
[~, ~, groupIndex] = unique(group, 'stable');
groupCount = accumarray(groupIndex, 1);
expectedByNode = (groupCount(groupIndex) - 1) ./ (n - 1);
expected = mean(expectedByNode(sourceMask), 'omitnan');
end

function one = local_null_row(scope, scaleIndex, chunkSec, compositionId, ...
    nNodes, nEdges, observed, expected, sourceRole, arenaUsed)
one = table();
one.source_scope = string(scope);
one.source_scale_index = scaleIndex;
one.source_chunk_sec = chunkSec;
one.composition_id = string(compositionId);
one.n_source_nodes = nNodes;
one.n_audited_edges = nEdges;
one.observed_neighbor_fraction = observed;
one.random_mixing_expected_fraction = expected;
one.observed_over_random_ratio = observed ./ max(expected, eps);
one.null_rule = "node_count_preserving_random_target_assignment_excluding_self";
one.composition_source_role = string(sourceRole);
one.audit_role = "postfit_null_normalized_neighbor_composition_not_graph_input";
one.labels_used_for_graph = "none";
one.arena_used_for_graph = false;
one.arena_used_for_postfit_audit = logical(arenaUsed);
one.condition_used_for_audit = false;
end
