function [quality, definition] = compute_run09_partition_quality(Adapter, membership, objective, resolution)
%COMPUTE_RUN09_PARTITION_QUALITY Independently evaluate weighted objectives.

membership = canonicalize_run09_membership(membership);
s = Adapter.source_index;
t = Adapter.target_index;
w = Adapter.weights;
objective = lower(string(objective));
resolution = double(resolution);

switch objective
    case "cpm"
        internalWeight = sum(w(membership(s) == membership(t)));
        sizes = accumarray(membership, 1);
        rawQuality = internalWeight - resolution * ...
            sum(sizes .* (sizes - 1) / 2);
        quality = rawQuality / Adapter.total_edge_weight;
        definition = ...
            "(internal_weight-gamma*sum_choose2_candidate_size)/total_edge_weight";
    case "modularity"
        totalWeight = Adapter.total_edge_weight;
        internalMask = membership(s) == membership(t);
        internalByCandidate = accumarray(membership(s(internalMask)), ...
            w(internalMask), [max(membership) 1]);
        strengthByCandidate = accumarray(membership, ...
            Adapter.weighted_degree, [max(membership) 1]);
        quality = sum(internalByCandidate / totalWeight - ...
            resolution * (strengthByCandidate / (2 * totalWeight)).^2);
        definition = ...
            "weighted_generalized_modularity_with_canonical_gamma";
    otherwise
        error('compute_run09_partition_quality:BadObjective', ...
            'Unsupported objective: %s', objective);
end
assert(isfinite(quality), ...
    'compute_run09_partition_quality:NonfiniteQuality', ...
    'Independent partition quality must be finite.');
end
