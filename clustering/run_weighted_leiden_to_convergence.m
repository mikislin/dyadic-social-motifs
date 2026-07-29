function Result = run_weighted_leiden_to_convergence( ...
        Adapter, Bridge, objective, resolution, seed, params, initialMembership)
%RUN_WEIGHTED_LEIDEN_TO_CONVERGENCE Run one deterministic weighted replicate.
%
% The native igraph negative-iteration mode can loop indefinitely on some
% graphs. This adapter implements the same documented stopping rule safely:
% call one native Leiden iteration at a time and stop after an iteration does
% not change the canonical membership, with a configured fail-closed cap.

assert(Bridge.ready && string(which('mexIgraphDispatcher')) == ...
    string(Bridge.bridge_mex), ...
    'run_weighted_leiden_to_convergence:BridgeNotReady', ...
    'The pinned weighted native bridge is not active.');
assert(isequal(string(Adapter.graph.Edges.Properties.VariableNames), ...
    ["EndNodes", "Weight"]), ...
    'run_weighted_leiden_to_convergence:BadGraphSchema', ...
    'Backend graph must contain only EndNodes and Weight edge fields.');
assert(isempty(Adapter.graph.Nodes.Properties.VariableNames), ...
    'run_weighted_leiden_to_convergence:NodeMetadataLeak', ...
    'Backend graph must not contain node attributes.');

objective = lower(string(objective));
resolution = double(resolution);
seed = double(seed);
switch objective
    case "cpm"
        nativeResolution = resolution;
        metric = 'cpm';
    case "modularity"
        nativeResolution = resolution / (2 * Adapter.total_edge_weight);
        metric = 'modularity';
    otherwise
        error('run_weighted_leiden_to_convergence:BadObjective', ...
            'Objective must be cpm or modularity.');
end

graphOpts = struct('isdirected', false, 'isweighted', true, ...
    'weight', 'Weight');
if nargin < 7 || isempty(initialMembership)
    current = (1:numel(Adapter.node_ids))';
    initializationRule = "singleton_initialization";
else
    current = double(canonicalize_run09_membership(initialMembership));
    assert(numel(current) == numel(Adapter.node_ids), ...
        'run_weighted_leiden_to_convergence:BadInitialization', ...
        'Initial membership length must match the frozen graph.');
    initializationRule = "provided_condition_blind_consensus";
end
converged = false;
iterations = 0;
igraph.rng(seed, 'mt19937');
timer = tic;
for iteration = 1:params.maximum_iterations
    methodOpts = struct('resolution', nativeResolution, ...
        'randomness', params.leiden_randomness, ...
        'nIterations', 1, 'metric', metric, ...
        'initial', current');
    next = mexIgraphDispatcher('cluster', Adapter.graph, 'leiden', ...
        graphOpts, methodOpts);
    next = canonicalize_run09_membership(next(:));
    iterations = iteration;
    if isequal(next, current)
        converged = true;
        current = next;
        break
    end
    current = next;
end
runtimeSeconds = toc(timer);

[quality, qualityDefinition] = compute_run09_partition_quality( ...
    Adapter, current, objective, resolution);
membershipHash = compute_run09_membership_sha256(current);

Result = struct();
Result.membership = uint32(current);
Result.objective = objective;
Result.resolution = resolution;
Result.native_resolution = nativeResolution;
Result.seed = seed;
Result.randomness = params.leiden_randomness;
Result.iteration_mode = params.iteration_mode;
Result.n_iterations = iterations;
Result.converged = converged;
Result.quality = quality;
Result.quality_definition = qualityDefinition;
Result.runtime_seconds = runtimeSeconds;
Result.n_candidates = numel(unique(current));
Result.membership_hash = membershipHash;
Result.status = i_status(converged);
Result.initialization_rule = initializationRule;
end

function value = i_status(converged)
if converged
    value = "success";
else
    value = "maximum_iterations_without_unchanged_membership";
end
end
