function Ensemble = execute_run09_partition_ensemble(Adapter, Bridge, params, opts)
%EXECUTE_RUN09_PARTITION_ENSEMBLE Run or resume the weighted ensemble.
%
% Every completed replicate is stored in its own atomically replaced MAT
% checkpoint. A checkpoint is reused only after its frozen-input contract,
% metadata, canonical membership hash, dimensions, independently recomputed
% quality, and finite outputs pass validation.

if nargin < 4 || isempty(opts)
    opts = struct();
end
opts = i_apply_defaults(opts, params);
if ~isfolder(opts.checkpoint_dir)
    mkdir(opts.checkpoint_dir);
end

[Specs, cpmResolutions, modResolutions] = i_build_specs(Adapter, Bridge, params);
nReplicates = height(Specs);
nNodes = numel(Adapter.node_ids);
Memberships = zeros(nNodes, nReplicates, 'uint32');
ReplicateAudit = i_empty_audit(Specs, params);

for cursor = 1:nReplicates
    spec = Specs(cursor, :);
    fprintf(['run_09 %s replicate %d/%d | resolution %.12g | ' ...
        'seed %d\n'], upper(spec.objective), cursor, nReplicates, ...
        spec.resolution, spec.seed);
    contract = i_make_contract(spec, Adapter, Bridge, params);
    checkpointPath = fullfile(opts.checkpoint_dir, ...
        spec.replicate_id + ".mat");
    try
        if isfile(checkpointPath)
            [Result, checkpointHash] = i_load_checkpoint( ...
                checkpointPath, contract, Adapter);
            executionSource = "resumed_checkpoint";
            fprintf('  resumed validated checkpoint %s\n', ...
                spec.replicate_id);
        else
            Result = run_weighted_leiden_to_convergence( ...
                Adapter, Bridge, spec.objective, spec.resolution, ...
                spec.seed, params);
            checkpointHash = i_save_checkpoint_atomic( ...
                checkpointPath, contract, Result);
            executionSource = "computed";
            fprintf('  checkpointed %s\n', spec.replicate_id);
        end
        Memberships(:, cursor) = Result.membership;
        ReplicateAudit = i_record_success(ReplicateAudit, cursor, ...
            Result, checkpointPath, checkpointHash, ...
            contract.contract_sha256, executionSource);
    catch ME
        ReplicateAudit.status(cursor) = "error";
        ReplicateAudit.error_message(cursor) = string(getReport( ...
            ME, 'extended', 'hyperlinks', 'off'));
    end
    if opts.replicate_audit_path ~= ""
        write_run09_table_atomic(ReplicateAudit(1:cursor, :), ...
            opts.replicate_audit_path);
    end
end

Ensemble = struct();
Ensemble.memberships = Memberships;
Ensemble.node_ids = uint64(Adapter.node_ids);
Ensemble.replicate_audit = ReplicateAudit;
Ensemble.cpm_resolution_grid = cpmResolutions;
Ensemble.modularity_resolution_grid = modResolutions;
Ensemble.weighted_graph_density = Adapter.weighted_density;
Ensemble.all_replicates_converged = all(ReplicateAudit.converged) && ...
    all(ReplicateAudit.status == "success");
Ensemble.checkpoint_dir = string(opts.checkpoint_dir);
end

function opts = i_apply_defaults(opts, params)
if ~isfield(opts, 'checkpoint_dir') || ...
        strlength(string(opts.checkpoint_dir)) == 0
    opts.checkpoint_dir = fullfile(string(params.output_dir), ...
        string(params.checkpoint_dir));
end
if ~isfield(opts, 'replicate_audit_path') || ...
        strlength(string(opts.replicate_audit_path)) == 0
    opts.replicate_audit_path = "";
end
opts.checkpoint_dir = string(opts.checkpoint_dir);
opts.replicate_audit_path = string(opts.replicate_audit_path);
end

function [Specs, cpmResolutions, modResolutions] = ...
        i_build_specs(Adapter, Bridge, params)
cpmResolutions = Adapter.weighted_density .* ...
    params.active_cpm_density_multipliers(:)';
cpmMultipliers = params.active_cpm_density_multipliers(:)';
cpmSeeds = params.active_cpm_seeds(:)';
modResolutions = params.active_modularity_resolutions(:)';
modSeeds = params.active_modularity_seeds(:)';

replicateId = strings(0, 1);
analysisRole = strings(0, 1);
backendId = strings(0, 1);
objective = strings(0, 1);
resolutionRule = strings(0, 1);
resolutionIndex = zeros(0, 1);
densityMultiplier = zeros(0, 1);
resolution = zeros(0, 1);
seed = zeros(0, 1);
for r = 1:numel(cpmResolutions)
    for s = 1:numel(cpmSeeds)
        replicateId(end+1, 1) = compose( ... %#ok<AGROW>
            "run09_cpm_g%02d_seed%06d", r, cpmSeeds(s));
        analysisRole(end+1, 1) = "primary_weighted_leiden_cpm"; %#ok<AGROW>
        backendId(end+1, 1) = string(Bridge.backend_id); %#ok<AGROW>
        objective(end+1, 1) = "cpm"; %#ok<AGROW>
        resolutionRule(end+1, 1) = params.cpm_resolution_rule; %#ok<AGROW>
        resolutionIndex(end+1, 1) = r; %#ok<AGROW>
        densityMultiplier(end+1, 1) = cpmMultipliers(r); %#ok<AGROW>
        resolution(end+1, 1) = cpmResolutions(r); %#ok<AGROW>
        seed(end+1, 1) = cpmSeeds(s); %#ok<AGROW>
    end
end
if logical(params.enable_modularity_challenger)
    for r = 1:numel(modResolutions)
        for s = 1:numel(modSeeds)
            replicateId(end+1, 1) = compose( ... %#ok<AGROW>
                "run09_mod_g%02d_seed%06d", r, modSeeds(s));
            analysisRole(end+1, 1) = ... %#ok<AGROW>
                "audit_only_weighted_leiden_modularity";
            backendId(end+1, 1) = string(Bridge.backend_id); %#ok<AGROW>
            objective(end+1, 1) = "modularity"; %#ok<AGROW>
            resolutionRule(end+1, 1) = ... %#ok<AGROW>
                "canonical_weighted_modularity_gamma";
            resolutionIndex(end+1, 1) = r; %#ok<AGROW>
            densityMultiplier(end+1, 1) = NaN; %#ok<AGROW>
            resolution(end+1, 1) = modResolutions(r); %#ok<AGROW>
            seed(end+1, 1) = modSeeds(s); %#ok<AGROW>
        end
    end
end
membershipColumn = (1:numel(replicateId))';
Specs = table(replicateId, analysisRole, backendId, objective, ...
    resolutionRule, resolutionIndex, densityMultiplier, resolution, ...
    seed, membershipColumn, ...
    'VariableNames', {'replicate_id', 'analysis_role', 'backend_id', ...
    'objective', 'resolution_rule', 'resolution_index', ...
    'density_multiplier', 'resolution', 'seed', ...
    'membership_mat_column'});
end

function Audit = i_empty_audit(Specs, params)
n = height(Specs);
Audit = Specs;
Audit.native_resolution = nan(n, 1);
Audit.randomness = repmat(params.leiden_randomness, n, 1);
Audit.iteration_mode = repmat(params.iteration_mode, n, 1);
Audit.n_iterations = nan(n, 1);
Audit.converged = false(n, 1);
Audit.quality = nan(n, 1);
Audit.quality_definition = strings(n, 1);
Audit.runtime_seconds = nan(n, 1);
Audit.candidate_count = nan(n, 1);
Audit.membership_sha256 = strings(n, 1);
Audit.checkpoint_contract_sha256 = strings(n, 1);
Audit.checkpoint_file = strings(n, 1);
Audit.checkpoint_sha256 = strings(n, 1);
Audit.execution_source = repmat("not_started", n, 1);
Audit.status = repmat("not_started", n, 1);
Audit.error_message = strings(n, 1);
end

function contract = i_make_contract(spec, Adapter, Bridge, params)
contract = struct();
contract.schema = string(params.checkpoint_schema);
contract.replicate_id = spec.replicate_id;
contract.freeze_id = string(Adapter.freeze_id);
contract.node_input_sha256 = string(params.expected_node_hash);
contract.edge_input_sha256 = string(params.expected_edge_hash);
contract.backend_id = string(Bridge.backend_id);
contract.dispatcher_sha256 = string(params.expected_dispatcher_sha256);
contract.objective = spec.objective;
contract.resolution = double(spec.resolution);
contract.seed = double(spec.seed);
contract.randomness = double(params.leiden_randomness);
contract.iteration_mode = string(params.iteration_mode);
contract.maximum_iterations = double(params.maximum_iterations);
contract.n_nodes = double(numel(Adapter.node_ids));
contract.n_edges = double(numel(Adapter.weights));
payload = strjoin([ ...
    "schema=" + contract.schema
    "replicate_id=" + contract.replicate_id
    "freeze_id=" + contract.freeze_id
    "node_input_sha256=" + contract.node_input_sha256
    "edge_input_sha256=" + contract.edge_input_sha256
    "backend_id=" + contract.backend_id
    "dispatcher_sha256=" + contract.dispatcher_sha256
    "objective=" + contract.objective
    "resolution=" + compose('%.17g', contract.resolution)
    "seed=" + compose('%.0f', contract.seed)
    "randomness=" + compose('%.17g', contract.randomness)
    "iteration_mode=" + contract.iteration_mode
    "maximum_iterations=" + compose('%.0f', contract.maximum_iterations)
    "n_nodes=" + compose('%.0f', contract.n_nodes)
    "n_edges=" + compose('%.0f', contract.n_edges)], newline);
contract.contract_sha256 = compute_run09_sha256_text(payload);
end

function checkpointHash = i_save_checkpoint_atomic(pathText, contract, Result)
checkpoint = struct();
checkpoint.contract = contract;
checkpoint.membership = uint32(Result.membership);
checkpoint.native_resolution = double(Result.native_resolution);
checkpoint.n_iterations = double(Result.n_iterations);
checkpoint.converged = logical(Result.converged);
checkpoint.quality = double(Result.quality);
checkpoint.quality_definition = string(Result.quality_definition);
checkpoint.runtime_seconds = double(Result.runtime_seconds);
checkpoint.n_candidates = double(Result.n_candidates);
checkpoint.membership_sha256 = string(Result.membership_hash);
checkpoint.status = string(Result.status);
save_run09_mat_atomic(pathText, struct('checkpoint', checkpoint), false);
checkpointHash = compute_file_sha256(pathText);
end

function [Result, checkpointHash] = ...
        i_load_checkpoint(pathText, expectedContract, Adapter)
loaded = load(pathText, 'checkpoint');
assert(isfield(loaded, 'checkpoint') && isstruct(loaded.checkpoint), ...
    'execute_run09_partition_ensemble:BadCheckpoint', ...
    'Checkpoint lacks its required checkpoint struct: %s', pathText);
C = loaded.checkpoint;
assert(isfield(C, 'contract') && isstruct(C.contract) && ...
    isfield(C.contract, 'contract_sha256') && ...
    string(C.contract.contract_sha256) == ...
    string(expectedContract.contract_sha256) && ...
    isequaln(orderfields(C.contract), orderfields(expectedContract)), ...
    'execute_run09_partition_ensemble:CheckpointContractMismatch', ...
    'Checkpoint contract mismatch: %s', pathText);
required = {'membership','native_resolution','n_iterations','converged', ...
    'quality','quality_definition','runtime_seconds','n_candidates', ...
    'membership_sha256','status'};
assert(all(isfield(C, required)), ...
    'execute_run09_partition_ensemble:BadCheckpoint', ...
    'Checkpoint fields are incomplete: %s', pathText);
membership = uint32(C.membership(:));
assert(numel(membership) == numel(Adapter.node_ids) && ...
    isequal(membership, canonicalize_run09_membership(membership)), ...
    'execute_run09_partition_ensemble:BadCheckpointMembership', ...
    'Checkpoint membership is not canonical or has the wrong length.');
observedMembershipHash = compute_run09_membership_sha256(membership);
assert(observedMembershipHash == string(C.membership_sha256), ...
    'execute_run09_partition_ensemble:CheckpointMembershipHashMismatch', ...
    'Checkpoint membership SHA-256 mismatch: %s', pathText);
[observedQuality, observedDefinition] = ...
    compute_run09_partition_quality(Adapter, membership, ...
    expectedContract.objective, expectedContract.resolution);
tolerance = 1e-12 * max(1, abs(observedQuality));
assert(isfinite(C.quality) && ...
    abs(double(C.quality) - observedQuality) <= tolerance && ...
    string(C.quality_definition) == observedDefinition, ...
    'execute_run09_partition_ensemble:CheckpointQualityMismatch', ...
    'Checkpoint quality failed independent recomputation: %s', pathText);
assert(isfinite(C.native_resolution) && isfinite(C.n_iterations) && ...
    C.n_iterations >= 1 && isfinite(C.runtime_seconds) && ...
    C.runtime_seconds >= 0 && isfinite(C.n_candidates) && ...
    C.n_candidates == numel(unique(membership)) && ...
    any(string(C.status) == [ ...
    "success", "maximum_iterations_without_unchanged_membership"]), ...
    'execute_run09_partition_ensemble:BadCheckpointMetrics', ...
    'Checkpoint contains invalid or nonfinite metrics: %s', pathText);

Result = struct('membership', membership, ...
    'native_resolution', double(C.native_resolution), ...
    'n_iterations', double(C.n_iterations), ...
    'converged', logical(C.converged), ...
    'quality', observedQuality, ...
    'quality_definition', observedDefinition, ...
    'runtime_seconds', double(C.runtime_seconds), ...
    'n_candidates', double(C.n_candidates), ...
    'membership_hash', observedMembershipHash, ...
    'status', string(C.status));
checkpointHash = compute_file_sha256(pathText);
end

function Audit = i_record_success(Audit, row, Result, checkpointPath, ...
        checkpointHash, contractHash, executionSource)
Audit.native_resolution(row) = Result.native_resolution;
Audit.n_iterations(row) = Result.n_iterations;
Audit.converged(row) = Result.converged;
Audit.quality(row) = Result.quality;
Audit.quality_definition(row) = Result.quality_definition;
Audit.runtime_seconds(row) = Result.runtime_seconds;
Audit.candidate_count(row) = Result.n_candidates;
Audit.membership_sha256(row) = Result.membership_hash;
Audit.checkpoint_contract_sha256(row) = contractHash;
Audit.checkpoint_file(row) = string(checkpointPath);
Audit.checkpoint_sha256(row) = checkpointHash;
Audit.execution_source(row) = executionSource;
Audit.status(row) = Result.status;
Audit.error_message(row) = "";
end
