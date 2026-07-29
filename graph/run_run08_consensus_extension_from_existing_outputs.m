function outputs = run_run08_consensus_extension_from_existing_outputs(repoRoot, opts)
%RUN_RUN08_CONSENSUS_EXTENSION_FROM_EXISTING_OUTPUTS Continue a completed run-08.
%
% This avoids rebuilding the already validated primary graph. It reads the
% exact numeric matrix and node manifest from the run-08 MAT artifact and adds
% only the consensus-neighborhood and run-09 handoff files.

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
assert(logical(params.consensus_enabled), ...
    'run_run08_consensus_extension_from_existing_outputs:ConsensusDisabled', ...
    'consensus_enabled must be true.');
outRoot = resolve_repo_path(repoRoot, params.output_dir);
modelPath = fullfile(outRoot, char(params.graph_model_mat_file));
eventPath = fullfile(outRoot, 'graph_event_node_audit.csv');
assert(isfile(modelPath) && isfile(eventPath), ...
    'run_run08_consensus_extension_from_existing_outputs:MissingRun08Output', ...
    'A completed run-08 model and event-node audit are required in %s.', outRoot);
S = load(modelPath, 'GraphModel');
M = S.GraphModel;
required = {'X','nodeManifest','Edges'};
assert(all(isfield(M, required)), ...
    'run_run08_consensus_extension_from_existing_outputs:IncompleteGraphModel', ...
    'The run-08 model lacks X, nodeManifest, or Edges.');
Primary = struct('Edges', M.Edges, 'k', params.k_neighbors, ...
    'n_nodes', height(M.nodeManifest), 'n_dims', size(M.X, 2));
Event = local_read(eventPath);
Audit = build_condition_blind_consensus_neighborhood( ...
    M.X, M.nodeManifest, Primary, Event, params, outRoot);
outputs = Audit;
outputs.output_root = string(outRoot);
end

function T = local_read(pathText)
opts = detectImportOptions(pathText, 'FileType', 'text', ...
    'Delimiter', ',', 'TextType', 'string');
T = readtable(pathText, opts);
end
