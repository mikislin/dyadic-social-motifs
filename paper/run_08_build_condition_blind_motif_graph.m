%RUN_08_BUILD_CONDITION_BLIND_MOTIF_GRAPH Build condition-blind graph layer.
%
% Parameters live in config/multiscale_graph_config.csv.
% This script must remain an entry point only; reusable logic lives in graph/.
%
% Interactive examples:
%   Smoke/default:
%     setenv('RUN08_GRAPH_RUN_MODE',''); setenv('RUN08_GRAPH_OUTPUT_DIR','');
%     setenv('RUN08_EMBEDDING_INPUT_DIR','');
%     run('paper/run_08_build_condition_blind_motif_graph.m')
%
%   Full, after full run_07 has produced derived/embedding_motif_discovery:
%     setenv('RUN08_GRAPH_RUN_MODE','full');
%     setenv('RUN08_GRAPH_OUTPUT_DIR','');
%     setenv('RUN08_EMBEDDING_INPUT_DIR','');
%     run('paper/run_08_build_condition_blind_motif_graph.m')
%
%   Fast debug without writing figure files:
%     setenv('RUN08_FIGURE_EXPORT_PNG','false');
%     setenv('RUN08_FIGURE_EXPORT_PDF','false');

repoRoot = fileparts(fileparts(mfilename('fullpath')));
cd(repoRoot);
addpath(genpath(repoRoot));

configPath = fullfile(repoRoot, 'config', 'multiscale_graph_config.csv');
run_condition_blind_motif_graph_build(repoRoot, struct('configPath', configPath));
