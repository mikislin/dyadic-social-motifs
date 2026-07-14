%RUN_07_BUILD_CONDITION_BLIND_EMBEDDING Build condition-blind embedding layer.
%
% Parameters live in config/multiscale_embedding_config.csv.
% This script must remain an entry point only; reusable logic lives in embedding/.
%
% Interactive examples:
%   Smoke/default:
%     setenv('RUN07_EMBEDDING_RUN_MODE',''); setenv('RUN07_EMBEDDING_OUTPUT_DIR','');
%     run('paper/run_07_build_condition_blind_embedding.m')
%
%   Full:
%     setenv('RUN07_EMBEDDING_RUN_MODE','full');
%     setenv('RUN07_EMBEDDING_OUTPUT_DIR','');
%     run('paper/run_07_build_condition_blind_embedding.m')
%
%   Fast debug without writing figure files:
%     setenv('RUN07_FIGURE_EXPORT_PNG','false');
%     setenv('RUN07_FIGURE_EXPORT_PDF','false');

repoRoot = fileparts(fileparts(mfilename('fullpath')));
cd(repoRoot);
addpath(genpath(repoRoot));

configPath = fullfile(repoRoot, 'config', 'multiscale_embedding_config.csv');
run_condition_blind_embedding_build(repoRoot, struct('configPath', configPath));
