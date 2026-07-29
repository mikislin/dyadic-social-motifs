%RUN_09_DISCOVER_CONDITION_BLIND_MOTIF_CANDIDATES Partition frozen run_08 graph.
%
% Parameters live in config/motif_candidate_discovery_config.csv.
% This script remains a thin entry point; reusable logic lives in clustering/.
%
% Smoke/default:
%   setenv('RUN09_CANDIDATE_RUN_MODE','');
%   setenv('RUN09_CANDIDATE_OUTPUT_DIR','');
%   run('paper/run_09_discover_condition_blind_motif_candidates.m')
%
% Full:
%   setenv('RUN09_CANDIDATE_RUN_MODE','full');
%   setenv('RUN09_CANDIDATE_OUTPUT_DIR','');
%   run('paper/run_09_discover_condition_blind_motif_candidates.m')

repoRoot = fileparts(fileparts(mfilename('fullpath')));
cd(repoRoot);
addpath(genpath(repoRoot));

configPath = fullfile(repoRoot, 'config', 'motif_candidate_discovery_config.csv');
run_condition_blind_motif_candidate_discovery(repoRoot, struct('configPath', configPath));
