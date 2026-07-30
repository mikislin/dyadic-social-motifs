%RUN_10_VALIDATE_MOTIF_CANDIDATES_BEHAVIORALLY Validate frozen candidates.
%
% Parameters and the frozen status rule live in:
%   config/motif_candidate_behavioral_validation_config.csv
%
% Smoke/default:
%   setenv('RUN10_VALIDATION_RUN_MODE','');
%   setenv('RUN10_VALIDATION_PHASE','');
%   setenv('RUN10_VALIDATION_OUTPUT_DIR','');
%   run('paper/run_10_validate_motif_candidates_behaviorally.m')
%
% Full:
%   setenv('RUN10_VALIDATION_RUN_MODE','full');
%   setenv('RUN10_VALIDATION_PHASE','full');
%   setenv('RUN10_VALIDATION_OUTPUT_DIR','');
%   run('paper/run_10_validate_motif_candidates_behaviorally.m')

repoRoot = fileparts(fileparts(mfilename('fullpath')));
cd(repoRoot);
addpath(genpath(repoRoot));

configPath = fullfile(repoRoot, 'config', ...
    'motif_candidate_behavioral_validation_config.csv');
run_condition_blind_candidate_behavioral_validation( ...
    repoRoot, struct('configPath', configPath));
