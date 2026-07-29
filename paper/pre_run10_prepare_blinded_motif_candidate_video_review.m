%PRE_RUN10_PREPARE_BLINDED_MOTIF_CANDIDATE_VIDEO_REVIEW Prepare review media.
%
% This unregistered pre-run_10 utility creates blinded validation material.
% It is not a final motif figure stage and does not alter run_09 outputs.
%
% Smoke:
%   setenv('PRE_RUN10_VIDEO_REVIEW_MODE', 'smoke');
%   run('paper/pre_run10_prepare_blinded_motif_candidate_video_review.m')
%
% Full:
%   setenv('PRE_RUN10_VIDEO_REVIEW_MODE', 'full');
%   run('paper/pre_run10_prepare_blinded_motif_candidate_video_review.m')

repoRoot = fileparts(fileparts(mfilename('fullpath')));
cd(repoRoot);
addpath(genpath(repoRoot));

configPath = fullfile(repoRoot, 'config', ...
    'motif_candidate_video_review_config.csv');
prepare_blinded_motif_candidate_video_review( ...
    repoRoot, struct('configPath', configPath));
