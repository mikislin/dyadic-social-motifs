%RUN_06_BUILD_MULTISCALE_CHUNKS_AND_SELECT_SCALES Build condition-blind chunks.
%
% Parameters live in config/multiscale_chunking_config.csv.
% This script must remain an entry point only; reusable logic lives in chunks/.

% setenv('RUN06_CHUNK_RUN_MODE','full');
% setenv('RUN06_CHUNK_OUTPUT_DIR','derived/chunks_motif_discovery');
% run('run_06_build_multiscale_chunks_and_select_scales.m')

% setenv('RUN06_CHUNK_RUN_MODE','full');
% setenv('RUN06_CHUNK_OUTPUT_DIR','');
% setenv('RUN06_SAVE_CHUNK_MAT_ARTIFACT','false');
% setenv('RUN06_USE_SCALE_SUMMARY_SHARDS','true');
% setenv('RUN06_REUSE_SCALE_SUMMARY_SHARDS','true');
% run('run_06_build_multiscale_chunks_and_select_scales.m')


repoRoot = fileparts(fileparts(mfilename('fullpath')));
cd(repoRoot);
addpath(genpath(repoRoot));

configPath = fullfile(repoRoot, 'config', 'multiscale_chunking_config.csv');
run_multiscale_chunking_and_scale_selection(repoRoot, struct('configPath', configPath));
