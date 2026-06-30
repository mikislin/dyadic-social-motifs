%RUN_05_EXTRACT_MOTIF_DYAD_FEATURES Paper entry point for dyad features.
%
% Parameters live in config/motif_dyad_feature_extraction_config.csv.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
cd(repoRoot);
addpath(genpath(repoRoot));

configPath = fullfile(repoRoot, 'config', 'motif_dyad_feature_extraction_config.csv');
run_motif_dyad_feature_extraction(repoRoot, struct('configPath', configPath));
