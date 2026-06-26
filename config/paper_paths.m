function paths = paper_paths(varargin)
%PAPER_PATHS Resolve local paths
%
% Raw data are intentionally not committed. Point the pipeline at a local
% data directory with one of these environment variables:
%
%   DYADIC_SOCIAL_MOTIFS_DATA_DIR      directory containing ALLPAIRS-Feb2025.mat and files2run/
%   DYADIC_SOCIAL_MOTIFS_DBASE_PATH    full path to ALLPAIRS-Feb2025.mat
%   DYADIC_SOCIAL_MOTIFS_FILES2RUN_DIR full path to files2run/
%
% Users can also create config/paper_paths_local.m, which should
% accept and return a paths struct.

p = inputParser;
p.addParameter('requireRawData', false, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});
P = p.Results;

repoRoot = fileparts(fileparts(mfilename('fullpath')));

paths = struct();
paths.repoRoot = repoRoot;
paths.workspaceRoot = fileparts(repoRoot);
paths.configDir = fullfile(repoRoot, 'config');
paths.dataDir = fullfile(repoRoot, 'data');
paths.derivedDir = fullfile(repoRoot, 'derived');
paths.manifestPath = fullfile(paths.dataDir, 'session_manifest.csv');
paths.preprocessingPipelineConfigPath = fullfile(paths.configDir, 'preprocessing_pipeline_config.csv');
paths.experimentalGroupStylesPath = fullfile(paths.configDir, 'experimental_group_styles.csv');
paths.rawDataRoot = "";
paths.dbasePath = "";
paths.files2runDir = "";
paths.rawDataConfigured = false;
paths.rawDataInstructions = ...
    "Set DYADIC_SOCIAL_MOTIFS_DATA_DIR to the directory containing " + ...
    "ALLPAIRS-Feb2025.mat and files2run, or set " + ...
    "DYADIC_SOCIAL_MOTIFS_DBASE_PATH and DYADIC_SOCIAL_MOTIFS_FILES2RUN_DIR.";

dataRoot = string(strtrim(getenv('DYADIC_SOCIAL_MOTIFS_DATA_DIR')));
if dataRoot == ""
    dataRoot = string(strtrim(getenv('DYADIC_SOCIAL_MOTIFS_RAW_ROOT')));
end
dbasePath = string(strtrim(getenv('DYADIC_SOCIAL_MOTIFS_DBASE_PATH')));
files2runDir = string(strtrim(getenv('DYADIC_SOCIAL_MOTIFS_FILES2RUN_DIR')));

if dataRoot ~= ""
    paths.rawDataRoot = dataRoot;
    if dbasePath == ""
        dbasePath = fullfile(dataRoot, 'ALLPAIRS-Feb2025.mat');
    end
    if files2runDir == ""
        files2runDir = fullfile(dataRoot, 'files2run');
    end
end

if dbasePath ~= ""
    paths.dbasePath = dbasePath;
    if paths.rawDataRoot == ""
        paths.rawDataRoot = string(fileparts(char(dbasePath)));
    end
end
if files2runDir ~= ""
    paths.files2runDir = files2runDir;
    if paths.rawDataRoot == ""
        paths.rawDataRoot = string(fileparts(char(files2runDir)));
    end
end

if exist('paper_paths_local', 'file') == 2
    paths = paper_paths_local(paths);
end

paths.rawDataConfigured = strlength(string(paths.dbasePath)) > 0 && ...
    strlength(string(paths.files2runDir)) > 0;

if P.requireRawData
    assert(paths.rawDataConfigured, 'paper_paths:RawDataNotConfigured', ...
        '%s', paths.rawDataInstructions);
    assert(isfile(paths.dbasePath), 'paper_paths:MissingDbase', ...
        'Raw database not found: %s\n%s', paths.dbasePath, paths.rawDataInstructions);
    assert(isfolder(paths.files2runDir), 'paper_paths:MissingFiles2Run', ...
        'Session directory not found: %s\n%s', paths.files2runDir, ...
        paths.rawDataInstructions);
end
end
