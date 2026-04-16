function run_preprocessing_array(inDir, outDir, taskId)
%RUN_PREPROCESSING_ARRAY Preprocess a single exported session based on task id.

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'preprocessing'));
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'io'));

if nargin < 3 || isempty(taskId)
    taskId = str2double(getenv('SLURM_ARRAY_TASK_ID'));
end
assert(~isnan(taskId) && taskId >= 1, 'Valid task id required');

P = default_preprocessing_params();
P.output.return_raw = false;
P.output.store_full_masks = false;
P.output.make_plots = false;
P.debug.enabled = false;

files = dir(fullfile(inDir, 'session_*.mat'));
assert(taskId <= numel(files), 'Task id exceeds number of session files');

if ~isfolder(outDir)
    mkdir(outDir);
end

inFile = fullfile(files(taskId).folder, files(taskId).name);
S = load(inFile, 'sessionRaw');
out = preprocess_session(S.sessionRaw, P);
[~, name] = fileparts(files(taskId).name);
save(fullfile(outDir, [name '_preproc.mat']), 'out', '-v7.3');
exit
end
