function run_preprocessing_batch(inDir, outDir, P)
%RUN_PREPROCESSING_BATCH Preprocess exported session MAT files.

if nargin < 3 || isempty(P)
    P = default_preprocessing_params();
end
P = validate_preprocessing_params(P);
assert(isfolder(inDir), 'Input directory not found');
if ~isfolder(outDir)
    mkdir(outDir);
end

files = dir(fullfile(inDir, 'session_*.mat'));
assert(~isempty(files), 'No session_*.mat files found in input directory');

for i = 1:numel(files)
    inFile = fullfile(files(i).folder, files(i).name);
    S = load(inFile, 'sessionRaw');
    out = preprocess_session(S.sessionRaw, P);
    [~, name] = fileparts(files(i).name);
    save(fullfile(outDir, [name '_preproc.mat']), 'out', '-v7.3');
end
end
