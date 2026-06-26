%RUN_01_MANIFEST Build the paper session manifest.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
cd(repoRoot);
addpath(genpath(repoRoot));

paths = paper_paths('requireRawData', true);

manifest = build_session_manifest(paths.files2runDir, paths.dbasePath, ...
    'writePath', paths.manifestPath, ...
    'verbose', true);

fprintf('Wrote manifest: %s\n', paths.manifestPath);
fprintf('Rows: %d\n', height(manifest));
