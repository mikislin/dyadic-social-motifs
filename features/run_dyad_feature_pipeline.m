function results = run_dyad_feature_pipeline(dbase, fps, nodeMap, opts)
%RUN_DYAD_FEATURE_PIPELINE End-to-end dyadic feature pipeline across sessions.

arguments
    dbase (1,:) struct
    fps (1,1) double {mustBePositive}
    nodeMap struct
    opts.windowSec (1,:) double = [0.2 0.5 1 2]
    opts.computeOpts struct = struct()
    opts.windowOpts struct = struct()
    opts.selectionOpts struct = struct()
    opts.readinessOpts struct = struct()
end

n = numel(dbase);
sessions = repmat(struct('dyad', [], 'windows', []), 1, n);
for i = 1:n
    cOpts = opts.computeOpts;
    if isfield(dbase(i), 'badframes')
        cOpts.badframes = dbase(i).badframes;
    end
    cArgs = namedargs2cell(cOpts);
    sessions(i).dyad = compute_dyad_features(dbase(i).tracks, fps, nodeMap, cArgs{:});

    wOpts = opts.windowOpts;
    wOpts.sessionIdx = i;
    wArgs = namedargs2cell(wOpts);
    sessions(i).windows = summarize_windows(sessions(i).dyad, fps, opts.windowSec, wArgs{:});
end

matrix = build_feature_matrix([sessions.windows]);
inspection = inspect_feature_table(matrix.X, matrix.names);
sArgs = namedargs2cell(opts.selectionOpts);
rArgs = namedargs2cell(opts.readinessOpts);
selection = select_features_for_clustering(matrix.X, matrix.names, sArgs{:});
readiness = assess_clustering_readiness(matrix.X, matrix.table, selection, rArgs{:});

results = struct();
results.sessions = sessions;
results.matrix = matrix;
results.inspection = inspection;
results.selection = selection;
results.readiness = readiness;
results.nodeMap = nodeMap;
results.fps = fps;
end
