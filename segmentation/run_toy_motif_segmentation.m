%RUN_TOY_MOTIF_SEGMENTATION Minimal toy example for build_motif_segments.
clear; close all;

fps = 10;
N = 300;
time = (0:N-1)' ./ fps;
trueLabels = [ones(60,1); 2*ones(80,1); ones(20,1); 3*ones(90,1); 2*ones(50,1)];
rawLabels = trueLabels;
rawLabels([55 56 120 121 122 210]) = [2 2 1 1 3 1]; % flicker/noise
K = 3;
post = 0.05 * ones(N,K);
for i = 1:N
    post(i, rawLabels(i)) = 0.90;
end
post = post ./ sum(post,2);

anchorTable = table(ones(N,1), (1:N)', time, ...
    'VariableNames', {'session_index','anchor_frame','anchor_time_s'});
Cluster = struct();
Cluster.Data.anchorTable = anchorTable;
Cluster.labels = rawLabels;
Cluster.posteriors = post;
Cluster.maxPosterior = max(post,[],2);
Cluster.NumClusters = K;

Segments = build_motif_segments(Cluster, ...
    'SmoothWindowAnchors', 5, ...
    'MinRunAnchors', 4, ...
    'MergeGapAnchors', 3, ...
    'MinSegmentDurationSec', 0.5, ...
    'FrameRate', fps, ...
    'Verbose', true);

SegReport = validate_motif_segments(Segments, Cluster);
Fig = plot_motif_segment_diagnostics(Segments, 'ExampleSession', 1);
