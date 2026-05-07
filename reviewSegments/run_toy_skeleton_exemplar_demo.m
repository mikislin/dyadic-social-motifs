%% Toy demo for skeleton-only motif exemplars
% This demo does not need raw video. It creates two simple mouse-like
% skeletons, a toy Cluster/Segments object, and renders exemplar clips.

clear; close all;
fps = 30;
T = 900;
nParts = 13;
nAnimals = 2;
tracks = nan(T, nParts, 2, nAnimals);
baseX = linspace(100, 500, T)';
baseY = 250 + 40*sin((1:T)'/60);
body = [0 0; 20 0; 10 -10; 10 10; 20 -20; 20 20; -10 -20; -10 20; -20 0; -40 0; -55 0; -70 0; -10 0];
for t = 1:T
    for m = 1:nAnimals
        offset = [baseX(t) + 70*(m-1), baseY(t) + 25*cos(t/45 + m)];
        theta = 0.4*sin(t/80 + m);
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        pts = body*R' + offset;
        tracks(t,:,:,m) = pts;
    end
end

anchorFrames = (30:10:870)';
nA = numel(anchorFrames);
labels = repelem((1:3)', ceil(nA/3)); labels = labels(1:nA);
Xpca = [sin(anchorFrames/100), cos(anchorFrames/80), labels + 0.1*randn(nA,1)];
Cluster = struct();
Cluster.NumClusters = 3;
Cluster.labels = labels;
Cluster.maxPosterior = 0.95 + 0.05*rand(nA,1);
Cluster.Xpca = Xpca;
Cluster.Data.anchorTable = table(ones(nA,1), anchorFrames, anchorFrames/fps, ...
    'VariableNames', {'session_index','anchor_frame','anchor_time_s'});

segStart = [1; 20; 40; 60];
segEnd = [19; 39; 59; nA];
segMotif = [1;2;3;2];
Segments = struct();
Segments.table = table((1:4)', ones(4,1), segMotif, anchorFrames(segStart), anchorFrames(segEnd), ...
    anchorFrames(segStart)/fps, anchorFrames(segEnd)/fps, ...
    'VariableNames', {'segment_id','session_index','motif_id','start_frame','end_frame','start_time_s','end_time_s'});

Exemplars = find_motif_video_exemplars(Cluster, Segments, ...
    'TopN', 2, ...
    'WindowSec', 1, ...
    'FrameRate', fps, ...
    'PoseData', {tracks}, ...
    'MakeSkeletonClips', true, ...
    'OutputDir', fullfile(tempdir, 'toy_motif_skeleton_exemplars'), ...
    'Verbose', true);

disp(Exemplars.table(:, {'motif_id','rank','session_index','center_frame','clip_file'}));
