repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));

cfg = load_preprocessing_pipeline_config( ...
    fullfile(repoRoot, 'config', 'preprocessing_pipeline_config.csv'));
bodypointMeta = cfg.bodypoints.bodypoint_metadata;
part_names = cellstr(bodypointMeta.bodypart_name)';

nodeMap = struct();
for iNode = 1:height(bodypointMeta)
    nodeMap.(char(bodypointMeta.matlab_field(iNode))) = bodypointMeta.node_index(iNode);
end

% Aliases commonly expected by different code versions.
nodeMap.left_ear  = nodeMap.leftEar;
nodeMap.right_ear = nodeMap.rightEar;
nodeMap.lf_paw    = nodeMap.lfPaw;
nodeMap.rf_paw    = nodeMap.rfPaw;
nodeMap.lh_paw    = nodeMap.lhPaw;
nodeMap.rh_paw    = nodeMap.rhPaw;
nodeMap.body_pos  = nodeMap.body;
nodeMap.tail_base = nodeMap.tailBase;
nodeMap.mid_tail  = nodeMap.tailMid;
nodeMap.tail_tip  = nodeMap.tailTip;
nodeMap.mid_body  = nodeMap.midBody;

% Helpful semantic aliases
nodeMap.head        = nodeMap.nose;
nodeMap.body_center = nodeMap.midBody;
nodeMap.center      = nodeMap.midBody;

% Feature extraction options
opts = struct();

% General
opts.sessionIdx = 1;
opts.requireBothAnimals = true;

% Window summarization defaults
opts.requireFullWindow = true;
opts.maxMissingFrac = 0.10;

% Contact thresholds in pixels
% Adjust after checking your camera calibration / preprocessing scale
opts.contactThreshPx = 30;
opts.noseToBodyContactThreshPx = 35;
opts.headToHeadThreshPx = 25;

% Approach / withdrawal thresholds
opts.approachSpeedThresh = 5;     % pixels/s
opts.retreatSpeedThresh  = -5;    % pixels/s

% Facing thresholds in degrees
opts.facingAngleThreshDeg = 60;
opts.parallelAngleThreshDeg = 30;
opts.antiparallelAngleThreshDeg = 150;

% If you have mm-per-pixel available, store it here too
opts.pixel_size_mm = 1/1.97;   % from your earlier preprocessing code
opts.fps = 80;

fps = opts.fps;
i = 31;

badframes1d = any(dbase(i).badframes ~= 0, 2);

dyad = compute_dyad_features(dbase(i).tracks, fps, nodeMap, ...
    'badframes', badframes1d, ...
    'pixelSizeMM', opts.pixel_size_mm);

W = summarize_windows(dyad, fps, [0.2 0.5 1.0 2.0], ...
    'requireFullWindow', true, ...
    'maxMissingFrac', 0.10, ...
    'sessionIdx', i);

