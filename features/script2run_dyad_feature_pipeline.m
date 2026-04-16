% Part names from SLEAP
part_names = { ...
    'nose','neck','left_ear','right_ear','LF_paw','RF_paw', ...
    'LH_paw','RH_paw','body_pos','tail_base','mid_tail', ...
    'tail_tip','mid_body_pos'};

nodeMap = struct();

% Primary nodes
nodeMap.nose      = lookupNode(part_names, 'nose');
nodeMap.neck      = lookupNode(part_names, 'neck');
nodeMap.left_ear  = lookupNode(part_names, 'left_ear');
nodeMap.right_ear = lookupNode(part_names, 'right_ear');
nodeMap.lf_paw    = lookupNode(part_names, 'LF_paw');
nodeMap.rf_paw    = lookupNode(part_names, 'RF_paw');
nodeMap.lh_paw    = lookupNode(part_names, 'LH_paw');
nodeMap.rh_paw    = lookupNode(part_names, 'RH_paw');
nodeMap.body_pos  = lookupNode(part_names, 'body_pos');
nodeMap.tail_base = lookupNode(part_names, 'tail_base');
nodeMap.mid_tail  = lookupNode(part_names, 'mid_tail');
nodeMap.tail_tip  = lookupNode(part_names, 'tail_tip');
nodeMap.mid_body  = lookupNode(part_names, 'mid_body_pos');

% Aliases commonly expected by different code versions
nodeMap.body        = nodeMap.body_pos;
nodeMap.tailBase    = nodeMap.tail_base;
nodeMap.tailTip     = nodeMap.tail_tip;
nodeMap.midBody     = nodeMap.mid_body;
nodeMap.leftEar     = nodeMap.left_ear;
nodeMap.rightEar    = nodeMap.right_ear;
nodeMap.lfPaw       = nodeMap.lf_paw;
nodeMap.rfPaw       = nodeMap.rf_paw;
nodeMap.lhPaw       = nodeMap.lh_paw;
nodeMap.rhPaw       = nodeMap.rh_paw;

% Helpful semantic aliases
nodeMap.head        = nodeMap.nose;
nodeMap.body_center = nodeMap.mid_body;
nodeMap.center      = nodeMap.mid_body;

function idx = lookupNode(part_names, name)
    idx = find(strcmp(part_names, name));
    assert(~isempty(idx), 'Missing node name: %s', name);
    assert(isscalar(idx), 'Duplicate node name: %s', name);
end

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

