function tests = test_preprocess_session
%TEST_PREPROCESS_SESSION MATLAB unit tests for dyadic preprocessing.
tests = functiontests(localfunctions);
end

function setupOnce(~)
rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(rootDir, 'preprocessing'));
addpath(fullfile(rootDir, 'io'));
addpath(fullfile(rootDir, 'cluster'));
end

function testDefaultParamsValidate(testCase)
P = default_preprocessing_params();
P = validate_preprocessing_params(P);
verifyEqual(testCase, P.data.fps, 80);
verifyEqual(testCase, P.confidence.threshold, 0.10, 'AbsTol', 1e-12);
verifyEqual(testCase, P.jump.max_disp_mm_per_frame_by_node(10), 12);
end

function testSessionValidationAddsDefaults(testCase)
P = default_preprocessing_params();
session = makeSyntheticSession(100, 13, 2, false, false, false);
session = rmfield(session, 'time');
session = rmfield(session, 'excludedFrames');
sessionOut = validate_session_inputs(session, P);
verifySize(testCase, sessionOut.time, [100 1]);
verifySize(testCase, sessionOut.excludedFrames, [100 1]);
verifyTrue(testCase, all(sessionOut.excludedFrames == 0));
end

function testConfidenceMaskApplies(testCase)
tracks = rand(20, 13, 2);
scores = ones(20,13);
scores(5,2) = 0.1;
[tracksMasked, mask] = apply_confidence_mask(tracks, scores, 0.2);
verifyTrue(testCase, mask(5,2));
verifyTrue(testCase, all(isnan(squeeze(tracksMasked(5,2,:)))));
end

function testConfidenceMaskNodeSpecific(testCase)
tracks = rand(20, 13, 2);
scores = 0.3 * ones(20,13);
thr = 0.1 * ones(1,13);
thr(4) = 0.5;
[tracksMasked, mask] = apply_confidence_mask(tracks, scores, thr);
verifyTrue(testCase, any(mask(:,4)));
verifyFalse(testCase, any(mask(:,3)));
verifyTrue(testCase, all(isnan(tracksMasked(mask(:,4),4,1))));
end

function testInterpolateShort2DGap(testCase)
time = (0:19)';
xy = [(1:20)' (101:120)'];
xy(5:6,:) = NaN;
xy(10:15,:) = NaN;
gp = struct('method','linear','max_gap_frames',2,'fill_ends',false,'min_points_for_interp',2);
[xyi, mask] = interpolate_short_gaps_2d(xy, time, gp);
verifyTrue(testCase, all(all(isfinite(xyi(5:6,:)))));
verifyTrue(testCase, all(any(isnan(xyi(10:15,:)),2)));
verifyEqual(testCase, sum(mask), 2);
end

function testJumpFilter2DReducesTeleport(testCase)
time = (0:99)'/80;
xy = [sin(time) cos(time)];
xy(50,:) = xy(49,:) + [100 100];
params = struct('diff_quantile',0.95,'max_iter',10,'min_valid_points',5, ...
    'max_disp_mm_per_frame',25,'max_disp_px_per_frame',50,'method','linear','dilate_frames',0);
[xyOut, jumpMask, meta] = iterative_jump_filter_2d(xy, time, true(size(time)), params);
verifyTrue(testCase, any(jumpMask));
verifyLessThan(testCase, norm(xyOut(50,:) - xyOut(49,:)), 20);
verifyGreaterThanOrEqual(testCase, meta.n_iter, 1);
end

function testSmoothTrackBreaksAtRepairs(testCase)
xy = [(1:100)' (1:100)'];
xy(50,:) = [500 500];
sp = struct('method','lowpass','low_pass_frq',4,'n_pad',2,'min_segment_length',5,'break_at_repairs',true,'break_dilate_frames',1);
[xy2, meta] = smooth_track_2d_segments(xy, sp, 80, [false(49,1); true; false(50,1)]);
verifyGreaterThan(testCase, meta.n_breakpoints, 0);
verifyEqual(testCase, xy2(50,:), xy(50,:));
end

function testGeometryBlameAssignsOneNode(testCase)
tracks = nan(100,13,2);
for n = 1:13
    tracks(:,n,1) = n;
    tracks(:,n,2) = 0;
end
tracks(40,1,1) = 500;
scores = ones(100,13);
scores(40,1) = 0.05;
scores(40,10) = 0.9;
P = default_preprocessing_params();
[mask, info] = detect_geometry_outliers(tracks, [1 10], [1 99], scores, P);
verifyTrue(testCase, mask(40,1));
verifyFalse(testCase, mask(40,10));
verifyEqual(testCase, info.blameNode(40,1), 1);
end

function testSingleAnimalSessionRuns(testCase)
P = default_preprocessing_params();
P.smooth.enabled = false;
P.output.return_raw = false;

session = makeSyntheticSession(120, 13, 1, true, true, true);
out = preprocess_session(session, P);

verifySize(testCase, out.clean.tracks, [120 13 2 1]);
verifySize(testCase, out.qc.badframes, [120 1]);
verifyEqual(testCase, numel(out.qc.animals), 1);
end

function testPreprocessSessionReturnsExpectedFields(testCase)
P = default_preprocessing_params();
P.smooth.enabled = false;
P.output.return_raw = true;
session = makeSyntheticSession(120, 13, 2, true, true, true);
out = preprocess_session(session, P);
verifyTrue(testCase, isfield(out, 'clean'));
verifyTrue(testCase, isfield(out.clean, 'tracks'));
verifySize(testCase, out.clean.tracks, size(session.SLEAPtracks));
verifySize(testCase, out.qc.badframes, [120 2]);
verifyTrue(testCase, isfield(out.debug, 'frameHasAnyNaN'));
end

function testExcludedFramesPropagateToBadframes(testCase)
P = default_preprocessing_params();
P.smooth.enabled = false;
session = makeSyntheticSession(80, 13, 2, false, false, false);
session.excludedFrames(10:15) = true;
out = preprocess_session(session, P);
verifyTrue(testCase, all(out.qc.badframes(10:15,1)));
verifyTrue(testCase, all(out.qc.badframes(10:15,2)));
end

function testArenaMaskCatchesExtremeExcursion(testCase)
P = default_preprocessing_params();
P.smooth.enabled = false;
P.qc.arena_percentile = [5 95];
P.qc.arena_margin_px = 5;
session = makeSyntheticSession(150, 13, 2, false, false, false);
session.SLEAPtracks(100,10,1,1) = 1000;
session.SLEAPtracks(100,10,2,1) = 1000;
out = preprocess_session(session, P);
verifyTrue(testCase, full(out.qc.animals(1).arenaMask(100,10)) || out.qc.badframes(100,1));
end

function testLowConfidenceMaskCreatesBadframeSignal(testCase)
P = default_preprocessing_params();
P.smooth.enabled = false;
session = makeSyntheticSession(100, 13, 2, false, false, true);
session.SLEAPscores(30, :, :) = 0.05;
out = preprocess_session(session, P);
verifyGreaterThan(testCase, out.qc.frames.fracLowConf(30,1), 0.5);
end

function testDistortionAuditPlotRuns(testCase)
P = default_preprocessing_params();
P.smooth.enabled = false;
P.output.return_raw = true;
session = makeSyntheticSession(80, 13, 2, true, true, true);
out = preprocess_session(session, P);
fig = plot_preprocessing_qc(out);
verifyTrue(testCase, isgraphics(fig));
close(fig);
end

function session = makeSyntheticSession(T, nNodes, nAnimals, addShortGap, addTeleport, addScores)
t = (0:T-1)' ./ 80;
tracks = nan(T, nNodes, 2, nAnimals);
for m = 1:nAnimals
    ctr = [100 + 30*m, 150 + 25*m];
    bodyLen = 25 + 2*m;
    for n = 1:nNodes
        phase = 2*pi*(0.01*n);
        offset = [2*n, -n];
        x = ctr(1) + offset(1) + 15*sin(2*pi*0.2*t + phase);
        y = ctr(2) + offset(2) + 12*cos(2*pi*0.2*t + phase);
        if n == 1
            x = x + bodyLen;
        elseif n == 10
            x = x - bodyLen;
        end
        tracks(:,n,1,m) = x;
        tracks(:,n,2,m) = y;
    end
end

if addShortGap
    tracks(20:23, 2, :, 1) = NaN;
end
if addTeleport
    tracks(60, 3, 1, 1) = tracks(59, 3, 1, 1) + 120;
    tracks(60, 3, 2, 1) = tracks(59, 3, 2, 1) - 90;
    if nAnimals >= 2
        tracks(60, 3, 1, 2) = tracks(59, 3, 1, 2) + 150;
        tracks(60, 3, 2, 2) = tracks(59, 3, 2, 2) - 120;
    end
end

session = struct();
session.SLEAPtracks = tracks;
session.time = t;
session.excludedFrames = false(T,1);
session.session_id = 'synthetic_session';
if addScores
    session.SLEAPscores = ones(T, nNodes, nAnimals);
else
    session.SLEAPscores = [];
end
end
