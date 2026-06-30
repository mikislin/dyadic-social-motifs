function dyad = compute_dyad_features(tracks, fps, nodeMap, opts)
%COMPUTE_DYAD_FEATURES Compute frame-level dyadic social features.
%
% Inputs
%   tracks   : [T x N x 2 x 2] cleaned tracks for two animals.
%   fps      : frame rate in Hz.
%   nodeMap  : canonical SLEAP node map from default_sleap_node_map().
%              If omitted or empty, the canonical map is loaded from config.
%   opts     : name-value options
%              .contactThresholdMM = 20
%              .closeThresholdMM   = 50
%              .pixelSizeMM        = 1
%              .smoothSpanFrames   = 5
%              .badframes          = [] or [T x 1] / [T x 2] logical
%
% Output
%   dyad is a struct with:
%       .time_s, .frameMask, .badframeMask, .coreNanMask
%       .X                [T x F] numeric matrix; invalid frames are NaN
%       .featureNames     {1 x F}
%       .featureMeta      canonical feature dictionary table
%       .raw              struct with named feature vectors
%       .nodeMetadata     canonical SLEAP bodypart metadata table
%       .params           feature extraction parameters

arguments
    tracks (:,:,:,:) double
    fps (1,1) double {mustBePositive}
    nodeMap struct = struct()
    opts.contactThresholdMM (1,1) double {mustBePositive} = 20
    opts.closeThresholdMM (1,1) double {mustBePositive} = 40
    opts.pixelSizeMM (1,1) double {mustBePositive} = 1
    opts.smoothSpanFrames (1,1) double {mustBeInteger,mustBePositive} = 5
    opts.badframes = []
end

assert(ndims(tracks) == 4, 'tracks must be [T x N x 2 x 2].');
assert(size(tracks,3) == 2, 'tracks third dimension must be x/y.');
assert(size(tracks,4) == 2, 'tracks must contain exactly 2 animals.');

T = size(tracks,1);
time_s = (0:T-1)' ./ fps;
px2mm = opts.pixelSizeMM;

[nodeMap, partNames, nodeMetadata] = local_canonical_node_map(nodeMap);
local_validate_tracks_cover_nodes(tracks, nodeMap);

[badframeMask, badframesByAnimal] = local_normalize_badframes(opts.badframes, T);
coreNanMask = local_core_nan_mask(tracks, nodeMap);

% Positions [T x 2]. Bad frames are made missing before feature computation
% so frame-level artifacts cannot leak through consumers that ignore masks.
p1 = local_extract_points(tracks, nodeMap, 1, badframeMask);
p2 = local_extract_points(tracks, nodeMap, 2, badframeMask);

c1 = mean(cat(3, p1.body, p1.midBody, p1.tailBase), 3, 'omitnan');
c2 = mean(cat(3, p2.body, p2.midBody, p2.tailBase), 3, 'omitnan');

frameMask = ~(badframeMask | coreNanMask);

% Head/body axis: tailBase -> nose.
a1 = normalize_rows(p1.nose - p1.tailBase);
a2 = normalize_rows(p2.nose - p2.tailBase);

% Orthogonal axes for egocentric coordinates.
a1_lat = [-a1(:,2), a1(:,1)];
a2_lat = [-a2(:,2), a2(:,1)];

v12 = c2 - c1; % animal 2 relative to animal 1
v21 = -v12;

head_heading_diff = signed_angle_deg(a1, a2);
cos_head_alignment = row_dot(a1, a2);

% Egocentric partner coordinates.
partner_long_1 = row_dot(v12, a1) * px2mm;
partner_lat_1  = row_dot(v12, a1_lat) * px2mm;
partner_long_2 = row_dot(v21, a2) * px2mm;
partner_lat_2  = row_dot(v21, a2_lat) * px2mm;

nose2nose_dist      = row_norm(p1.nose     - p2.nose)     * px2mm;
head2head_dist      = nose2nose_dist;
centroid_dist       = row_norm(c1          - c2)          * px2mm;
body2body_dist      = row_norm(p1.body     - p2.body)     * px2mm;
tailbase2tailbase   = row_norm(p1.tailBase - p2.tailBase) * px2mm;
nose1_to_body2_dist = row_norm(p1.nose     - p2.body)     * px2mm;
nose2_to_body1_dist = row_norm(p2.nose     - p1.body)     * px2mm;
nose1_to_tail2_dist = row_norm(p1.nose     - p2.tailBase) * px2mm;
nose2_to_tail1_dist = row_norm(p2.nose     - p1.tailBase) * px2mm;

% Facing / line-of-sight scores.
facing_1_to_2 = row_dot(a1, normalize_rows(c2 - p1.nose));
facing_2_to_1 = row_dot(a2, normalize_rows(c1 - p2.nose));
mutual_facing = min(facing_1_to_2, facing_2_to_1);

% Relative kinematics in millimeters and seconds.
vel1 = gradient_by_time(c1 * px2mm, fps);
vel2 = gradient_by_time(c2 * px2mm, fps);
acc1 = gradient_by_time(vel1, fps);
acc2 = gradient_by_time(vel2, fps);

rel_vel = vel2 - vel1;
rhat12 = normalize_rows(v12);
radial_speed_12 = row_dot(rel_vel, rhat12);
tangential_speed_12 = row_dot(rel_vel, [-rhat12(:,2), rhat12(:,1)]);
approach_speed_1 = row_dot(vel1 - vel2, rhat12);
approach_speed_2 = row_dot(vel2 - vel1, normalize_rows(v21));

speed_alignment = cosine_similarity_rows(vel1, vel2);
accel_alignment = cosine_similarity_rows(acc1, acc2);

% Relative bearing of nose and body targets.
nose_bearing_1 = signed_angle_deg(a1, normalize_rows(p2.nose - p1.nose));
nose_bearing_2 = signed_angle_deg(a2, normalize_rows(p1.nose - p2.nose));
body_bearing_1 = signed_angle_deg(a1, normalize_rows(p2.body - p1.nose));
body_bearing_2 = signed_angle_deg(a2, normalize_rows(p1.body - p2.nose));

% Contact / proximity.
in_contact = centroid_dist <= opts.contactThresholdMM;
head_close  = head2head_dist <= opts.contactThresholdMM;
body_close  = body2body_dist <= opts.contactThresholdMM;
close_pair  = centroid_dist <= opts.closeThresholdMM;

% Simple interaction-state indicators to aid inspection.
mutual_approach = (approach_speed_1 > 0) & (approach_speed_2 > 0);
withdrawal      = (approach_speed_1 < 0) & (approach_speed_2 < 0);
asym_investigate = (nose1_to_body2_dist < nose2_to_body1_dist) - ...
                   (nose2_to_body1_dist < nose1_to_body2_dist);

num = @(x) smooth_if_numeric(x, opts.smoothSpanFrames);
circ = @(x) smooth_circular_deg(x, opts.smoothSpanFrames);

raw = struct();
raw.centroid_dist = num(centroid_dist);
raw.body2body_dist = num(body2body_dist);
raw.head2head_dist = num(head2head_dist);
raw.tailbase2tailbase_dist = num(tailbase2tailbase);
raw.nose1_to_body2_dist = num(nose1_to_body2_dist);
raw.nose2_to_body1_dist = num(nose2_to_body1_dist);
raw.nose1_to_tail2_dist = num(nose1_to_tail2_dist);
raw.nose2_to_tail1_dist = num(nose2_to_tail1_dist);
raw.partner_long_1 = num(partner_long_1);
raw.partner_lat_1 = num(partner_lat_1);
raw.partner_long_2 = num(partner_long_2);
raw.partner_lat_2 = num(partner_lat_2);
raw.facing_1_to_2 = num(facing_1_to_2);
raw.facing_2_to_1 = num(facing_2_to_1);
raw.mutual_facing = num(mutual_facing);
raw.heading_diff_deg = circ(head_heading_diff);
raw.cos_head_alignment = num(cos_head_alignment);
raw.radial_speed_12 = num(radial_speed_12);
raw.tangential_speed_12 = num(tangential_speed_12);
raw.approach_speed_1 = num(approach_speed_1);
raw.approach_speed_2 = num(approach_speed_2);
raw.speed_alignment = num(speed_alignment);
raw.accel_alignment = num(accel_alignment);
raw.nose_bearing_1_deg = circ(nose_bearing_1);
raw.nose_bearing_2_deg = circ(nose_bearing_2);
raw.body_bearing_1_deg = circ(body_bearing_1);
raw.body_bearing_2_deg = circ(body_bearing_2);
raw.in_contact = double(in_contact);
raw.head_close = double(head_close);
raw.body_close = double(body_close);
raw.close_pair = double(close_pair);
raw.mutual_approach = double(mutual_approach);
raw.withdrawal = double(withdrawal);
raw.asym_investigate = double(asym_investigate);

[featureNames, featureMeta] = default_dyad_feature_metadata();
X = zeros(T, numel(featureNames));
for k = 1:numel(featureNames)
    X(:,k) = raw.(featureNames{k});
end

X(~frameMask, :) = NaN;
raw = local_mask_raw_features(raw, featureNames, frameMask);

dyad = struct();
dyad.time_s = time_s;
dyad.frameMask = frameMask;
dyad.badframeMask = badframeMask;
dyad.badframesByAnimal = badframesByAnimal;
dyad.coreNanMask = coreNanMask;
dyad.X = X;
dyad.featureNames = featureNames;
dyad.featureMeta = featureMeta;
dyad.raw = raw;
dyad.nodeMap = nodeMap;
dyad.partNames = partNames;
dyad.nodeMetadata = nodeMetadata;
dyad.fps = fps;
dyad.params = struct( ...
    'fps', fps, ...
    'pixelSizeMM', px2mm, ...
    'contactThresholdMM', opts.contactThresholdMM, ...
    'closeThresholdMM', opts.closeThresholdMM, ...
    'smoothSpanFrames', opts.smoothSpanFrames, ...
    'featureMetadataSource', 'features/default_dyad_feature_metadata.m', ...
    'nodeMetadataSource', 'config/sleap_bodypart_metadata.csv');
dyad.maskAudit = struct( ...
    'nFrames', T, ...
    'nValidFrames', nnz(frameMask), ...
    'validFrameFraction', mean(frameMask), ...
    'nBadframeRows', nnz(badframeMask), ...
    'badframeFraction', mean(badframeMask), ...
    'nCoreNanFrames', nnz(coreNanMask), ...
    'coreNanFraction', mean(coreNanMask), ...
    'nCoreNanOnlyFrames', nnz(coreNanMask & ~badframeMask), ...
    'coreNanOnlyFraction', mean(coreNanMask & ~badframeMask), ...
    'featureMatrixNanFraction', mean(isnan(X), 'all'));
end

function [nodeMap, partNames, nodeMetadata] = local_canonical_node_map(nodeMapIn)
[canonicalMap, partNames, nodeMetadata] = default_sleap_node_map();
canonicalFields = string(nodeMetadata.matlab_field(:));

if isempty(fieldnames(nodeMapIn))
    nodeMap = canonicalMap;
    return
end

missing = canonicalFields(~isfield(nodeMapIn, cellstr(canonicalFields)));
assert(isempty(missing), 'compute_dyad_features:NoncanonicalNodeMap', ...
    'nodeMap must contain canonical SLEAP fields from default_sleap_node_map. Missing: %s', ...
    strjoin(missing, ', '));

for i = 1:numel(canonicalFields)
    f = char(canonicalFields(i));
    assert(isequal(nodeMapIn.(f), canonicalMap.(f)), ...
        'compute_dyad_features:NoncanonicalNodeMap', ...
        'nodeMap.%s=%d does not match canonical SLEAP index %d.', ...
        f, nodeMapIn.(f), canonicalMap.(f));
end
nodeMap = canonicalMap;
end

function local_validate_tracks_cover_nodes(tracks, nodeMap)
required = {'nose','neck','body','midBody','tailBase','tailMid'};
for i = 1:numel(required)
    f = required{i};
    assert(isfield(nodeMap, f), 'compute_dyad_features:MissingNode', ...
        'Canonical nodeMap.%s is required.', f);
    assert(nodeMap.(f) >= 1 && nodeMap.(f) <= size(tracks,2), ...
        'compute_dyad_features:NodeIndexOutOfRange', ...
        'nodeMap.%s=%d exceeds track node count %d.', f, nodeMap.(f), size(tracks,2));
end
end

function [badframeMask, badframesByAnimal] = local_normalize_badframes(badframes, T)
if isempty(badframes)
    badframesByAnimal = false(T, 0);
    badframeMask = false(T, 1);
    return
end

bf = logical(badframes);
if isvector(bf)
    bf = bf(:);
end
assert(size(bf,1) == T, 'compute_dyad_features:BadframeSize', ...
    'opts.badframes must have T rows.');

if size(bf,2) >= 2
    badframesByAnimal = bf(:,1:2);
else
    badframesByAnimal = bf(:,1);
end
badframeMask = any(badframesByAnimal, 2);
end

function coreNanMask = local_core_nan_mask(tracks, nodeMap)
p1nose = squeeze(tracks(:, nodeMap.nose, :, 1));
p2nose = squeeze(tracks(:, nodeMap.nose, :, 2));
p1body = squeeze(tracks(:, nodeMap.body, :, 1));
p2body = squeeze(tracks(:, nodeMap.body, :, 2));
p1tail = squeeze(tracks(:, nodeMap.tailBase, :, 1));
p2tail = squeeze(tracks(:, nodeMap.tailBase, :, 2));
p1mid = squeeze(tracks(:, nodeMap.midBody, :, 1));
p2mid = squeeze(tracks(:, nodeMap.midBody, :, 2));
c1 = mean(cat(3, p1body, p1mid, p1tail), 3, 'omitnan');
c2 = mean(cat(3, p2body, p2mid, p2tail), 3, 'omitnan');
corePts = [c1 c2 p1nose p2nose p1body p2body p1tail p2tail];
coreNanMask = any(isnan(corePts), 2);
end

function p = local_extract_points(tracks, nodeMap, animal, badframeMask)
p = struct();
p.nose     = local_node_xy(tracks, nodeMap.nose,     animal, badframeMask);
p.neck     = local_node_xy(tracks, nodeMap.neck,     animal, badframeMask);
p.body     = local_node_xy(tracks, nodeMap.body,     animal, badframeMask);
p.midBody  = local_node_xy(tracks, nodeMap.midBody,  animal, badframeMask);
p.tailBase = local_node_xy(tracks, nodeMap.tailBase, animal, badframeMask);
p.tailMid  = local_node_xy(tracks, nodeMap.tailMid,  animal, badframeMask);
end

function xy = local_node_xy(tracks, nodeIdx, animal, badframeMask)
xy = squeeze(tracks(:, nodeIdx, :, animal));
xy(badframeMask, :) = NaN;
end

function raw = local_mask_raw_features(raw, featureNames, frameMask)
for k = 1:numel(featureNames)
    fn = featureNames{k};
    x = raw.(fn);
    x(~frameMask) = NaN;
    raw.(fn) = x;
end
end

function x = smooth_if_numeric(x, span)
if islogical(x) || all(ismember(unique(x(~isnan(x))), [0 1]))
    return
end
if span > 1
    x = smoothdata(x, 1, 'movmean', span, 'omitnan');
end
end

function x = smooth_circular_deg(x, span)
if span <= 1
    return
end
r = deg2rad(x);
c = smoothdata(cos(r), 1, 'movmean', span, 'omitnan');
s = smoothdata(sin(r), 1, 'movmean', span, 'omitnan');
x = rad2deg(atan2(s, c));
x(~isfinite(c) | ~isfinite(s)) = NaN;
end

function out = normalize_rows(x)
d = sqrt(sum(x.^2, 2));
d(d == 0) = NaN;
out = x ./ d;
end

function y = row_norm(x)
y = sqrt(sum(x.^2, 2));
end

function y = row_dot(a, b)
y = sum(a .* b, 2);
end

function ang = signed_angle_deg(a, b)
a = normalize_rows(a);
b = normalize_rows(b);
cross2d = a(:,1).*b(:,2) - a(:,2).*b(:,1);
dotp = row_dot(a, b);
ang = atan2d(cross2d, dotp);
end

function c = cosine_similarity_rows(a, b)
da = row_norm(a);
db = row_norm(b);
den = da .* db;
den(den == 0) = NaN;
c = row_dot(a, b) ./ den;
end

function g = gradient_by_time(x, fps)
if isvector(x)
    x = x(:);
end
g = [zeros(1, size(x,2)); diff(x, 1, 1)] .* fps;
end
