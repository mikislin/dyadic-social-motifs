function dyad = compute_dyad_features(tracks, fps, nodeMap, opts)
%COMPUTE_DYAD_FEATURES Compute frame-level dyadic social features.
%   dyad = compute_dyad_features(tracks, fps, nodeMap)
%
% Inputs
%   tracks   : [T x N x 2 x 2] cleaned tracks for two animals
%   fps      : frame rate
%   nodeMap  : struct with node indices, e.g.
%              nose, neck, body, tailBase, tailMid, midBody
%   opts     : optional struct fields
%              .contactThresholdMM = 30
%              .closeThresholdMM   = 60
%              .pixelSizeMM        = 1
%              .smoothSpanFrames   = 5
%              .badframes          = [] or [T x 1] / [T x 2] logical
%
% Output
%   dyad is a struct with fields:
%       .time_s
%       .frameMask
%       .X                [T x F] numeric matrix
%       .featureNames     {1 x F}
%       .featureMeta      table describing each feature
%       .raw              struct with named vectors
%
% The package is intentionally dyad-centered. Single-animal kinematics are
% assumed to be available elsewhere and are not recomputed here.

arguments
    tracks (:,:,:,:) double
    fps (1,1) double {mustBePositive}
    nodeMap struct
    opts.contactThresholdMM (1,1) double {mustBePositive} = 30
    opts.closeThresholdMM (1,1) double {mustBePositive} = 60
    opts.pixelSizeMM (1,1) double {mustBePositive} = 1
    opts.smoothSpanFrames (1,1) double {mustBeInteger,mustBePositive} = 5
    opts.badframes = []
end

assert(ndims(tracks) == 4, 'tracks must be [T x N x 2 x 2].');
assert(size(tracks,3) == 2, 'tracks third dimension must be x/y.');
assert(size(tracks,4) == 2, 'tracks must contain exactly 2 animals.');

T = size(tracks,1);
time_s = (0:T-1)' ./ fps;

reqFields = {'nose','body','tailBase'};
for i = 1:numel(reqFields)
    assert(isfield(nodeMap, reqFields{i}), 'nodeMap.%s is required.', reqFields{i});
end

% Optional fallback fields.
if ~isfield(nodeMap, 'neck');    nodeMap.neck    = nodeMap.body;     end
if ~isfield(nodeMap, 'midBody'); nodeMap.midBody = nodeMap.body;     end
if ~isfield(nodeMap, 'tailMid'); nodeMap.tailMid = nodeMap.tailBase; end

px2mm = opts.pixelSizeMM;

% Positions [T x 2].
p1.nose     = squeeze(tracks(:, nodeMap.nose,     :, 1));
p1.neck     = squeeze(tracks(:, nodeMap.neck,     :, 1));
p1.body     = squeeze(tracks(:, nodeMap.body,     :, 1));
p1.midBody  = squeeze(tracks(:, nodeMap.midBody,  :, 1));
p1.tailBase = squeeze(tracks(:, nodeMap.tailBase, :, 1));
p1.tailMid  = squeeze(tracks(:, nodeMap.tailMid,  :, 1));

p2.nose     = squeeze(tracks(:, nodeMap.nose,     :, 2));
p2.neck     = squeeze(tracks(:, nodeMap.neck,     :, 2));
p2.body     = squeeze(tracks(:, nodeMap.body,     :, 2));
p2.midBody  = squeeze(tracks(:, nodeMap.midBody,  :, 2));
p2.tailBase = squeeze(tracks(:, nodeMap.tailBase, :, 2));
p2.tailMid  = squeeze(tracks(:, nodeMap.tailMid,  :, 2));

c1 = nanmean(cat(3, p1.body, p1.midBody, p1.tailBase), 3);
c2 = nanmean(cat(3, p2.body, p2.midBody, p2.tailBase), 3);

% Head/body axis: tailBase -> nose.
a1 = normalize_rows(p1.nose - p1.tailBase);
a2 = normalize_rows(p2.nose - p2.tailBase);

% Orthogonal axes for egocentric coordinates.
a1_lat = [ -a1(:,2), a1(:,1) ];
a2_lat = [ -a2(:,2), a2(:,1) ];

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

% Relative kinematics.
vel1 = gradient_by_time(c1 * px2mm, fps);
vel2 = gradient_by_time(c2 * px2mm, fps);
acc1 = gradient_by_time(vel1, fps);
acc2 = gradient_by_time(vel2, fps);

speed1 = row_norm(vel1);
speed2 = row_norm(vel2);
rel_vel = vel2 - vel1;

rhat12 = normalize_rows(v12);
radial_speed_12 = row_dot(rel_vel, rhat12);
tangential_speed_12 = row_dot(rel_vel, [ -rhat12(:,2), rhat12(:,1) ]);
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

% Smooth numeric features lightly to reduce frame jitter.
num = @(x) smooth_if_numeric(x, opts.smoothSpanFrames);

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
raw.heading_diff_deg = num(head_heading_diff);
raw.cos_head_alignment = num(cos_head_alignment);
raw.radial_speed_12 = num(radial_speed_12);
raw.tangential_speed_12 = num(tangential_speed_12);
raw.approach_speed_1 = num(approach_speed_1);
raw.approach_speed_2 = num(approach_speed_2);
raw.speed_alignment = num(speed_alignment);
raw.accel_alignment = num(accel_alignment);
raw.nose_bearing_1_deg = num(nose_bearing_1);
raw.nose_bearing_2_deg = num(nose_bearing_2);
raw.body_bearing_1_deg = num(body_bearing_1);
raw.body_bearing_2_deg = num(body_bearing_2);
raw.in_contact = double(in_contact);
raw.head_close = double(head_close);
raw.body_close = double(body_close);
raw.close_pair = double(close_pair);
raw.mutual_approach = double(mutual_approach);
raw.withdrawal = double(withdrawal);
raw.asym_investigate = double(asym_investigate);

featureNames = fieldnames(raw)';
X = zeros(T, numel(featureNames));
for k = 1:numel(featureNames)
    X(:,k) = raw.(featureNames{k});
end

% Frame mask: valid if enough coordinates exist and not marked bad.
frameMask = true(T,1);
corePts = [c1 c2 p1.nose p2.nose p1.body p2.body p1.tailBase p2.tailBase];
frameMask(any(isnan(corePts),2)) = false;

if ~isempty(opts.badframes)
    bf = opts.badframes;
    if ismatrix(bf) && size(bf,1) == T
        frameMask(any(logical(bf),2)) = false;
    else
        error('opts.badframes must have T rows.');
    end
end

families = {
    'distance','distance','distance','distance', ...
    'distance','distance','distance','distance', ...
    'egocentric','egocentric','egocentric','egocentric', ...
    'orientation','orientation','orientation','orientation', ...
    'kinematics','kinematics','kinematics','kinematics','kinematics', ...
    'coupling','coupling', ...
    'orientation','orientation','orientation','orientation', ...
    'contact','contact','contact','contact', ...
    'interaction_logic','interaction_logic','interaction_logic'}';

isCircular = contains(featureNames, 'bearing_')' | contains(featureNames, 'heading_diff')';
isBoolean  = ismember(featureNames, {'in_contact','head_close','body_close','close_pair','mutual_approach','withdrawal'});
isDirected = contains(featureNames, '_1') | contains(featureNames, '_2') | contains(featureNames, '1_to_') | contains(featureNames, '2_to_');

transformHint = repmat("none", numel(featureNames), 1);
for k = 1:numel(featureNames)
    if contains(featureNames{k}, 'dist')
        transformHint(k) = "log1p";
    elseif isBoolean(k)
        transformHint(k) = "binary";
    elseif isCircular(k)
        transformHint(k) = "circular";
    else
        transformHint(k) = "zscore";
    end
end

featureMeta = table(featureNames(:), families, isDirected(:), isCircular(:), isBoolean(:), transformHint, ...
    'VariableNames', {'Name','Family','IsDirected','IsCircular','IsBoolean','TransformHint'});

dyad = struct();
dyad.time_s = time_s;
dyad.frameMask = frameMask;
dyad.X = X;
dyad.featureNames = featureNames;
dyad.featureMeta = featureMeta;
dyad.raw = raw;
end

function x = smooth_if_numeric(x, span)
if islogical(x) || all(ismember(unique(x(~isnan(x))), [0 1]))
    return;
end
if span > 1
    x = smoothdata(x, 1, 'movmean', span, 'omitnan');
end
end

function out = normalize_rows(x)
d = sqrt(sum(x.^2, 2));
d(d == 0) = NaN;
out = x ./ d;
end

function y = row_norm(x)
y = sqrt(sum(x.^2,2));
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
c = row_dot(a,b) ./ den;
end

function g = gradient_by_time(x, fps)
if isvector(x)
    x = x(:);
end
g = [zeros(1,size(x,2)); diff(x,1,1)] .* fps;
end
