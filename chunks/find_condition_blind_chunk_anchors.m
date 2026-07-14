function Anchor = find_condition_blind_chunk_anchors(Seq, scaleSec, opts)
%FIND_CONDITION_BLIND_CHUNK_ANCHORS Find mask-valid multiscale anchors.
%
% Anchor eligibility uses only time geometry, Seq.validMask, transformed
% feature finite fractions, and caller-provided thresholds. It does not read
% condition, cohort, arena, drug, genotype, or outcome labels.

arguments
    Seq struct
    scaleSec (:,1) double {mustBePositive}
    opts.strideSec (1,1) double {mustBePositive} = 0.25
    opts.anchorMode (1,1) string {mustBeMember(opts.anchorMode, ["center","past"])} = "center"
    opts.minValidFrac (1,1) double {mustBeGreaterThanOrEqual(opts.minValidFrac,0),mustBeLessThanOrEqual(opts.minValidFrac,1)} = 0.85
    opts.minFeatureFiniteFrac (1,1) double {mustBeGreaterThanOrEqual(opts.minFeatureFiniteFrac,0),mustBeLessThanOrEqual(opts.minFeatureFiniteFrac,1)} = 0.95
    opts.requireAnchorFrameValid (1,1) logical = true
end

assert(isfield(Seq, 'time') && isfield(Seq, 'validMask') && isfield(Seq, 'X') && isfield(Seq, 'fps'), ...
    'find_condition_blind_chunk_anchors:BadSeq', ...
    'Seq must contain time, validMask, X, and fps.');

time = Seq.time(:);
validMask = logical(Seq.validMask(:));
X = Seq.X;
T = numel(time);
assert(size(X, 1) == T && numel(validMask) == T, ...
    'find_condition_blind_chunk_anchors:LengthMismatch', ...
    'Seq time, X, and validMask lengths must match.');

fps = double(Seq.fps);
strideFrames = max(1, round(opts.strideSec * fps));
scaleSec = scaleSec(:);
scaleFrames = max(1, round(scaleSec .* fps));
maxFrames = max(scaleFrames);
[leftMax, rightMax] = local_chunk_geometry(maxFrames, opts.anchorMode);

if T < maxFrames
    Anchor = local_empty_anchor_table();
    return
end

anchors = ((leftMax + 1):strideFrames:(T - rightMax))';
if isempty(anchors)
    Anchor = local_empty_anchor_table();
    return
end

validNum = double(validMask);
validCs = [0; cumsum(validNum)];
rowFiniteFrac = mean(isfinite(X), 2);
rowFiniteFrac(~isfinite(rowFiniteFrac)) = 0;
finiteCs = [0; cumsum(rowFiniteFrac)];

nAnchors = numel(anchors);
nScales = numel(scaleSec);
validByScale = nan(nAnchors, nScales);
finiteByScale = nan(nAnchors, nScales);
inBoundsByScale = false(nAnchors, nScales);
ok = true(nAnchors, 1);
if opts.requireAnchorFrameValid
    ok = ok & validMask(anchors);
end

for s = 1:nScales
    L = scaleFrames(s);
    [left, right] = local_chunk_geometry(L, opts.anchorMode);
    st = anchors - left;
    en = anchors + right;
    inBounds = st >= 1 & en <= T;
    inBoundsByScale(:, s) = inBounds;

    vf = nan(nAnchors, 1);
    ff = nan(nAnchors, 1);
    vf(inBounds) = (validCs(en(inBounds) + 1) - validCs(st(inBounds))) ./ L;
    ff(inBounds) = (finiteCs(en(inBounds) + 1) - finiteCs(st(inBounds))) ./ L;
    validByScale(:, s) = vf;
    finiteByScale(:, s) = ff;
    ok = ok & inBounds & vf >= opts.minValidFrac & ff >= opts.minFeatureFiniteFrac;
end

Anchor = table();
Anchor.anchor_frame = anchors(ok);
Anchor.anchor_time_s = time(anchors(ok));
Anchor.n_scales_tested = repmat(nScales, nnz(ok), 1);
Anchor.min_scale_valid_frac = min(validByScale(ok, :), [], 2, 'omitnan');
Anchor.min_scale_feature_finite_frac = min(finiteByScale(ok, :), [], 2, 'omitnan');
Anchor.all_scales_in_bounds = all(inBoundsByScale(ok, :), 2);
Anchor.anchor_frame_valid = validMask(anchors(ok));
Anchor.anchor_candidate_rule = repmat( ...
    "frame_mask_and_feature_availability_only_no_condition_labels", nnz(ok), 1);
end

function [leftFrames, rightFrames] = local_chunk_geometry(L, anchorMode)
switch string(anchorMode)
    case "center"
        leftFrames = floor((L - 1) / 2);
        rightFrames = ceil((L - 1) / 2);
    case "past"
        leftFrames = L - 1;
        rightFrames = 0;
end
end

function T = local_empty_anchor_table()
T = table('Size', [0 8], ...
    'VariableTypes', {'double','double','double','double','double','logical','logical','string'}, ...
    'VariableNames', {'anchor_frame','anchor_time_s','n_scales_tested', ...
    'min_scale_valid_frac','min_scale_feature_finite_frac', ...
    'all_scales_in_bounds','anchor_frame_valid','anchor_candidate_rule'});
end
