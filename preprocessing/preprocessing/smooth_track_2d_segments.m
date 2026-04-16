function [xySmooth, smoothMeta] = smooth_track_2d_segments(xyIn, smoothParams, fps, repairMask)
%SMOOTH_TRACK_2D_SEGMENTS Smooth contiguous finite segments, optionally split at repairs.

if nargin < 4 || isempty(repairMask)
    repairMask = false(size(xyIn,1),1);
end

xySmooth = xyIn;
smoothMeta = struct('n_segments_smoothed', 0, 'segment_lengths', []);
if strcmpi(char(smoothParams.method), 'none')
    return;
end

breakMask = ~all(isfinite(xyIn), 2);

breakAtRepairs = false;
if isfield(smoothParams, 'break_at_repaired_points')
    breakAtRepairs = logical(smoothParams.break_at_repaired_points);
elseif isfield(smoothParams, 'break_at_repairs')
    breakAtRepairs = logical(smoothParams.break_at_repairs);
end

if breakAtRepairs
    repairMask = logical(repairMask(:));
    halfWidth = 1;
    if isfield(smoothParams, 'repair_break_halfwidth')
        halfWidth = max(0, round(smoothParams.repair_break_halfwidth));
    end
    repairBreak = dilate_mask_1d(repairMask, halfWidth);
    breakMask = breakMask | repairBreak;
end

runs = find_runs(~breakMask);
if isempty(runs)
    return;
end

for i = 1:size(runs,1)
    s = runs(i,1);
    e = runs(i,2);
    segLen = e - s + 1;
    if segLen < smoothParams.min_segment_length
        continue;
    end

    try
        xySmooth(s:e,1) = segmentwise_lowpass(xyIn(s:e,1), smoothParams.low_pass_frq, fps, smoothParams.n_pad);
        xySmooth(s:e,2) = segmentwise_lowpass(xyIn(s:e,2), smoothParams.low_pass_frq, fps, smoothParams.n_pad);
        smoothMeta.n_segments_smoothed = smoothMeta.n_segments_smoothed + 1;
        smoothMeta.segment_lengths(end+1,1) = segLen; 
    catch ME
        warning('smooth_track_2d_segments:SegmentFailed', ...
            'Skipping smoothing for segment [%d %d]: %s', s, e, ME.message);
    end
end
end
