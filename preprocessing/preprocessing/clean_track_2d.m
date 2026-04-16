function [xyClean, qc2d] = clean_track_2d(xyRaw, time, excludedFrames, P)
time = time(:);
excludedFrames = logical(excludedFrames(:));
xyRaw = double(xyRaw);

qc2d.nanMaskOriginal = any(~isfinite(xyRaw), 2);
qc2d.jumpMask = false(size(xyRaw,1),1);
qc2d.interpMask = false(size(xyRaw,1),1);
qc2d.finalNanMask = false(size(xyRaw,1),1);
qc2d.jumpMeta = struct('cutoff_quantile_px', NaN, 'cutoff_final_px', NaN, 'n_iter', 0, 'converged', true, 'nFlagged', 0);
qc2d.smoothMeta = struct('n_segments_smoothed', 0, 'segment_lengths', []);

xyWork = xyRaw;
validMask = ~excludedFrames;
if P.jump.enabled
    [xyWork, qc2d.jumpMask, qc2d.jumpMeta] = iterative_jump_filter_2d(xyWork, time, validMask, P.jump);
end
xyWork(qc2d.jumpMask,:) = NaN;

if P.gaps.enabled
    [xyWork, qc2d.interpMask] = interpolate_short_gaps_2d(xyWork, time, P.gaps);
end

qc2d.finalNanMask = any(~isfinite(xyWork),2);
xyClean = xyWork;
end
