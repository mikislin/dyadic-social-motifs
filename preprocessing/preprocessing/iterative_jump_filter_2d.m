function [xyClean, jumpMask, meta] = iterative_jump_filter_2d(xyRaw, time, validMask, jumpParams)
%ITERATIVE_JUMP_FILTER_2D Detect 2D teleports and replace them by interpolation.

xyClean = double(xyRaw);
time = time(:);
validMask = logical(validMask(:));
T = size(xyRaw,1);
jumpMask = false(T,1);
meta = struct('cutoff_quantile_px', NaN, 'cutoff_final_px', NaN, 'n_iter', 0, ...
    'converged', true, 'nFlagged', 0);

finiteMask = all(isfinite(xyRaw),2) & validMask;
if nnz(finiteMask) < jumpParams.min_valid_points
    meta.converged = false;
    return
end

baseDisp = sqrt(sum(diff(xyRaw,1,1).^2, 2));
baseDisp = baseDisp(isfinite(baseDisp));
if isempty(baseDisp)
    return
end

cutoffQuant = quantile(baseDisp, jumpParams.diff_quantile);
cutoff = min(cutoffQuant, jumpParams.max_disp_px_per_frame);
meta.cutoff_quantile_px = cutoffQuant;
meta.cutoff_final_px = cutoff;

for k = 1:jumpParams.max_iter
    dispNow = sqrt(sum(diff(xyClean,1,1).^2, 2));
    badIdx = find(dispNow > cutoff) + 1;
    badIdx = badIdx(validMask(badIdx));
    badIdx = badIdx(all(isfinite(xyClean(badIdx,:)),2));

    if isempty(badIdx)
        meta.n_iter = k - 1;
        meta.nFlagged = nnz(jumpMask);
        return
    end

    if jumpParams.dilate_frames > 0
        localMask = false(T,1);
        localMask(badIdx) = true;
        localMask = dilate_mask_1d(localMask, jumpParams.dilate_frames);
        badIdx = find(localMask & validMask);
    end

    goodMask = all(isfinite(xyClean),2) & validMask;
    goodMask(badIdx) = false;
    if nnz(goodMask) < 2
        meta.n_iter = k;
        meta.converged = false;
        meta.nFlagged = nnz(jumpMask);
        return
    end

    xyClean(badIdx,1) = interp1(time(goodMask), xyClean(goodMask,1), time(badIdx), jumpParams.method, NaN);
    xyClean(badIdx,2) = interp1(time(goodMask), xyClean(goodMask,2), time(badIdx), jumpParams.method, NaN);
    jumpMask(badIdx) = true;
    meta.n_iter = k;
end

meta.converged = false;
meta.nFlagged = nnz(jumpMask);
end
