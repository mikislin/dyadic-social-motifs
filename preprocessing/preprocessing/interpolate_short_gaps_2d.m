function [xyOut, interpMask] = interpolate_short_gaps_2d(xyIn, time, gapParams)
%INTERPOLATE_SHORT_GAPS_2D Interpolate only short NaN runs in x/y together.

xyOut = xyIn;
time = time(:);
interpMask = false(size(xyIn,1),1);
missing = any(~isfinite(xyIn), 2);
% Early exit: no missing data → skip everything
if ~any(missing)
    return
end
runs = find_runs(missing);

for i = 1:size(runs,1)
    if isempty(runs)
        return
    end
    s = runs(i,1);
    e = runs(i,2);
    runLen = e - s + 1;
    if runLen > gapParams.max_gap_frames
        continue
    end
    if (~gapParams.fill_ends) && (s == 1 || e == size(xyIn,1))
        continue
    end

    goodMask = all(isfinite(xyOut),2);
    if nnz(goodMask) < gapParams.min_points_for_interp
        continue
    end

    tq = time(s:e);
    xq = interp1(time(goodMask), xyOut(goodMask,1), tq, gapParams.method, NaN);
    yq = interp1(time(goodMask), xyOut(goodMask,2), tq, gapParams.method, NaN);
    goodFill = isfinite(xq) & isfinite(yq);
    xyOut(s:e,1) = xq;
    xyOut(s:e,2) = yq;
    interpMask(s:e) = goodFill;
end
end
