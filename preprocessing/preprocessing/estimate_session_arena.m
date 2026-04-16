function arena = estimate_session_arena(anchorXY, validMask, pct, marginPx)
%ESTIMATE_SESSION_ARENA Estimate session-specific arena bounds from anchor positions.

if nargin < 2 || isempty(validMask)
    validMask = true(size(anchorXY,1),1);
end
validMask = logical(validMask(:));
xy = anchorXY(validMask,:);
xy = xy(all(isfinite(xy),2),:);

arena = struct('xlim', [-Inf Inf], 'ylim', [-Inf Inf], 'isValid', false);
if size(xy,1) < 50
    return
end

arena.xlim = prctile(xy(:,1), pct) + [-marginPx marginPx];
arena.ylim = prctile(xy(:,2), pct) + [-marginPx marginPx];
arena.isValid = true;
end
