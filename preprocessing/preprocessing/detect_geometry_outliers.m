function [geomMaskNode, info] = detect_geometry_outliers(tracksClean, nodePairs, prctileRange, scores)
% Detect geometry outliers and assign blame to one node in each offending pair.
[T,nNodes,~] = size(tracksClean);
geomMaskNode = false(T,nNodes);
info = struct();
info.nodePairs = nodePairs;
info.thresholds = nan(size(nodePairs,1),2);

if nargin < 4
    scores = [];
end

for p = 1:size(nodePairs,1)
    n1 = nodePairs(p,1); n2 = nodePairs(p,2);
    xy1 = squeeze(tracksClean(:,n1,:));
    xy2 = squeeze(tracksClean(:,n2,:));
    dd = sqrt(sum((xy1-xy2).^2,2));
    valid = isfinite(dd);
    if nnz(valid) < 20, continue; end
    lims = prctile(dd(valid), prctileRange);
    info.thresholds(p,:) = lims;
    bad = valid & (dd < lims(1) | dd > lims(2));
    if ~any(bad), continue; end

    d1 = [0; sqrt(sum(diff(xy1).^2,2))];
    d2 = [0; sqrt(sum(diff(xy2).^2,2))];
    for t = find(bad)'
        blame = n1;
        if ~isempty(scores)
            s1 = scores(t, min(n1,size(scores,2)));
            s2 = scores(t, min(n2,size(scores,2)));
            if isfinite(s1) && isfinite(s2)
                if s2 < s1, blame = n2; end
                if s1 < s2, blame = n1; end
            elseif isfinite(s2) && ~isfinite(s1)
                blame = n1;
            elseif isfinite(s1) && ~isfinite(s2)
                blame = n2;
            end
            if s1 == s2
                if d2(t) > d1(t), blame = n2; else, blame = n1; end
            end
        else
            if d2(t) > d1(t), blame = n2; else, blame = n1; end
        end
        geomMaskNode(t, blame) = true;
    end
end
end
