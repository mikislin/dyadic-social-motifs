function [tracksOut, lowConfMask, thrByNode, info] = apply_confidence_mask(tracksIn, scores, confParams)
%APPLY_CONFIDENCE_MASK Mask low-confidence node samples as NaN.
% tracksIn: [T x nodes x 2]
% scores:   [T x nodes]

[T,nNodes,nCoords] = size(tracksIn);
assert(nCoords == 2, 'tracksIn must be [T x nodes x 2]');
assert(isequal(size(scores), [T nNodes]), 'scores must match [T x nodes]');

[thrByNode, info] = derive_confidence_thresholds(scores, confParams);
lowConfMask = false(T,nNodes);
tracksOut = tracksIn;

for node = 1:nNodes
    lowConfMask(:,node) = scores(:,node) < thrByNode(node);
    idx = lowConfMask(:,node);
    if any(idx)
        tracksOut(idx,node,:) = NaN;
    end
end
end
