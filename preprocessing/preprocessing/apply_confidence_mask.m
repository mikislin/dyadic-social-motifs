function [tracksOut, lowConfMask, thrByNode, info] = apply_confidence_mask(tracksIn, scores, confParams)
%APPLY_CONFIDENCE_MASK Mask low-confidence node samples as NaN.
% tracksIn: [T x nodes x 2]
% scores:   [T x nodes]

[T,nNodes,nCoords] = size(tracksIn);
assert(nCoords == 2, 'tracksIn must be [T x nodes x 2]');
assert(isequal(size(scores), [T nNodes]), 'scores must match [T x nodes]');

if nargin < 3 || isempty(confParams)
    P = default_preprocessing_params();
    confParams = P.confidence;
elseif isnumeric(confParams)
    threshold = confParams;
    P = default_preprocessing_params();
    confParams = P.confidence;
    if isscalar(threshold)
        confParams.mode = 'fixed';
        confParams.threshold = threshold;
    else
        assert(numel(threshold) >= nNodes, ...
            'Numeric confidence threshold vector must provide one value per node.');
        confParams.mode = 'fixed_by_node';
        confParams.threshold_by_node = threshold(:)';
    end
    confParams.max_mask_fraction = 1;
end

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
