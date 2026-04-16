function [thrByNode, info] = derive_confidence_thresholds(scores, confParams)
%DERIVE_CONFIDENCE_THRESHOLDS Compute per-node score thresholds.
% scores is [T x nodes].

scores = double(scores);
assert(ismatrix(scores), 'scores must be [T x nodes]');

nNodes = size(scores,2);
thrByNode = nan(1,nNodes);
info = struct();
info.maskFractionInitial = nan(1,nNodes);
info.maskFractionFinal = nan(1,nNodes);
info.wasRelaxed = false(1,nNodes);

for node = 1:nNodes
    s = scores(:,node);
    s = s(isfinite(s));

    switch lower(char(confParams.mode))
        case 'adaptive_per_node'
            if isempty(s)
                thr = confParams.min_threshold;
            else
                thr = prctile(s, confParams.prctile);
                thr = max(confParams.min_threshold, thr);
            end

        case 'fixed'
            thr = confParams.threshold;

        case 'fixed_by_node'
            assert(~isempty(confParams.threshold_by_node), ...
                'confParams.threshold_by_node is required for fixed_by_node mode');
            assert(numel(confParams.threshold_by_node) >= node, ...
                'confParams.threshold_by_node must provide one threshold per node');
            thr = confParams.threshold_by_node(node);

        otherwise
            error('Unknown confidence mode: %s', confParams.mode);
    end

    thrByNode(node) = thr;
    info.maskFractionInitial(node) = mean(scores(:,node) < thrByNode(node), 'omitnan');

    if info.maskFractionInitial(node) > confParams.max_mask_fraction
        thrByNode(node) = min(thrByNode(node), confParams.relaxed_threshold);
        info.wasRelaxed(node) = true;
    end

    info.maskFractionFinal(node) = mean(scores(:,node) < thrByNode(node), 'omitnan');
end
end
