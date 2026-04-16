function [badframes, qcFrame] = make_badframes(tracksClean, excludedFrames, qcAnimals, P)
%MAKE_BADFRAMES Build final badframe mask after preprocessing.

[T, ~, ~, nAnimals] = size(tracksClean);
excludedFrames = logical(excludedFrames(:));
badframes = false(T, nAnimals);
qcFrame = struct();
qcFrame.fracFinalNaN = zeros(T, nAnimals);
qcFrame.fracInterp = zeros(T, nAnimals);
qcFrame.fracJump = zeros(T, nAnimals);
qcFrame.fracGeom = zeros(T, nAnimals);
qcFrame.fracArena = zeros(T, nAnimals);
qcFrame.fracLowConf = zeros(T, nAnimals);
qcFrame.bodyLength = nan(T, nAnimals);
qcFrame.bodyLengthOutlier = false(T, nAnimals);
qcFrame.distortionAnchorDispRatio = nan(T, nAnimals);

expandSevere = 0;
if isfield(P.qc, 'expand_severe_badframes') && ~isempty(P.qc.expand_severe_badframes)
    expandSevere = P.qc.expand_severe_badframes;
end

for m = 1:nAnimals
    finalNanMask = full(qcAnimals(m).finalNanMask);
    interpMask = full(qcAnimals(m).interpMask);
    jumpMask = full(qcAnimals(m).jumpMask);
    geomMask = full(qcAnimals(m).geomMask);
    arenaMask = full(qcAnimals(m).arenaMask);
    lowConfMask = full(qcAnimals(m).lowConfMask);

    qcFrame.fracFinalNaN(:,m) = mean(finalNanMask, 2);
    qcFrame.fracInterp(:,m) = mean(interpMask, 2);
    qcFrame.fracJump(:,m) = mean(jumpMask, 2);
    qcFrame.fracGeom(:,m) = mean(geomMask, 2);
    qcFrame.fracArena(:,m) = mean(arenaMask, 2);
    qcFrame.fracLowConf(:,m) = mean(lowConfMask, 2);

    n1 = min(P.qc.body_length_nodes(1), size(tracksClean,2));
    n2 = min(P.qc.body_length_nodes(2), size(tracksClean,2));
    xy1 = squeeze(tracksClean(:, n1, :, m));
    xy2 = squeeze(tracksClean(:, n2, :, m));
    qcFrame.bodyLength(:,m) = sqrt(sum((xy1 - xy2).^2, 2));

    anchor = min(P.qc.arena_anchor_node, size(tracksClean,2));
    anchorXY = squeeze(tracksClean(:, anchor, :, m));
    anchorDisp = [NaN; sqrt(sum(diff(anchorXY,1,1).^2, 2))];
    bodyLen = qcFrame.bodyLength(:,m);
    qcFrame.distortionAnchorDispRatio(:,m) = anchorDisp ./ max(bodyLen, eps);

    finiteBL = isfinite(bodyLen);
    if nnz(finiteBL) >= 10
        pr = prctile(bodyLen(finiteBL), P.qc.body_length_prctile_range);
        qcFrame.bodyLengthOutlier(:,m) = bodyLen < pr(1) | bodyLen > pr(2);
    end

    badframes(:,m) = excludedFrames | ...
        (qcFrame.fracFinalNaN(:,m) > P.qc.frame_bad_node_frac_thresh) | ...
        (qcFrame.fracInterp(:,m) > P.qc.max_interp_frac_per_frame) | ...
        (qcFrame.fracJump(:,m) > P.qc.max_jump_frac_per_frame) | ...
        (qcFrame.fracGeom(:,m) > P.qc.max_geom_frac_per_frame) | ...
        (qcFrame.fracArena(:,m) > P.qc.max_arena_frac_per_frame) | ...
        (qcFrame.fracLowConf(:,m) > P.qc.max_lowconf_frac_per_frame) | ...
        qcFrame.bodyLengthOutlier(:,m) | ...
        (P.qc.require_finite_body_length & ~finiteBL);

    if expandSevere > 0
        severe = excludedFrames | (qcFrame.fracArena(:,m) > 0) | (qcFrame.fracGeom(:,m) > 0);
        badframes(:,m) = badframes(:,m) | dilate_mask_1d(severe, expandSevere);
    end
end

qcFrame.badframes_pair = any(badframes, 2);
end
