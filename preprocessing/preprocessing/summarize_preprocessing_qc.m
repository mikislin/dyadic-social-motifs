function stats = summarize_preprocessing_qc(sessionPreproc)
%SUMMARIZE_PREPROCESSING_QC Session-level QC summaries.

badframes = sessionPreproc.qc.badframes;
qcFrame = sessionPreproc.qc.frames;
nAnimals = size(badframes,2);
stats = struct();
for m = 1:nAnimals
    stats.animal(m).pctBadframes = 100 * mean(badframes(:,m));
    stats.animal(m).pctInterpFrames = 100 * mean(qcFrame.fracInterp(:,m) > 0);
    stats.animal(m).pctJumpFrames = 100 * mean(qcFrame.fracJump(:,m) > 0);
    stats.animal(m).pctGeomFrames = 100 * mean(qcFrame.fracGeom(:,m) > 0);
    stats.animal(m).pctArenaFrames = 100 * mean(qcFrame.fracArena(:,m) > 0);
    stats.animal(m).pctLowConfFrames = 100 * mean(qcFrame.fracLowConf(:,m) > 0);
    stats.animal(m).pctJumpSamples = 100 * mean(qcFrame.fracJump(:,m));
    stats.animal(m).pctInterpSamples = 100 * mean(qcFrame.fracInterp(:,m));
    stats.animal(m).medianBodyLength = median(qcFrame.bodyLength(isfinite(qcFrame.bodyLength(:,m)),m));
    stats.animal(m).medianBodyLengthRaw = NaN;
    stats.animal(m).medianDistortionAnchorDispRatio = NaN;

    anchor = min(sessionPreproc.params.qc.arena_anchor_node, size(sessionPreproc.clean.tracks,2));
    tracksCleanAll = sessionPreproc.clean.tracks;
    if ndims(tracksCleanAll) == 3
        tracksCleanAll = reshape(tracksCleanAll, size(tracksCleanAll,1), size(tracksCleanAll,2), size(tracksCleanAll,3), 1);
    end

    if ~isempty(sessionPreproc.raw) && isfield(sessionPreproc.raw, 'SLEAPtracks') && ~isempty(sessionPreproc.raw.SLEAPtracks)
        tracksRawAll = sessionPreproc.raw.SLEAPtracks;
        if ndims(tracksRawAll) == 3
            tracksRawAll = reshape(tracksRawAll, size(tracksRawAll,1), size(tracksRawAll,2), size(tracksRawAll,3), 1);
        end

        tracksRaw = reshape(tracksRawAll(:,:,:,m), size(tracksRawAll,1), size(tracksRawAll,2), size(tracksRawAll,3));
        rawXY = reshape(tracksRaw(:, anchor, :), size(tracksRaw,1), size(tracksRaw,3));
        rawSpeed = sqrt(sum(diff(rawXY,1,1).^2,2));
        stats.animal(m).medianCentroidSpeedRaw = median(rawSpeed(isfinite(rawSpeed)));

        n1 = sessionPreproc.params.qc.body_length_nodes(1);
        n2 = sessionPreproc.params.qc.body_length_nodes(2);
        rawBL = sqrt(sum((squeeze(tracksRaw(:,n1,:)) - squeeze(tracksRaw(:,n2,:))).^2, 2));
        stats.animal(m).medianBodyLengthRaw = median(rawBL(isfinite(rawBL)));
    else
        stats.animal(m).medianCentroidSpeedRaw = NaN;
    end

    tracksClean = reshape(tracksCleanAll(:,:,:,m), size(tracksCleanAll,1), size(tracksCleanAll,2), size(tracksCleanAll,3));
    cleanXY = reshape(tracksClean(:, anchor, :), size(tracksClean,1), size(tracksClean,3));
    cleanSpeed = sqrt(sum(diff(cleanXY,1,1).^2,2));
    stats.animal(m).medianCentroidSpeedClean = median(cleanSpeed(isfinite(cleanSpeed)));

    if isfield(qcFrame, 'distortionAnchorDispRatio')
        vv = qcFrame.distortionAnchorDispRatio(:,m);
        stats.animal(m).medianDistortionAnchorDispRatio = median(vv(isfinite(vv)));
    end
end
end
