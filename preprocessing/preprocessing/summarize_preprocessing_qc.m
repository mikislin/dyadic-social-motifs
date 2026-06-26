function stats = summarize_preprocessing_qc(sessionPreproc)
%SUMMARIZE_PREPROCESSING_QC Session-level QC summaries.

badframes = sessionPreproc.qc.badframes;
qcFrame = sessionPreproc.qc.frames;
nAnimals = size(badframes,2);
predictionIssueByAnimal = false(size(badframes));
interpolatedIssueByAnimal = false(size(badframes));
usableIssueByAnimal = false(size(badframes));
repairedIssueByAnimal = false(size(badframes));
unresolvedIssueByAnimal = false(size(badframes));
stats = struct();
for m = 1:nAnimals
    predictionIssueFrame = qcFrame.fracLowConf(:,m) > 0 | ...
        qcFrame.fracJump(:,m) > 0 | ...
        qcFrame.fracGeom(:,m) > 0;
    interpolatedIssueFrame = predictionIssueFrame & qcFrame.fracInterp(:,m) > 0;
    finalNanFrame = qcFrame.fracFinalNaN(:,m) > 0;
    usableIssueFrame = predictionIssueFrame & ~badframes(:,m) & ~finalNanFrame;
    repairedIssueFrame = interpolatedIssueFrame & usableIssueFrame;
    unresolvedIssueFrame = predictionIssueFrame & (badframes(:,m) | finalNanFrame);

    predictionIssueByAnimal(:,m) = predictionIssueFrame;
    interpolatedIssueByAnimal(:,m) = interpolatedIssueFrame;
    usableIssueByAnimal(:,m) = usableIssueFrame;
    repairedIssueByAnimal(:,m) = repairedIssueFrame;
    unresolvedIssueByAnimal(:,m) = unresolvedIssueFrame;

    stats.animal(m).pctBadframes = 100 * mean(badframes(:,m));
    stats.animal(m).pctInterpFrames = 100 * mean(qcFrame.fracInterp(:,m) > 0);
    stats.animal(m).pctJumpFrames = 100 * mean(qcFrame.fracJump(:,m) > 0);
    stats.animal(m).pctGeomFrames = 100 * mean(qcFrame.fracGeom(:,m) > 0);
    stats.animal(m).pctArenaFrames = 100 * mean(qcFrame.fracArena(:,m) > 0);
    stats.animal(m).pctLowConfFrames = 100 * mean(qcFrame.fracLowConf(:,m) > 0);
    stats.animal(m).pctJumpSamples = 100 * mean(qcFrame.fracJump(:,m));
    stats.animal(m).pctInterpSamples = 100 * mean(qcFrame.fracInterp(:,m));
    stats.animal(m).nPredictionIssueFrames = nnz(predictionIssueFrame);
    stats.animal(m).nInterpolatedPredictionIssueFrames = nnz(interpolatedIssueFrame);
    stats.animal(m).nUsablePredictionIssueFrames = nnz(usableIssueFrame);
    stats.animal(m).nRepairedPredictionIssueFrames = nnz(repairedIssueFrame);
    stats.animal(m).nUnresolvedPredictionIssueFrames = nnz(unresolvedIssueFrame);
    stats.animal(m).pctPredictionIssueFrames = 100 * mean(predictionIssueFrame);
    stats.animal(m).pctPredictionIssueUsable = local_percent(nnz(usableIssueFrame), nnz(predictionIssueFrame));
    stats.animal(m).pctPredictionIssueRepaired = local_percent(nnz(repairedIssueFrame), nnz(predictionIssueFrame));
    stats.animal(m).pctPredictionIssueUnresolved = local_percent(nnz(unresolvedIssueFrame), nnz(predictionIssueFrame));
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

stats.badframeFraction = mean(badframes(:));
stats.nPredictionIssueAnimalFrames = nnz(predictionIssueByAnimal);
stats.nInterpolatedPredictionIssueAnimalFrames = nnz(interpolatedIssueByAnimal);
stats.nUsablePredictionIssueAnimalFrames = nnz(usableIssueByAnimal);
stats.nRepairedPredictionIssueAnimalFrames = nnz(repairedIssueByAnimal);
stats.nUnresolvedPredictionIssueAnimalFrames = nnz(unresolvedIssueByAnimal);
stats.predictionIssueFraction = mean(predictionIssueByAnimal(:));
stats.interpolatedPredictionIssueFraction = mean(interpolatedIssueByAnimal(:));
stats.usablePredictionIssueFraction = mean(usableIssueByAnimal(:));
stats.repairedPredictionIssueFraction = mean(repairedIssueByAnimal(:));
stats.unresolvedPredictionIssueFraction = mean(unresolvedIssueByAnimal(:));
stats.predictionIssueUsableRate = local_fraction(nnz(usableIssueByAnimal), nnz(predictionIssueByAnimal));
stats.predictionIssueRepairRate = local_fraction(nnz(repairedIssueByAnimal), nnz(predictionIssueByAnimal));
stats.predictionIssueUnresolvedRate = local_fraction(nnz(unresolvedIssueByAnimal), nnz(predictionIssueByAnimal));

stats.nPredictionIssueFramesAnyAnimal = nnz(any(predictionIssueByAnimal, 2));
stats.nRepairedPredictionIssueFramesAnyAnimal = nnz(any(repairedIssueByAnimal, 2));
stats.nUnresolvedPredictionIssueFramesAnyAnimal = nnz(any(unresolvedIssueByAnimal, 2));
stats.predictionIssueAnyAnimalFraction = mean(any(predictionIssueByAnimal, 2));
stats.repairedPredictionIssueAnyAnimalFraction = mean(any(repairedIssueByAnimal, 2));
stats.unresolvedPredictionIssueAnyAnimalFraction = mean(any(unresolvedIssueByAnimal, 2));
end

function p = local_percent(num, den)
p = 100 * local_fraction(num, den);
end

function f = local_fraction(num, den)
if den > 0
    f = num ./ den;
else
    f = NaN;
end
end
