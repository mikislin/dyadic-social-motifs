function [tracksClean, qcAnimal] = preprocess_animal_tracks(tracksRaw, time, excludedFrames, scores, P)
%PREPROCESS_ANIMAL_TRACKS Preprocess one animal's [T x nodes x 2] tracks.

if nargin < 4 || isempty(scores)
    scores = [];
end

[T, nNodes, nCoords] = size(tracksRaw);
assert(nCoords == 2, 'preprocess_animal_tracks expects [T x nodes x 2]');

tracksWork = double(tracksRaw);
qcAnimal = struct();
qcAnimal.confidenceThresholds = nan(1,nNodes);
qcAnimal.confidenceInfo = struct();

if P.confidence.enabled && ~isempty(scores)
    [tracksWork, lowConfMask, thrByNode, confInfo] = apply_confidence_mask(tracksWork, scores, P.confidence);
    qcAnimal.lowConfMask = sparse(lowConfMask);
    qcAnimal.confidenceThresholds = thrByNode;
    qcAnimal.confidenceInfo = confInfo;
else
    qcAnimal.lowConfMask = sparse(false(T, nNodes));
end

tracksClean = tracksWork;
qcAnimal.jumpMask = sparse(false(T, nNodes));
qcAnimal.interpMask = sparse(false(T, nNodes));
qcAnimal.geomMask = sparse(false(T, nNodes));
qcAnimal.arenaMask = sparse(false(T, nNodes));
qcAnimal.finalNanMask = sparse(false(T, nNodes));
qcAnimal.jumpMeta = cell(nNodes, 1);
qcAnimal.smoothMeta = cell(nNodes, 1);
qcAnimal.geometryInfo = struct();
qcAnimal.arena = struct();

for node = 1:nNodes
    xyRaw = reshape(tracksWork(:, node, :), T, 2);
    Pnode = P;
    if isfield(P.gaps, 'max_gap_frames_by_node') && numel(P.gaps.max_gap_frames_by_node) >= node
        Pnode.gaps.max_gap_frames = P.gaps.max_gap_frames_by_node(node);
    end
    if isfield(P.jump, 'max_disp_px_per_frame_by_node') && numel(P.jump.max_disp_px_per_frame_by_node) >= node
        Pnode.jump.max_disp_px_per_frame = P.jump.max_disp_px_per_frame_by_node(node);
    end
    if isfield(P.smooth, 'low_pass_frq_by_node') && numel(P.smooth.low_pass_frq_by_node) >= node
        Pnode.smooth.low_pass_frq = P.smooth.low_pass_frq_by_node(node);
    end

    [xyCleanNode, qc2d] = clean_track_2d(xyRaw, time, excludedFrames, Pnode);
    tracksClean(:, node, :) = xyCleanNode;
    qcAnimal.jumpMask(:, node) = sparse(qc2d.jumpMask);
    qcAnimal.interpMask(:, node) = sparse(qc2d.interpMask);
    qcAnimal.jumpMeta{node} = qc2d.jumpMeta;
    qcAnimal.smoothMeta{node} = qc2d.smoothMeta;
end

[geomMaskNode, geomInfo] = detect_geometry_outliers(tracksClean, P.qc.geometry_node_pairs, P.qc.geometry_prctile_range, scores);
qcAnimal.geometryInfo = geomInfo;
if any(geomMaskNode(:))
    for node = 1:nNodes
        idx = geomMaskNode(:,node);
        if any(idx)
            tracksClean(idx,node,:) = NaN;
        end
    end
    qcAnimal.geomMask = sparse(geomMaskNode);
end

anchorNode = min(P.qc.arena_anchor_node, nNodes);
anchorXY = reshape(tracksClean(:, anchorNode, :), T, 2);
validForArena = ~excludedFrames & ~geomMaskNode(:,anchorNode);
arena = estimate_session_arena(anchorXY, validForArena, P.qc.arena_percentile, P.qc.arena_margin_px);
qcAnimal.arena = arena;

for node = 1:nNodes
    nodeXY = reshape(tracksClean(:, node, :), T, 2);
    outside = outside_arena_mask(nodeXY, arena);
    if any(outside)
        tracksClean(outside, node, :) = NaN;
        qcAnimal.arenaMask(:,node) = sparse(outside);
    end
end

if P.gaps.enabled
    for node = 1:nNodes
        nodeXY = reshape(tracksClean(:, node, :), T, 2);
        gp = P.gaps;
        if isfield(P.gaps, 'max_gap_frames_by_node') && numel(P.gaps.max_gap_frames_by_node) >= node
            gp.max_gap_frames = P.gaps.max_gap_frames_by_node(node);
        end
        [nodeXY, interpMask2] = interpolate_short_gaps_2d(nodeXY, time, gp);
        tracksClean(:, node, :) = nodeXY;
        qcAnimal.interpMask(:, node) = sparse(full(qcAnimal.interpMask(:, node)) | interpMask2);
    end
end

if P.smooth.enabled
    for node = 1:nNodes
        nodeXY = reshape(tracksClean(:, node, :), T, 2);
        sp = P.smooth;
        if isfield(P.smooth, 'low_pass_frq_by_node') && numel(P.smooth.low_pass_frq_by_node) >= node
            sp.low_pass_frq = P.smooth.low_pass_frq_by_node(node);
        end
        repairMask = full(qcAnimal.jumpMask(:,node)) | full(qcAnimal.interpMask(:,node));
        [nodeXY, smoothMeta2] = smooth_track_2d_segments(nodeXY, sp, P.data.fps, repairMask);
        tracksClean(:, node, :) = nodeXY;
        qcAnimal.smoothMeta{node} = smoothMeta2;
    end
end

qcAnimal.finalNanMask = sparse(any(~isfinite(tracksClean), 3));
qcAnimal.nodeSummary = summarize_node_qc(qcAnimal);
end
