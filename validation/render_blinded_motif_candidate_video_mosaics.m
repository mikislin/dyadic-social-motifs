function RenderAudit = render_blinded_motif_candidate_video_mosaics( ...
        repoRoot, Selection, CandidateSummary, outRoot, params)
%RENDER_BLINDED_MOTIF_CANDIDATE_VIDEO_MOSAICS Render neutral 3x3 pose clips.
%
% Rendered tiles expose only a frozen candidate code, a neutral tile letter,
% and relative time. Selection stratum, session, arena, annotation, and all
% experimental metadata remain absent from the videos and contact sheets.

repoRoot = string(repoRoot);
outRoot = string(outRoot);
videoRoot = fullfile(outRoot, params.video_directory);
contactRoot = fullfile(outRoot, params.contact_sheet_directory);
i_make_dir(videoRoot);
i_make_dir(contactRoot);

[nodeMap, ~] = default_sleap_node_map();
edges = i_skeleton_edges(nodeMap);
excludedNodes = zeros(0, 1);
if params.exclude_ear_nodes
    excludedNodes = [nodeMap.leftEar; nodeMap.rightEar];
end

RenderAudit = table();
for iCandidate = 1:height(CandidateSummary)
    candidateId = string(CandidateSummary.motif_candidate_id(iCandidate));
    examples = Selection( ...
        string(Selection.motif_candidate_id) == candidateId, :);
    examples = sortrows(examples, 'display_tile_index');
    status = "success";
    errorIdentifier = "";
    errorMessage = "";
    videoRendered = false;
    outputFrameCount = 0;
    videoPath = fullfile(videoRoot, ...
        i_short_candidate_id(candidateId) + "_blinded_3x3.mp4");
    contactPath = fullfile(contactRoot, ...
        i_short_candidate_id(candidateId) + ...
        "_blinded_anchor_contact_sheet.png");
    fig = [];
    writer = [];
    try
        assert(height(examples) == params.tiles_per_candidate && ...
            numel(unique(examples.graph_node_id)) == ...
            params.tiles_per_candidate, ...
            'render_blinded_motif_candidate_video_mosaics:BadSelection', ...
            'Candidate %s does not have nine unique examples.', candidateId);
        [clips, relativeTimeSec] = i_load_clips( ...
            repoRoot, examples, params);
        [fig, axesHandles, mouseHandles, titleHandle] = ...
            i_initialize_figure(examples, candidateId, edges, ...
            excludedNodes, params);
        anchorIndex = find(abs(relativeTimeSec) == ...
            min(abs(relativeTimeSec)), 1);
        i_update_figure(clips, anchorIndex, relativeTimeSec(anchorIndex), ...
            axesHandles, mouseHandles, titleHandle, edges, ...
            excludedNodes, candidateId);
        exportgraphics(fig, contactPath, ...
            'Resolution', params.contact_sheet_dpi);

        if params.render_videos
            writer = VideoWriter(videoPath, char(params.video_profile));
            writer.FrameRate = params.output_frame_rate_hz;
            open(writer);
            for iFrame = 1:numel(relativeTimeSec)
                i_update_figure(clips, iFrame, relativeTimeSec(iFrame), ...
                    axesHandles, mouseHandles, titleHandle, edges, ...
                    excludedNodes, candidateId);
                writeVideo(writer, getframe(fig));
            end
            outputFrameCount = numel(relativeTimeSec);
            close(writer);
            videoRendered = true;
        else
            videoPath = "";
        end
        close(fig);
    catch ME
        status = "failed";
        errorIdentifier = string(ME.identifier);
        errorMessage = string(ME.message);
        if ~isempty(writer)
            try
                close(writer);
            catch
            end
        end
        if ~isempty(fig) && isgraphics(fig)
            close(fig);
        end
    end

    [videoHash, videoBytes] = i_hash_if_present(videoPath);
    [contactHash, contactBytes] = i_hash_if_present(contactPath);
    row = table(candidateId, ...
        CandidateSummary.candidate_local_index(iCandidate), ...
        height(examples), numel(unique(examples.session_index)), ...
        CandidateSummary.clip_half_window_sec(iCandidate), ...
        params.output_frame_rate_hz, outputFrameCount, ...
        2 * CandidateSummary.clip_half_window_sec(iCandidate), ...
        videoRendered, string(videoPath), videoBytes, videoHash, ...
        string(contactPath), contactBytes, contactHash, status, ...
        errorIdentifier, errorMessage, ...
        params.expected_membership_sha256, ...
        params.expected_candidate_freeze_id, "none", ...
        'VariableNames', {'motif_candidate_id', ...
        'candidate_local_index', 'selected_example_count', ...
        'selected_session_count', 'clip_half_window_sec', ...
        'output_frame_rate_hz', 'output_frame_count', ...
        'nominal_clip_duration_sec', 'video_rendered', 'video_file', ...
        'video_file_bytes', 'video_sha256', 'contact_sheet_file', ...
        'contact_sheet_file_bytes', 'contact_sheet_sha256', ...
        'render_status', 'render_error_identifier', ...
        'render_error_message', 'frozen_membership_sha256', ...
        'candidate_freeze_id', 'experimental_labels_used'});
    RenderAudit = [RenderAudit; row]; %#ok<AGROW>
end
end

function [clips, relativeTimeSec] = i_load_clips(repoRoot, examples, params)
halfWindowSec = examples.clip_half_window_sec(1);
outputFrameCount = round(2 * halfWindowSec * ...
    params.output_frame_rate_hz) + 1;
relativeTimeSec = linspace(-halfWindowSec, halfWindowSec, ...
    outputFrameCount)';
relativeFrameOffset = round(relativeTimeSec * ...
    params.input_frame_rate_hz);
clips = cell(height(examples), 1);
for i = 1:height(examples)
    pathText = resolve_repo_path(repoRoot, ...
        examples.preprocess_output_file(i));
    loaded = load(pathText, 'sessionPreproc');
    tracks = loaded.sessionPreproc.clean.tracks;
    assert(ndims(tracks) == 4 && size(tracks, 4) >= 2, ...
        'render_blinded_motif_candidate_video_mosaics:BadTrackShape', ...
        'Pose file must contain frame-by-node-by-coordinate-by-animal tracks.');
    assert(size(tracks, 1) == examples.source_frame_count(i), ...
        'render_blinded_motif_candidate_video_mosaics:FrameCountMismatch', ...
        'Pose and preprocessing-audit frame counts differ.');
    anchorFrame = round(examples.anchor_frame(i));
    sourceFrames = min(max(anchorFrame + relativeFrameOffset, 1), ...
        size(tracks, 1));
    clip = struct();
    clip.tracks = tracks(sourceFrames, :, :, 1:2);
    clip.source_frames = sourceFrames;
    clips{i} = clip;
    clear loaded tracks
end
end

function [fig, axesHandles, mouseHandles, titleHandle] = ...
        i_initialize_figure(examples, candidateId, edges, ...
        excludedNodes, params)
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
    'Position', [50 50 params.figure_pixels params.figure_pixels]);
layout = tiledlayout(fig, 3, 3, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
axesHandles = gobjects(params.tiles_per_candidate, 1);
mouseHandles = cell(params.tiles_per_candidate, 2);
colors = [0.00 0.45 0.70; 0.80 0.40 0.00];
for iTile = 1:params.tiles_per_candidate
    ax = nexttile(layout, iTile);
    axesHandles(iTile) = ax;
    hold(ax, 'on');
    for iMouse = 1:2
        lineHandle = plot(ax, nan, nan, '-', ...
            'Color', colors(iMouse, :), 'LineWidth', 2);
        pointHandle = scatter(ax, nan, nan, 20, ...
            colors(iMouse, :), 'filled', ...
            'MarkerEdgeColor', [0.95 0.95 0.95], ...
            'LineWidth', 0.4);
        mouseHandles{iTile, iMouse} = struct( ...
            'line', lineHandle, 'points', pointHandle);
    end
    hold(ax, 'off');
    axis(ax, 'equal');
    xlim(ax, params.camera_x_limits);
    ylim(ax, params.camera_y_limits);
    set(ax, 'YDir', 'reverse', 'XTick', [], 'YTick', [], ...
        'Box', 'on', 'Color', [0.98 0.98 0.98]);
    title(ax, examples.display_tile_label(iTile), ...
        'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 11);
end
titleHandle = sgtitle(layout, ...
    i_short_candidate_id(candidateId) + " | blinded candidate examples", ...
    'Interpreter', 'none', 'FontWeight', 'normal');

% Verify helper inputs during figure construction, before video encoding.
assert(~isempty(edges) && all(excludedNodes >= 1), ...
    'render_blinded_motif_candidate_video_mosaics:BadSkeletonMetadata', ...
    'Skeleton edges and excluded-node indices must be valid.');
end

function i_update_figure(clips, frameIndex, relativeTimeSec, ...
        axesHandles, mouseHandles, titleHandle, edges, ...
        excludedNodes, candidateId)
for iTile = 1:numel(clips)
    pose = squeeze(clips{iTile}.tracks(frameIndex, :, :, :));
    for iMouse = 1:2
        xy = squeeze(pose(:, :, iMouse));
        [edgeX, edgeY] = i_edge_coordinates(xy, edges);
        set(mouseHandles{iTile, iMouse}.line, ...
            'XData', edgeX, 'YData', edgeY);
        finiteNode = all(isfinite(xy), 2);
        finiteNode(excludedNodes(excludedNodes <= numel(finiteNode))) = false;
        set(mouseHandles{iTile, iMouse}.points, ...
            'XData', xy(finiteNode, 1), ...
            'YData', xy(finiteNode, 2));
    end
    if frameIndex == 1
        axesHandles(iTile).Toolbar.Visible = 'off';
    end
end
titleHandle.String = sprintf('%s | blinded candidate examples | t = %+.2f s', ...
    i_short_candidate_id(candidateId), relativeTimeSec);
drawnow;
end

function [x, y] = i_edge_coordinates(xy, edges)
nEdges = size(edges, 1);
x = nan(3 * nEdges, 1);
y = nan(3 * nEdges, 1);
for iEdge = 1:nEdges
    nodes = edges(iEdge, :);
    if all(nodes <= size(xy, 1)) && all(all(isfinite(xy(nodes, :))))
        loc = 3 * (iEdge - 1) + (1:2);
        x(loc) = xy(nodes, 1);
        y(loc) = xy(nodes, 2);
    end
end
end

function edges = i_skeleton_edges(nodeMap)
edges = [ ...
    nodeMap.nose, nodeMap.neck
    nodeMap.neck, nodeMap.midBody
    nodeMap.midBody, nodeMap.body
    nodeMap.body, nodeMap.tailBase
    nodeMap.tailBase, nodeMap.tailMid
    nodeMap.tailMid, nodeMap.tailTip
    nodeMap.midBody, nodeMap.lfPaw
    nodeMap.midBody, nodeMap.rfPaw
    nodeMap.body, nodeMap.lhPaw
    nodeMap.body, nodeMap.rhPaw];
end

function value = i_short_candidate_id(candidateId)
value = regexprep(string(candidateId), '^MC_M[0-9]+_', '');
end

function [hash, bytes] = i_hash_if_present(pathText)
if strlength(string(pathText)) > 0 && isfile(pathText)
    [hash, bytes] = compute_file_sha256(pathText);
else
    hash = "";
    bytes = NaN;
end
end

function i_make_dir(pathText)
if ~isfolder(pathText)
    mkdir(pathText);
end
end
