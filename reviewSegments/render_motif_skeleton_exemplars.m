function [fileTable, notes] = render_motif_skeleton_exemplars(PoseData, exemplarTable, varargin)
% render_motif_skeleton_exemplars
% Export skeleton-only dyadic motif exemplar clips without raw video.
%
% Robust to MATLAB installations where MPEG-4 VideoWriter is unavailable.
% Falls back to Motion JPEG AVI, then Uncompressed AVI.
%
% Inputs
%   PoseData       struct, DyadSet-like struct, cell, or numeric track array
%   exemplarTable  table returned by find_motif_video_exemplars
%
% Name-value options
%   OutputDir      directory for clips
%   FrameRate      frames/s for output video
%   WindowSec      clip half-window when table does not include start/end frames
%   Edges          Nedge x 2 skeleton edges, default mouse-like 13-node graph
%   NodeNames      optional node labels
%   Profile        preferred VideoWriter profile: 'MPEG-4', 'Motion JPEG AVI', etc.
%   FileFormat     'auto', 'mp4', or 'avi'
%   Visible        true/false figure visibility
%   MarkerSize     skeleton marker size
%   LineWidth      skeleton edge line width
%   ArenaPadding   plot padding in coordinate units
%   FlipY          true/false; use true if video-coordinate y axis should point down
%
% Output
%   fileTable      table with motif/exemplar rows and clip paths
%   notes          struct with writer profile and render diagnostics

p = inputParser;
p.addParameter('OutputDir', fullfile(pwd, 'motif_skeleton_exemplars'), @(x)ischar(x)||isstring(x));
p.addParameter('FrameRate', 80, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('WindowSec', 1.0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('Edges', [], @(x) isempty(x) || (isnumeric(x) && size(x,2)==2));
% Compatibility aliases used by find_motif_video_exemplars.
p.addParameter('SkeletonEdges', [], @(x) isempty(x) || (isnumeric(x) && size(x,2)==2));
p.addParameter('NodeNames', string.empty(1,0));
p.addParameter('BodypartNames', string.empty(1,0));
p.addParameter('Profile', 'auto', @(x)ischar(x)||isstring(x));
p.addParameter('VideoProfile', 'auto', @(x)ischar(x)||isstring(x));
p.addParameter('FileFormat', 'auto', @(x)ischar(x)||isstring(x));
p.addParameter('VideoFormat', 'auto', @(x)ischar(x)||isstring(x));
p.addParameter('Visible', false, @(x)islogical(x)||isnumeric(x));
p.addParameter('MarkerSize', 28, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('LineWidth', 2, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('ArenaPadding', 25, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('FlipY', true, @(x)islogical(x)||isnumeric(x));
p.addParameter('MaxFramesPerClip', 240, @(x)isnumeric(x)&&isscalar(x)&&x>=10);
p.addParameter('Verbose', true, @(x)islogical(x)||isnumeric(x));
p.parse(varargin{:});
P = p.Results;
P.OutputDir = char(P.OutputDir);
P.Profile = char(P.Profile);
P.VideoProfile = char(P.VideoProfile);
P.FileFormat = char(P.FileFormat);
P.VideoFormat = char(P.VideoFormat);
P.Verbose = logical(P.Verbose);

% Honor aliases without requiring callers to know the renderer's internal names.
if isempty(P.Edges) && ~isempty(P.SkeletonEdges)
    P.Edges = P.SkeletonEdges;
end
if isempty(P.NodeNames) && ~isempty(P.BodypartNames)
    P.NodeNames = P.BodypartNames;
end
if strcmpi(P.Profile, 'auto') && ~strcmpi(P.VideoProfile, 'auto')
    P.Profile = P.VideoProfile;
end
if strcmpi(P.FileFormat, 'auto') && ~strcmpi(P.VideoFormat, 'auto')
    P.FileFormat = P.VideoFormat;
end

if ~exist(P.OutputDir, 'dir')
    mkdir(P.OutputDir);
end

if isempty(P.Edges)
    nPartsGuess = local_guess_num_bodyparts(PoseData);
    P.Edges = local_default_mouse_edges(nPartsGuess);
end

tracksBySession = local_extract_tracks_by_session(PoseData);
if isempty(tracksBySession)
    error('render_motif_skeleton_exemplars:NoTracks', ...
        'Could not find pose tracks in PoseData. Expected tracks/SLEAPtracks/dyadCell or numeric [frames x nodes x xy x mice].');
end

nRows = height(exemplarTable);
outMotif = nan(nRows,1);
outRank = nan(nRows,1);
outSession = nan(nRows,1);
outStartFrame = nan(nRows,1);
outEndFrame = nan(nRows,1);
outAnchorFrame = nan(nRows,1);
outFile = strings(nRows,1);
outProfile = strings(nRows,1);
outStatus = strings(nRows,1);

writerInfo = local_choose_writer(P.Profile, P.FileFormat);

for r = 1:nRows
    row = exemplarTable(r,:);
    motifId = local_get_scalar(row, {'motif_id','motif','cluster_id'}, NaN);
    rankId = local_get_scalar(row, {'rank','exemplar_rank','neighbor_rank'}, r);
    sess = local_get_scalar(row, {'session_index','session'}, 1);
    anchorFrame = round(local_get_scalar(row, {'anchor_frame','medoid_anchor_frame','center_frame'}, NaN));

    startFrame = round(local_get_scalar(row, {'start_frame','clip_start_frame','segment_start_frame'}, NaN));
    endFrame = round(local_get_scalar(row, {'end_frame','clip_end_frame','segment_end_frame'}, NaN));
    if isnan(startFrame) || isnan(endFrame)
        if isnan(anchorFrame)
            startFrame = round(local_get_scalar(row, {'segment_start_frame','start_anchor_frame'}, 1));
            endFrame = round(local_get_scalar(row, {'segment_end_frame','end_anchor_frame'}, startFrame));
            anchorFrame = round((startFrame + endFrame) / 2);
        else
            halfWin = round(P.WindowSec * P.FrameRate);
            startFrame = anchorFrame - halfWin;
            endFrame = anchorFrame + halfWin;
        end
    end

    if sess < 1 || sess > numel(tracksBySession) || isempty(tracksBySession{sess})
        outStatus(r) = "missing_session_tracks";
        continue
    end

    tracks = tracksBySession{sess};
    nFrames = size(tracks,1);
    startFrame = max(1, min(nFrames, startFrame));
    endFrame = max(1, min(nFrames, endFrame));
    if endFrame < startFrame
        tmp = startFrame; startFrame = endFrame; endFrame = tmp;
    end
    frameIdx = startFrame:endFrame;
    if numel(frameIdx) > P.MaxFramesPerClip
        keep = round(linspace(1, numel(frameIdx), P.MaxFramesPerClip));
        frameIdx = frameIdx(keep);
    end

    baseName = sprintf('motif%02d_rank%02d_session%02d_frame%06d', motifId, rankId, sess, anchorFrame);
    clipFile = fullfile(P.OutputDir, [baseName writerInfo.extension]);

    try
        actualProfile = local_write_skeleton_clip(clipFile, tracks, frameIdx, P.Edges, row, P, writerInfo);
        outStatus(r) = "ok";
        outProfile(r) = string(actualProfile);
        outFile(r) = string(clipFile);
    catch ME
        outStatus(r) = "failed: " + string(ME.message);
        outFile(r) = string(clipFile);
    end

    outMotif(r) = motifId;
    outRank(r) = rankId;
    outSession(r) = sess;
    outStartFrame(r) = startFrame;
    outEndFrame(r) = endFrame;
    outAnchorFrame(r) = anchorFrame;
end

exemplarId = compose('motif%02d_rank%02d', outMotif, outRank);
clip_file = outFile;
fileTable = table(exemplarId, outMotif, outRank, outSession, outAnchorFrame, outStartFrame, outEndFrame, ...
    clip_file, outProfile, outStatus, ...
    'VariableNames', {'exemplar_id','motif_id','rank','session_index','anchor_frame','start_frame','end_frame','clip_file','writer_profile','status'});

notes = strings(0,1);
notes(end+1,1) = sprintf('render_motif_skeleton_exemplars | outputDir = %s', P.OutputDir);
notes(end+1,1) = sprintf('render_motif_skeleton_exemplars | requested profile = %s', P.Profile);
notes(end+1,1) = sprintf('render_motif_skeleton_exemplars | requested = %d | written = %d | first profile = %s', ...
    nRows, sum(outStatus == "ok"), writerInfo.profile);

if P.Verbose
    fprintf('render_motif_skeleton_exemplars | requested = %d | written = %d | first profile = %s\n', ...
        nRows, sum(outStatus == "ok"), writerInfo.profile);
end
end

function actualProfile = local_write_skeleton_clip(clipFile, tracks, frameIdx, edges, row, P, writerInfo)
profileList = writerInfo.profileList;
lastErr = [];
for ip = 1:numel(profileList)
    profile = profileList{ip};
    try
        [filepath, nameOnly, ~] = fileparts(clipFile);
        ext = local_profile_extension(profile);
        thisFile = fullfile(filepath, [nameOnly ext]);
        vw = VideoWriter(thisFile, profile);
        vw.FrameRate = P.FrameRate;
        if isprop(vw, 'Quality')
            vw.Quality = 90;
        end
        open(vw);
        cleanupObj = onCleanup(@() local_safe_close(vw)); %#ok<NASGU>

        fig = figure('Visible', local_onoff(P.Visible), 'Color', 'w', 'Position', [100 100 650 650]);
        ax = axes(fig);
        axis(ax, 'equal');
        hold(ax, 'on');
        box(ax, 'on');

        lims = local_axis_limits(tracks, frameIdx, P.ArenaPadding);
        xlim(ax, lims(1:2));
        ylim(ax, lims(3:4));
        if P.FlipY
            set(ax, 'YDir', 'reverse');
        end
        xlabel(ax, 'x'); ylabel(ax, 'y');

        motifId = local_get_scalar(row, {'motif_id','motif','cluster_id'}, NaN);
        rankId = local_get_scalar(row, {'rank','exemplar_rank','neighbor_rank'}, NaN);
        sess = local_get_scalar(row, {'session_index','session'}, NaN);

        for f = frameIdx
            cla(ax);
            title(ax, sprintf('Motif %g | exemplar %g | session %g | frame %d', motifId, rankId, sess, f), ...
                'Interpreter', 'none');
            local_draw_dyad_skeleton(ax, squeeze(tracks(f,:,:,:)), edges, P);
            xlim(ax, lims(1:2)); ylim(ax, lims(3:4));
            if P.FlipY
                set(ax, 'YDir', 'reverse');
            end
            drawnow;
            writeVideo(vw, getframe(fig));
        end
        close(vw);
        if isvalid(fig); close(fig); end
        actualProfile = profile;
        return
    catch ME
        lastErr = ME;
        try
            if exist('fig','var') && isvalid(fig); close(fig); end
        catch
        end
        continue
    end
end
if ~isempty(lastErr)
    rethrow(lastErr);
else
    error('No VideoWriter profiles were available.');
end
end

function local_draw_dyad_skeleton(ax, frameTracks, edges, P)
% frameTracks: nodes x xy x mice
if ndims(frameTracks) == 2
    frameTracks = reshape(frameTracks, size(frameTracks,1), size(frameTracks,2), 1);
end
nMice = size(frameTracks,3);
mouseColors = [0.1 0.35 0.9; 0.9 0.2 0.15; 0.25 0.65 0.25; 0.5 0.2 0.75];
for m = 1:nMice
    xy = frameTracks(:,:,m);
    if size(xy,2) < 2
        continue
    end
    c = mouseColors(1+mod(m-1,size(mouseColors,1)),:);
    for e = 1:size(edges,1)
        a = edges(e,1); b = edges(e,2);
        if a <= size(xy,1) && b <= size(xy,1) && all(isfinite(xy([a b],1))) && all(isfinite(xy([a b],2)))
            plot(ax, xy([a b],1), xy([a b],2), '-', 'Color', c, 'LineWidth', P.LineWidth);
        end
    end
    good = isfinite(xy(:,1)) & isfinite(xy(:,2));
    scatter(ax, xy(good,1), xy(good,2), P.MarkerSize, c, 'filled', 'MarkerFaceAlpha', 0.85);

    % Nose/body emphasis if nodes exist.
    if size(xy,1) >= 1 && all(isfinite(xy(1,1:2)))
        scatter(ax, xy(1,1), xy(1,2), P.MarkerSize*1.4, c, 'filled', 'MarkerEdgeColor', 'k');
    end
    if size(xy,1) >= 10 && all(isfinite(xy(10,1:2)))
        scatter(ax, xy(10,1), xy(10,2), P.MarkerSize*1.2, c, 'filled', 'MarkerEdgeColor', 'w');
    end
end
end

function writerInfo = local_choose_writer(requestedProfile, requestedFormat)
requestedProfile = char(requestedProfile);
requestedFormat = lower(char(requestedFormat));

if strcmpi(requestedProfile, 'auto') || isempty(requestedProfile)
    switch requestedFormat
        case 'mp4'
            candidates = {'MPEG-4','Motion JPEG AVI','Uncompressed AVI'};
        case 'avi'
            candidates = {'Motion JPEG AVI','Uncompressed AVI'};
        otherwise
            candidates = {'MPEG-4','Motion JPEG AVI','Uncompressed AVI'};
    end
else
    candidates = [{requestedProfile}, {'Motion JPEG AVI','Uncompressed AVI'}];
end

% Do not call VideoWriter just to test here, because some MATLAB installs throw
% only when opening. local_write_skeleton_clip will try in this order.
candidates = unique(candidates, 'stable');
writerInfo = struct();
writerInfo.profile = candidates{1};
writerInfo.profileList = candidates;
writerInfo.extension = local_profile_extension(candidates{1});
end

function ext = local_profile_extension(profile)
if strcmpi(profile, 'MPEG-4')
    ext = '.mp4';
else
    ext = '.avi';
end
end

function local_safe_close(vw)
try
    close(vw);
catch
end
end

function s = local_onoff(tf)
if logical(tf)
    s = 'on';
else
    s = 'off';
end
end

function edges = local_default_mouse_edges(nParts)
%LOCAL_DEFAULT_MOUSE_EDGES Mouse-like skeleton graph.
%
% Default assumes the common 13-node SLEAP ordering used elsewhere in this
% project, where node 1 is the nose/head, node 2 is neck/head-center, nodes
% 9/10/13 are body/tail-base/mid-body landmarks, and nodes 5:8 are paws.
% The core body chain is:
%   head/nose -> neck -> center/midbody -> tail_base -> tail/midtail -> tail_tip
% plus left/right forepaw and hindpaw edges to make the skeleton read as a mouse.
%
% If fewer keypoints are supplied, invalid edges are dropped automatically.
if nargin < 1 || isempty(nParts) || ~isfinite(nParts)
    nParts = 13;
end

edges13 = [
    1 2   % nose/head -> neck
    2 13  % neck -> body center/mid-body
    13 9  % body center -> lower body/hip
    9 10  % lower body/hip -> tail base
    10 11 % tail base -> tail/midtail
    11 12 % tail/midtail -> tail tip
    13 5  % body center -> left/front paw region
    13 6  % body center -> right/front paw region
    9 7   % hip -> left hind paw
    9 8   % hip -> right hind paw
    5 6   % front paw crossbar
    7 8   % hind paw crossbar
    2 3   % neck/head -> ear/side node
    2 4]; % neck/head -> ear/side node

edges = edges13(all(edges13 <= nParts, 2), :);
if isempty(edges) && nParts >= 2
    edges = [(1:nParts-1)' (2:nParts)'];
end
end

function lims = local_axis_limits(tracks, frameIdx, pad)
sub = tracks(frameIdx,:,:,:);
x = sub(:,:,1,:);
y = sub(:,:,2,:);
x = x(isfinite(x)); y = y(isfinite(y));
if isempty(x) || isempty(y)
    lims = [0 1024 0 1024];
    return
end
lims = [min(x)-pad max(x)+pad min(y)-pad max(y)+pad];
if lims(1) == lims(2); lims(1) = lims(1)-1; lims(2) = lims(2)+1; end
if lims(3) == lims(4); lims(3) = lims(3)-1; lims(4) = lims(4)+1; end
end


function nParts = local_guess_num_bodyparts(PoseData)
nParts = 13;
try
    tracksBySession = local_extract_tracks_by_session(PoseData);
    for ii = 1:numel(tracksBySession)
        if ~isempty(tracksBySession{ii})
            nParts = size(tracksBySession{ii}, 2);
            return
        end
    end
catch
    nParts = 13;
end
end

function tracksBySession = local_extract_tracks_by_session(PoseData)
tracksBySession = {};

if isnumeric(PoseData)
    tracksBySession = {local_standardize_tracks(PoseData)};
    return
end

if iscell(PoseData)
    tracksBySession = cell(size(PoseData));
    for i = 1:numel(PoseData)
        tracksBySession{i} = local_standardize_tracks(PoseData{i});
    end
    return
end

if ~isstruct(PoseData)
    return
end

% DyadSet-like: session array/cell
if isfield(PoseData, 'session')
    S = PoseData.session;
    if iscell(S)
        tracksBySession = cell(numel(S),1);
        for i = 1:numel(S)
            tracksBySession{i} = local_extract_one_session_tracks(S{i});
        end
    elseif isstruct(S)
        tracksBySession = cell(numel(S),1);
        for i = 1:numel(S)
            tracksBySession{i} = local_extract_one_session_tracks(S(i));
        end
    end
    return
end

if isfield(PoseData, 'dyadCell')
    C = PoseData.dyadCell;
    tracksBySession = cell(numel(C),1);
    for i = 1:numel(C)
        tracksBySession{i} = local_extract_one_session_tracks(C{i});
    end
    return
end

% dbase-like struct array
if numel(PoseData) > 1
    tracksBySession = cell(numel(PoseData),1);
    for i = 1:numel(PoseData)
        tracksBySession{i} = local_extract_one_session_tracks(PoseData(i));
    end
    return
end

tracksBySession = {local_extract_one_session_tracks(PoseData)};
end

function tracks = local_extract_one_session_tracks(S)
tracks = [];
if isnumeric(S)
    tracks = local_standardize_tracks(S);
    return
end
if ~isstruct(S)
    return
end
fields = {'tracks','SLEAPtracks','pose','points','keypoints','xy'};
for i = 1:numel(fields)
    f = fields{i};
    if isfield(S, f) && isnumeric(S.(f))
        tracks = local_standardize_tracks(S.(f));
        return
    end
end
if isfield(S, 'raw')
    tracks = local_extract_one_session_tracks(S.raw);
end
end

function T = local_standardize_tracks(T)
% Return frames x nodes x xy x mice.
if isempty(T)
    return
end
sz = size(T);
if ndims(T) == 4
    % Expected frames x nodes x xy x mice.
    if sz(3) == 2
        return
    end
    % Sometimes frames x nodes x mice x xy.
    if sz(4) == 2
        T = permute(T, [1 2 4 3]);
        return
    end
end
if ndims(T) == 3
    % frames x nodes x xy -> one mouse
    if sz(3) == 2
        T = reshape(T, sz(1), sz(2), 2, 1);
        return
    end
end
error('render_motif_skeleton_exemplars:BadTrackShape', ...
    'Track array must be [frames x nodes x 2 x mice] or [frames x nodes x 2].');
end

function val = local_get_scalar(row, names, defaultVal)
val = defaultVal;
if ~istable(row)
    return
end
for i = 1:numel(names)
    nm = names{i};
    if ismember(nm, row.Properties.VariableNames)
        x = row.(nm);
        if iscell(x); x = x{1}; end
        if isstring(x) || ischar(x); x = str2double(x); end
        if isnumeric(x) && ~isempty(x)
            val = double(x(1));
            return
        end
    end
end
end
