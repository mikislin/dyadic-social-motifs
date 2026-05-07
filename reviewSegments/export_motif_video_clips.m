function ClipTable = export_motif_video_clips(Exemplars, varargin)
%EXPORT_MOTIF_VIDEO_CLIPS Export video clips for exemplar bouts using VideoReader/VideoWriter.
%
% ClipTable = export_motif_video_clips(Exemplars, 'OutputDir', 'motif_clips')
%
% Requires Exemplars.table.video_path and clip start/stop times. If video
% files are unavailable, use find_motif_video_exemplars only and inspect its
% table.

p = inputParser;
p.addParameter('OutputDir', fullfile(pwd, 'motif_clips'), @(x)ischar(x) || isstring(x));
p.addParameter('FrameRate', [], @(x)isempty(x) || (isscalar(x) && x > 0));
p.addParameter('Verbose', true, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});
P = p.Results;

T = Exemplars.table;
if isempty(T)
    ClipTable = T;
    return
end
if ~ismember('video_path', T.Properties.VariableNames)
    error('export_motif_video_clips:MissingVideoPath', 'Exemplars.table lacks video_path.');
end
outDir = char(P.OutputDir);
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

output_path = strings(height(T),1);
export_ok = false(height(T),1);
error_message = strings(height(T),1);

for i = 1:height(T)
    try
        vp = char(T.video_path(i));
        if isempty(vp) || ~exist(vp, 'file')
            error('Video path not found: %s', vp);
        end
        vr = VideoReader(vp);
        if isempty(P.FrameRate)
            outFps = vr.FrameRate;
        else
            outFps = P.FrameRate;
        end
        startT = max(0, T.clip_start_time_s(i));
        stopT = min(vr.Duration, T.clip_stop_time_s(i));
        label = char(T.clip_label(i));
        outPath = fullfile(outDir, [label '.mp4']);
        vw = VideoWriter(outPath, 'MPEG-4');
        vw.FrameRate = outFps;
        open(vw);
        vr.CurrentTime = startT;
        while hasFrame(vr) && vr.CurrentTime <= stopT
            frame = readFrame(vr);
            writeVideo(vw, frame);
        end
        close(vw);
        output_path(i) = string(outPath);
        export_ok(i) = true;
    catch ME
        error_message(i) = string(ME.message);
        try close(vw); catch, end %#ok<CTCH>
    end
end

ClipTable = T;
ClipTable.output_path = output_path;
ClipTable.export_ok = export_ok;
ClipTable.error_message = error_message;

if logical(P.Verbose)
    fprintf('export_motif_video_clips | exported %d / %d clips to %s\n', sum(export_ok), height(T), outDir);
end
end
