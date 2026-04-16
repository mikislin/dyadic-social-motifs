function export_dbase_sessions_to_mat(dbase, outDir)
%EXPORT_DBASE_SESSIONS_TO_MAT Export one MAT file per session for cluster preprocessing.

if ~isfolder(outDir)
    mkdir(outDir);
end

for i = 1:numel(dbase)
    sessionRaw = struct();
    sessionRaw.SLEAPtracks = dbase(i).SLEAPtracks;
    if isfield(dbase, 'SLEAPscores') && ~isempty(dbase(i).SLEAPscores)
        sessionRaw.SLEAPscores = dbase(i).SLEAPscores;
    else
        sessionRaw.SLEAPscores = [];
    end
    if isfield(dbase, 'time') && ~isempty(dbase(i).time)
        sessionRaw.time = dbase(i).time;
    else
        sessionRaw.time = [];
    end
    if isfield(dbase, 'excludedFrames') && ~isempty(dbase(i).excludedFrames)
        sessionRaw.excludedFrames = dbase(i).excludedFrames;
    else
        sessionRaw.excludedFrames = [];
    end
    if isfield(dbase, 'fileID')
        sessionRaw.session_id = sprintf('session_%s', string(dbase(i).fileID));
    else
        sessionRaw.session_id = sprintf('session_%04d', i);
    end
    save(fullfile(outDir, sprintf('session_%04d.mat', i)), 'sessionRaw', '-v7.3');
end
end
