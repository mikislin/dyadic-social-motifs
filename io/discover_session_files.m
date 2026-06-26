function files = discover_session_files(sessionDir)
%DISCOVER_SESSION_FILES Return sorted session_*.mat files with raw indices.

assert(isfolder(sessionDir), 'discover_session_files:MissingDir', ...
    'Session directory does not exist: %s', sessionDir);

D = dir(fullfile(sessionDir, 'session_*.mat'));
assert(~isempty(D), 'discover_session_files:NoFiles', ...
    'No session_*.mat files found in %s.', sessionDir);

names = string({D.name})';
rawIndex = nan(numel(D),1);
for i = 1:numel(D)
    tok = regexp(D(i).name, '^session_(\d+)\.mat$', 'tokens', 'once');
    if ~isempty(tok)
        rawIndex(i) = str2double(tok{1});
    end
end

[~, order] = sort(rawIndex);
D = D(order);
rawIndex = rawIndex(order);

files = table();
files.raw_index = rawIndex;
files.session_file = string({D.name})';
files.session_path = string({D.name})';
files.session_rel_path = files.session_path;
files.session_source_dir = repmat("files2run", numel(D), 1);
files.file_bytes = [D.bytes]';
files.file_modified = datetime([D.datenum]', 'ConvertFrom', 'datenum');
end
