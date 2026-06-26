function [sessionRaw, manifestRow, animalQC] = load_session_raw(manifestOrPath, idxOrRow, varargin)
%LOAD_SESSION_RAW Load a sessionRaw struct and apply animal-identity QC.
%
% Usage
%   [sessionRaw, row, qc] = load_session_raw(manifestTable, rowIndex)
%   [sessionRaw, row, qc] = load_session_raw(sessionPath)

p = inputParser;
p.addParameter('applyAnimalQC', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('sessionRoot', '', @(x)ischar(x) || isstring(x));
p.parse(varargin{:});
P = p.Results;

manifestRow = table();
if istable(manifestOrPath)
    M = manifestOrPath;
    if nargin < 2 || isempty(idxOrRow)
        idxOrRow = 1;
    end
    manifestRow = M(idxOrRow,:);
    sessionPath = local_resolve_session_path(manifestRow, P.sessionRoot);
else
    sessionPath = string(manifestOrPath);
end

S = load(sessionPath, 'sessionRaw');
assert(isfield(S, 'sessionRaw'), 'load_session_raw:MissingSessionRaw', ...
    'File does not contain sessionRaw: %s', sessionPath);
sessionRaw = S.sessionRaw;

animalQC = qc_session_animals(sessionRaw);
if P.applyAnimalQC && ~isempty(animalQC.keep_indices)
    keep = animalQC.keep_indices;
    if ndims(sessionRaw.SLEAPtracks) == 4
        sessionRaw.SLEAPtracks = sessionRaw.SLEAPtracks(:,:,:,keep);
    end
    if isfield(sessionRaw, 'SLEAPscores') && ~isempty(sessionRaw.SLEAPscores)
        if ndims(sessionRaw.SLEAPscores) == 3
            sessionRaw.SLEAPscores = sessionRaw.SLEAPscores(:,:,keep);
        elseif ndims(sessionRaw.SLEAPscores) == 4
            sessionRaw.SLEAPscores = sessionRaw.SLEAPscores(:,:,:,keep);
        end
    end
    sessionRaw.animal_qc = animalQC;
end
end

function sessionPath = local_resolve_session_path(row, sessionRoot)
vars = string(row.Properties.VariableNames);
sessionRoot = string(sessionRoot);

if ismember("session_path", vars) && strlength(string(row.session_path)) > 0
    candidate = string(row.session_path);
elseif ismember("session_rel_path", vars) && strlength(string(row.session_rel_path)) > 0
    candidate = string(row.session_rel_path);
elseif ismember("session_file", vars)
    candidate = string(row.session_file);
else
    error('load_session_raw:MissingSessionPath', ...
        'Manifest row must contain session_path, session_rel_path, or session_file.');
end

if local_is_absolute_path(candidate)
    sessionPath = candidate;
    return
end

if sessionRoot == ""
    paths = paper_paths();
    if isfield(paths, 'files2runDir') && strlength(string(paths.files2runDir)) > 0
        sessionRoot = string(paths.files2runDir);
    end
end

assert(sessionRoot ~= "", 'load_session_raw:MissingSessionRoot', ...
    'Manifest session_path is relative (%s). Provide sessionRoot or configure raw data paths.', ...
    candidate);
sessionPath = string(fullfile(sessionRoot, candidate));
end

function tf = local_is_absolute_path(pathName)
pathName = string(pathName);
tf = startsWith(pathName, filesep) || startsWith(pathName, "/") || ...
    ~cellfun('isempty', regexp(cellstr(pathName), '^[A-Za-z]:[\\/]', 'once'));
tf = logical(tf(1));
end
