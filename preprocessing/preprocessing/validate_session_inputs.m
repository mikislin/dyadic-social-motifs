function sessionRaw = validate_session_inputs(sessionRaw, P)
%VALIDATE_SESSION_INPUTS Validate a raw session struct.

assert(isstruct(sessionRaw), 'sessionRaw must be a struct');
assert(isfield(sessionRaw, 'SLEAPtracks'), 'sessionRaw.SLEAPtracks is required');

tracks = double(sessionRaw.SLEAPtracks);
if ndims(tracks) == 3
    tracks = reshape(tracks, size(tracks,1), size(tracks,2), size(tracks,3), 1);
elseif ndims(tracks) ~= 4
    error('validate_session_inputs:BadTracksShape', ...
        'SLEAPtracks must be [T x nodes x coords] or [T x nodes x coords x animals]');
end

[T, nNodes, nCoords, nAnimals] = size(tracks);
assert(T >= 2, 'Need at least two frames');
assert(nNodes >= 1, 'Need at least one node');
assert(nCoords == P.data.n_coords_expected, 'Unexpected number of coordinates');
assert(ismember(nAnimals, [1 2]), 'This package currently supports 1 or 2 animals');

sessionRaw.SLEAPtracks = tracks;
sessionRaw.nAnimals = nAnimals;

if ~isfield(sessionRaw, 'time') || isempty(sessionRaw.time)
    sessionRaw.time = (0:T-1)' ./ P.data.fps;
else
    assert(numel(sessionRaw.time) == T, 'time length must match frame count');
    sessionRaw.time = double(sessionRaw.time(:));
    assert(all(diff(sessionRaw.time) > 0), 'time must be strictly increasing');
end

if ~isfield(sessionRaw, 'excludedFrames') || isempty(sessionRaw.excludedFrames)
    sessionRaw.excludedFrames = false(T,1);
else
    assert(numel(sessionRaw.excludedFrames) == T, 'excludedFrames length must match frame count');
    sessionRaw.excludedFrames = logical(sessionRaw.excludedFrames(:));
end

if isfield(sessionRaw, 'SLEAPscores') && ~isempty(sessionRaw.SLEAPscores)
    scores = double(sessionRaw.SLEAPscores);

    if ismatrix(scores)
        assert(size(scores,1) == T, 'SLEAPscores must match frame count');
        assert(size(scores,2) == nNodes, 'SLEAPscores must match node count');
        scores = repmat(scores, 1, 1, nAnimals);
    elseif ndims(scores) == 3
        assert(size(scores,1) == T, 'SLEAPscores must match frame count');
        assert(size(scores,2) == nNodes, 'SLEAPscores must match node count');
        assert(size(scores,3) == nAnimals, 'SLEAPscores third dim must match number of animals');
    elseif ndims(scores) == 4 && size(scores,3) == 1
        assert(size(scores,1) == T, 'SLEAPscores must match frame count');
        assert(size(scores,2) == nNodes, 'SLEAPscores must match node count');
        assert(size(scores,4) == nAnimals, 'SLEAPscores 4th dim must match number of animals');
        scores = squeeze(scores(:,:,1,:));
        if nAnimals == 1
            scores = reshape(scores, T, nNodes, 1);
        end
    else
        error('validate_session_inputs:BadScoresShape', ...
            'SLEAPscores must be [T x nodes], [T x nodes x animals], or [T x nodes x 1 x animals]');
    end

    sessionRaw.SLEAPscores = scores;
elseif P.confidence.require_scores
    error('validate_session_inputs:MissingScores', ...
        'SLEAPscores are required by current parameters');
else
    sessionRaw.SLEAPscores = [];
end

if ~isfield(sessionRaw, 'session_id')
    if isfield(sessionRaw, 'fileID')
        sessionRaw.session_id = sprintf('session_%s', string(sessionRaw.fileID));
    else
        sessionRaw.session_id = 'session_unknown';
    end
end
end
