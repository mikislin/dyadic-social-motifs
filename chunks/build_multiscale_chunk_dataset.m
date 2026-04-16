function ChunkSet = build_multiscale_chunk_dataset(dyadCell, varargin)
%BUILD_MULTISCALE_CHUNK_DATASET Build dense multi-scale chunks from dyad structs.
%
% This version supports either an explicit scale list via 'chunkSec' or an
% automatically generated dense 25-point log-spaced scale bank via
% 'useLogScaleBank', true.
%
% Required dyad fields
%   dyad.X, dyad.featureNames, dyad.featureMeta, dyad.frameMask
%   and dyad.time_s or dyad.time
%
% Output
%   ChunkSet.scale(s).X      [N x L x D] scaled transformed chunk tensor
%   ChunkSet.scale(s).Xraw   [N x L x D] transformed but unscaled chunk tensor
%   ChunkSet.scale(s).valid  [N x L]     frame-valid mask within chunk
%   ChunkSet.scale(s).meta   table       chunk metadata
%   ChunkSet.chunkTable      combined metadata across scales

p = inputParser;
p.addParameter('chunkSec', [0.4 0.8 1.5 5.0], @(x)isnumeric(x) && isvector(x) && all(x > 0));
p.addParameter('useLogScaleBank', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('minScaleSec', 0.2, @(x)isscalar(x) && x > 0);
p.addParameter('maxScaleSec', 8.0, @(x)isscalar(x) && x > 0);
p.addParameter('nScales', 25, @(x)isscalar(x) && x >= 3);
p.addParameter('strideSec', 0.1, @(x)isscalar(x) && x > 0);
p.addParameter('anchorMode', 'center', @(x)ischar(x) || isstring(x));
p.addParameter('minValidFrac', 0.90, @(x)isscalar(x) && x >= 0 && x <= 1);
p.addParameter('requireAllAnchorsValid', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('anchorMasks', {}, @(x)iscell(x) || islogical(x) || isnumeric(x));
p.addParameter('sessionIds', [], @(x)isnumeric(x) || iscell(x));
p.parse(varargin{:});
P = p.Results;

if P.useLogScaleBank
    chunkSec = make_logscale_chunk_bank('minSec', P.minScaleSec, 'maxSec', P.maxScaleSec, 'nScales', P.nScales);
else
    chunkSec = unique(round(sort(P.chunkSec(:)'), 4));
end

anchorMode = lower(string(P.anchorMode));
assert(any(anchorMode == ["center","past"]), 'anchorMode must be center or past.');

if ~iscell(dyadCell)
    dyadCell = num2cell(dyadCell);
end
nSess = numel(dyadCell);
anchorMasks = i_normalize_anchor_masks(dyadCell, P.anchorMasks);
sessionIds = i_normalize_session_ids(nSess, P.sessionIds);

% Prepare sequences with shared scaling model.
seqsRaw = cell(nSess,1);
Xall = [];
for i = 1:nSess
    seqsRaw{i} = prepare_dyad_timeseries(dyadCell{i}, 'anchorMask', anchorMasks{i});
    Xall = [Xall; seqsRaw{i}.X(seqsRaw{i}.validMask,:)]; %#ok<AGROW>
end
channelMeta = seqsRaw{1}.channelMeta;
stats = i_fit_global_stats(Xall, channelMeta);

seqs = cell(nSess,1);
for i = 1:nSess
    seqs{i} = prepare_dyad_timeseries(dyadCell{i}, 'anchorMask', anchorMasks{i}, 'stats', stats);
    seqs{i}.sessionId = sessionIds{i};
end

obsNames = seqs{1}.obsNames;
D = size(seqs{1}.Xscaled, 2);
Scale = repmat(struct('chunkSec',[],'nFrames',[],'X',[],'Xraw',[],'valid',[],'meta',table(),'sessionIndex',[]), numel(chunkSec), 1);
allMeta = table();
chunkIdOffset = 0;

for s = 1:numel(chunkSec)
    thisSec = chunkSec(s);
    L = max(1, round(thisSec * seqs{1}.fps));
    Xscale = zeros(0, L, D, 'single');
    XrawScale = zeros(0, L, D, 'single');
    Vscale = false(0, L);
    metaRows = table();
    sessIndex = zeros(0,1);

    for i = 1:nSess
        [Xc, Xr, Vc, meta] = extract_multiscale_chunks_from_seq(seqs{i}, ...
            'chunkSec', thisSec, ...
            'strideSec', P.strideSec, ...
            'anchorMode', anchorMode, ...
            'minValidFrac', P.minValidFrac, ...
            'requireAllAnchorsValid', P.requireAllAnchorsValid);
        if isempty(meta)
            continue
        end
        nNew = size(Xc,1);
        meta.session_index = repmat(i, nNew, 1);
        meta.session_id = repmat(string(seqs{i}.sessionId), nNew, 1);
        meta.scale_index = repmat(s, nNew, 1);
        meta.chunk_id = ((1:nNew)' + chunkIdOffset);
        chunkIdOffset = chunkIdOffset + nNew;

        Xscale(end+1:end+nNew,:,:) = single(Xc); %#ok<AGROW>
        XrawScale(end+1:end+nNew,:,:) = single(Xr); %#ok<AGROW>
        Vscale(end+1:end+nNew,:) = Vc; %#ok<AGROW>
        metaRows = [metaRows; meta]; %#ok<AGROW>
        sessIndex = [sessIndex; repmat(i, nNew, 1)]; %#ok<AGROW>
    end

    Scale(s).chunkSec = thisSec;
    Scale(s).nFrames = L;
    Scale(s).X = Xscale;
    Scale(s).Xraw = XrawScale;
    Scale(s).valid = Vscale;
    Scale(s).meta = metaRows;
    Scale(s).sessionIndex = sessIndex;
    if ~isempty(metaRows)
        allMeta = [allMeta; metaRows]; %#ok<AGROW>
    end
end

ChunkSet = struct();
ChunkSet.stats = stats;
ChunkSet.obsNames = obsNames;
ChunkSet.channelMeta = channelMeta;
ChunkSet.featureNames = seqs{1}.featureNames;
ChunkSet.featureMeta = seqs{1}.featureMeta;
ChunkSet.chunkSec = chunkSec;
ChunkSet.strideSec = P.strideSec;
ChunkSet.anchorMode = char(anchorMode);
ChunkSet.scale = Scale;
ChunkSet.chunkTable = allMeta;
ChunkSet.sessions = seqs;
ChunkSet.nSessions = nSess;
ChunkSet.nObs = D;
ChunkSet.scaleBankMeta = table((1:numel(chunkSec))', chunkSec(:), cellfun(@(x) i_scale_band_label(x), num2cell(chunkSec(:))), ...
    'VariableNames', {'scale_index','chunk_sec','temporal_band'});
end

function [Xchunks, XrawChunks, validChunks, meta] = extract_multiscale_chunks_from_seq(Seq, varargin)
p = inputParser;
p.addParameter('chunkSec', 0.8, @(x)isscalar(x) && x > 0);
p.addParameter('strideSec', 0.1, @(x)isscalar(x) && x > 0);
p.addParameter('anchorMode', "center", @(x)ischar(x) || isstring(x));
p.addParameter('minValidFrac', 0.90, @(x)isscalar(x) && x >= 0 && x <= 1);
p.addParameter('requireAllAnchorsValid', true, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});
P = p.Results;

anchorMode = lower(string(P.anchorMode));
X = Seq.Xscaled;
Xraw = Seq.X;
validMask = logical(Seq.validMask(:));
time = Seq.time(:);
T = size(X,1);
D = size(X,2);
fps = Seq.fps;
L = max(1, round(P.chunkSec * fps));
strideFrames = max(1, round(P.strideSec * fps));
[leftFrames, rightFrames, anchorOffset] = i_chunk_geometry(L, anchorMode);
anchorCandidates = (anchorOffset + 1):strideFrames:(T - rightFrames);
anchorCandidates = anchorCandidates(:);

if isempty(anchorCandidates)
    Xchunks = zeros(0, L, D, 'single'); XrawChunks = zeros(0, L, D, 'single'); validChunks = false(0, L); meta = table(); return
end
if P.requireAllAnchorsValid
    anchorCandidates = anchorCandidates(validMask(anchorCandidates));
end
nCand = numel(anchorCandidates);
Xchunks = nan(nCand, L, D, 'single');
XrawChunks = nan(nCand, L, D, 'single');
validChunks = false(nCand, L);
startFrame = zeros(nCand,1); stopFrame = zeros(nCand,1); anchorFrame = zeros(nCand,1);
validFrac = zeros(nCand,1); startTime = zeros(nCand,1); stopTime = zeros(nCand,1); anchorTime = zeros(nCand,1);
keep = false(nCand,1);
for r = 1:nCand
    a = anchorCandidates(r);
    idx = (a - leftFrames):(a + rightFrames);
    if idx(1) < 1 || idx(end) > T, continue, end
    v = validMask(idx);
    vf = mean(v);
    if vf < P.minValidFrac, continue, end
    keep(r) = true;
    startFrame(r) = idx(1); stopFrame(r) = idx(end); anchorFrame(r) = a;
    validFrac(r) = vf; startTime(r) = time(idx(1)); stopTime(r) = time(idx(end)); anchorTime(r) = time(a);
    validChunks(r,:) = v(:)';
    Xchunks(r,:,:) = single(X(idx,:));
    XrawChunks(r,:,:) = single(Xraw(idx,:));
end
Xchunks = Xchunks(keep,:,:); XrawChunks = XrawChunks(keep,:,:); validChunks = validChunks(keep,:);
meta = table();
if any(keep)
    meta.anchor_frame = anchorFrame(keep);
    meta.start_frame = startFrame(keep);
    meta.stop_frame = stopFrame(keep);
    meta.anchor_time_s = anchorTime(keep);
    meta.start_time_s = startTime(keep);
    meta.stop_time_s = stopTime(keep);
    meta.chunk_sec = repmat(P.chunkSec, nnz(keep), 1);
    meta.chunk_frames = repmat(L, nnz(keep), 1);
    meta.stride_sec = repmat(P.strideSec, nnz(keep), 1);
    meta.valid_frac = validFrac(keep);
    meta.anchor_mode = repmat(string(anchorMode), nnz(keep), 1);
    meta.temporal_band = repmat(i_scale_band_label(P.chunkSec), nnz(keep), 1);
end
end

function masks = i_normalize_anchor_masks(dyadCell, anchorMasks)
nSess = numel(dyadCell);
if isempty(anchorMasks)
    masks = cell(nSess,1);
    for i = 1:nSess
        if isfield(dyadCell{i}, 'frameMask') && ~isempty(dyadCell{i}.frameMask)
            masks{i} = logical(dyadCell{i}.frameMask(:));
        else
            T = i_num_frames_from_dyad(dyadCell{i});
            masks{i} = true(T,1);
        end
    end
elseif ~iscell(anchorMasks)
    masks = cell(nSess,1);
    for i = 1:nSess
        masks{i} = logical(anchorMasks(:));
    end
else
    assert(numel(anchorMasks) == nSess, 'anchorMasks length mismatch.');
    masks = anchorMasks;
    for i = 1:nSess, masks{i} = logical(masks{i}(:)); end
end
end

function sessionIds = i_normalize_session_ids(nSess, sessionIds)
if isempty(sessionIds)
    sessionIds = num2cell(1:nSess);
elseif ~iscell(sessionIds)
    assert(numel(sessionIds) == nSess, 'sessionIds length mismatch.');
    sessionIds = num2cell(sessionIds);
else
    assert(numel(sessionIds) == nSess, 'sessionIds length mismatch.');
end
end

function stats = i_fit_global_stats(X, channelMeta)
D = size(X,2); stats = struct();
stats.center = zeros(1,D); stats.scale = ones(1,D); stats.impute = zeros(1,D);
for d = 1:D
    x = X(:,d); ok = isfinite(x);
    if ~any(ok), continue, end
    if channelMeta.ChannelType(d) == "boolean"
        stats.center(d) = 0; stats.scale(d) = 1; stats.impute(d) = median(x(ok));
    else
        stats.center(d) = median(x(ok));
        s = iqr(x(ok)); if ~(isfinite(s) && s > 0), s = std(x(ok),0); end
        if ~(isfinite(s) && s > 0), s = 1; end
        stats.scale(d) = s; stats.impute(d) = stats.center(d);
    end
end
end

function [leftFrames, rightFrames, anchorOffset] = i_chunk_geometry(L, anchorMode)
switch char(anchorMode)
    case 'center'
        leftFrames = floor((L - 1) / 2); rightFrames = ceil((L - 1) / 2); anchorOffset = leftFrames;
    case 'past'
        leftFrames = L - 1; rightFrames = 0; anchorOffset = leftFrames;
    otherwise
        error('Unsupported anchorMode');
end
end

function T = i_num_frames_from_dyad(dyad)
if isfield(dyad,'X') && ~isempty(dyad.X)
    T = size(dyad.X,1);
elseif isfield(dyad,'time_s')
    T = numel(dyad.time_s);
elseif isfield(dyad,'time')
    T = numel(dyad.time);
else
    error('Could not infer number of frames from dyad.');
end
end

function label = i_scale_band_label(sec)
if sec < 0.8
    label = "micro";
elseif sec < 2.5
    label = "motif";
else
    label = "context";
end
end
