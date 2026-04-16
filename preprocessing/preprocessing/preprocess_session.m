function sessionPreproc = preprocess_session(sessionRaw, P)
%PREPROCESS_SESSION Preprocess one dyadic session.

if nargin < 2 || isempty(P)
    P = default_preprocessing_params();
end
P = validate_preprocessing_params(P);
sessionRaw = validate_session_inputs(sessionRaw, P);

tracksRaw = sessionRaw.SLEAPtracks;
time = sessionRaw.time;
T = size(tracksRaw,1);
nAnimals = size(tracksRaw,4);

if isfield(sessionRaw, 'excludedFrames') && ~isempty(sessionRaw.excludedFrames)
    excludedFrames = logical(sessionRaw.excludedFrames(:));
else
    excludedFrames = false(T,1);
end

frameHasAnyNaN = squeeze(any(any(any(isnan(tracksRaw),2),3),4));
frameHasAnyNaN = logical(frameHasAnyNaN(:));
frameAllMissingPerAnimal = false(T, nAnimals);
for m = 1:nAnimals
    tr = tracksRaw(:,:,:,m);
    frameAllMissingPerAnimal(:,m) = squeeze(all(all(isnan(tr),2),3));
end

qcTemplate = struct( ...
    'confidenceThresholds', [], ...
    'confidenceInfo', struct(), ...
    'lowConfMask', [], ...
    'jumpMask', [], ...
    'interpMask', [], ...
    'geomMask', [], ...
    'arenaMask', [], ...
    'finalNanMask', [], ...
    'jumpMeta', [], ...
    'smoothMeta', [], ...
    'geometryInfo', struct(), ...
    'arena', struct(), ...
    'nodeSummary', [] ...
);
qcAnimals = repmat(qcTemplate, nAnimals, 1);

tracksClean = nan(size(tracksRaw));

for m = 1:nAnimals
    tracksRawAnimal = reshape(tracksRaw(:,:,:,m), size(tracksRaw,1), size(tracksRaw,2), size(tracksRaw,3));
    if isempty(sessionRaw.SLEAPscores)
        scoresAnimal = [];
    else
        scoresAnimal = sessionRaw.SLEAPscores(:,:,min(m, size(sessionRaw.SLEAPscores,3)));
    end
    [tracksClean(:,:,:,m), qcOut] = preprocess_animal_tracks(tracksRawAnimal, time, excludedFrames, scoresAnimal, P);
    qcAnimals(m) = qcOut;
    if P.debug.enabled && P.debug.verbose
        fprintf('Animal %d: low-conf frac = %.4f, jump frac = %.4f, final NaN frac = %.4f\n', ...
            m, full(mean(qcOut.lowConfMask(:))), full(mean(qcOut.jumpMask(:))), full(mean(qcOut.finalNanMask(:))));
    end
end

[badframes, qcFrame] = make_badframes(tracksClean, excludedFrames, qcAnimals, P);

sessionPreproc = struct();
sessionPreproc.session_id = sessionRaw.session_id;
sessionPreproc.params = P;
sessionPreproc.raw = struct();
if P.output.return_raw
    sessionPreproc.raw.SLEAPtracks = tracksRaw;
else
    sessionPreproc.raw.SLEAPtracks = [];
end
sessionPreproc.raw.time = time;
sessionPreproc.raw.excludedFrames = excludedFrames;
sessionPreproc.clean = struct();
sessionPreproc.clean.tracks = tracksClean;
sessionPreproc.qc = struct();
sessionPreproc.qc.animals = qcAnimals;
sessionPreproc.qc.frames = qcFrame;
sessionPreproc.qc.badframes = badframes;
sessionPreproc.qc.sessionStats = summarize_preprocessing_qc(sessionPreproc);

sessionPreproc.debug = struct();
sessionPreproc.debug.excludedFramesInput = excludedFrames;
sessionPreproc.debug.frameHasAnyNaN = frameHasAnyNaN;
sessionPreproc.debug.frameAllMissingPerAnimal = frameAllMissingPerAnimal;
sessionPreproc.debug.paramsVersion = P.meta.version;

assert(isequal(size(sessionPreproc.clean.tracks), size(sessionRaw.SLEAPtracks)), 'Output tracks size mismatch');
assert(isequal(size(sessionPreproc.qc.badframes), [T nAnimals]), 'badframes size mismatch');
end
