function Data = build_anchor_multiscale_matrix(ChunkSet, EmbedModel, selectedScales, varargin)
%BUILD_ANCHOR_MULTISCALE_MATRIX Align selected chunk embeddings on anchors.
%
% Key features:
%   - tolerant scale matching
%   - exact or nearest anchor alignment across scales
%   - reference-scale anchoring to avoid sparse union artifacts
%   - explicit coverage audit before and after retention
%   - mean imputation for missing scale blocks
%
% Output fields
%   Data.X
%   Data.anchorTable
%   Data.presence
%   Data.coverage
%   Data.selectedScaleSec
%   Data.selectedScaleIdx
%   Data.retainedMask
%   Data.retention
%   Data.source

p = inputParser;
p.addParameter('AnchorSetMode', 'reference', @(x)ischar(x) || isstring(x));
p.addParameter('AnchorAlignment', 'nearest', @(x)ischar(x) || isstring(x));
p.addParameter('ReferenceScaleMode', 'median_selected', @(x)ischar(x) || isstring(x));
p.addParameter('ReferenceScaleSec', [], @(x)isempty(x) || isscalar(x));
p.addParameter('RequireAllSelectedScales', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('MinScalesPresent', [], @(x)isempty(x) || (isscalar(x) && x >= 1));
p.addParameter('MaxMissingScaleFraction', [], @(x)isempty(x) || (isscalar(x) && x >= 0 && x <= 1));
p.addParameter('RetentionRule', 'auto', @(x)ischar(x) || isstring(x));
p.addParameter('AddScalePresenceIndicators', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('ImputeMissingWithScaleMean', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('ScaleMatchToleranceSec', [], @(x)isempty(x) || (isscalar(x) && x > 0));
p.addParameter('ScaleMatchToleranceRel', 0.05, @(x)isscalar(x) && x > 0 && x < 1);
p.addParameter('AnchorToleranceFrames', [], @(x)isempty(x) || (isscalar(x) && x >= 0));
p.addParameter('Verbose', true, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});
P = p.Results;

scaleSecAll = local_get_scale_seconds(ChunkSet);
selectedScaleSecReq = local_parse_selected_scales(selectedScales);
[selectedScaleIdx, selectedScaleSec, scaleMatchInfo] = local_match_selected_scales( ...
    scaleSecAll, selectedScaleSecReq, P.ScaleMatchToleranceSec, P.ScaleMatchToleranceRel);
nSel = numel(selectedScaleIdx);
if nSel == 0
    error('build_anchor_multiscale_matrix:NoSelectedScales', ...
        'No selected scales could be matched to ChunkSet.');
end

scaleData = cell(nSel,1);
for i = 1:nSel
    sIdx = selectedScaleIdx(i);
    meta = local_get_scale_meta(ChunkSet, sIdx);
    Y = local_get_scale_embedding(EmbedModel, ChunkSet, sIdx);
    if size(Y,1) ~= height(meta)
        error('build_anchor_multiscale_matrix:RowMismatch', ...
            'Embedding rows (%d) do not match meta rows (%d) for scale %.4g s.', ...
            size(Y,1), height(meta), scaleSecAll(sIdx));
    end
    scaleData{i} = struct();
    scaleData{i}.scaleIndex = sIdx;
    scaleData{i}.scaleSec = scaleSecAll(sIdx);
    scaleData{i}.requestedScaleSec = scaleMatchInfo.requestedScaleSec(i);
    scaleData{i}.meta = meta;
    scaleData{i}.Y = Y;
    scaleData{i}.meanY = mean(Y, 1, 'omitnan');
    scaleData{i}.dim = size(Y,2);
end

refLocalIdx = local_choose_reference_scale(scaleData, P.ReferenceScaleMode, P.ReferenceScaleSec);
refData = scaleData{refLocalIdx};
refMeta = refData.meta;
refTol = local_choose_anchor_tolerance(scaleData, refLocalIdx, P.AnchorToleranceFrames);

anchorMode = lower(string(P.AnchorSetMode));
alignMode = lower(string(P.AnchorAlignment));
if anchorMode == "reference"
    anchorSession = refMeta.session_index(:);
    anchorFrame = refMeta.anchor_frame(:);
else
    [anchorSession, anchorFrame] = local_build_union_anchor_set(scaleData);
end

nAnchorAll = numel(anchorFrame);
presence = false(nAnchorAll, nSel);
rowMap = cell(nSel,1);
frameDelta = cell(nSel,1);

for i = 1:nSel
    thisMeta = scaleData{i}.meta;
    if anchorMode == "reference"
        [loc, dFrame] = local_match_rows_to_reference(anchorSession, anchorFrame, ...
            thisMeta.session_index(:), thisMeta.anchor_frame(:), alignMode, refTol(i));
    else
        [loc, dFrame] = local_match_rows_to_reference(anchorSession, anchorFrame, ...
            thisMeta.session_index(:), thisMeta.anchor_frame(:), alignMode, refTol(i));
    end
    rowMap{i} = loc;
    frameDelta{i} = dFrame;
    presence(:,i) = loc > 0;
end

nPresentAll = sum(presence, 2);
missingFracAll = 1 - nPresentAll ./ max(nSel,1);
retention = local_choose_retention_rule(nPresentAll, missingFracAll, P, nSel);
retainedMask = retention.retainedMask;

anchorTable = table(anchorSession(retainedMask), anchorFrame(retainedMask), ...
    nPresentAll(retainedMask), missingFracAll(retainedMask), ...
    'VariableNames', {'session_index','anchor_frame','n_scales_present','missing_scale_fraction'});

blockDim = sum(cellfun(@(s)s.dim, scaleData));
X = nan(nnz(retainedMask), blockDim);
colBlock = zeros(nSel,2);
col0 = 1;
for i = 1:nSel
    d = scaleData{i}.dim;
    colBlock(i,:) = [col0, col0 + d - 1];
    col0 = col0 + d;
end

for i = 1:nSel
    rows = rowMap{i}(retainedMask);
    cols = colBlock(i,1):colBlock(i,2);
    ok = rows > 0;
    if any(ok)
        X(ok, cols) = scaleData{i}.Y(rows(ok), :);
    end
    if P.ImputeMissingWithScaleMean
        miss = ~ok;
        if any(miss)
            X(miss, cols) = repmat(scaleData{i}.meanY, nnz(miss), 1);
        end
    end
end

presenceRet = presence(retainedMask,:);
if P.AddScalePresenceIndicators
    X = [X, double(presenceRet)]; %#ok<AGROW>
end

coverage = struct();
coverage.nAnchorsAll = nAnchorAll;
coverage.nAnchorsRetained = nnz(retainedMask);
coverage.nSelectedScales = nSel;
coverage.scaleSecRequested = scaleMatchInfo.requestedScaleSec(:)';
coverage.scaleSecMatched = selectedScaleSec(:)';
coverage.scaleMatchAbsError = scaleMatchInfo.absError(:)';
coverage.scaleSec = selectedScaleSec(:)';
coverage.referenceScaleLocalIdx = refLocalIdx;
coverage.referenceScaleSec = scaleData{refLocalIdx}.scaleSec;
coverage.anchorToleranceFrames = refTol(:)';
coverage.presenceAll = presence;
coverage.presence = presenceRet;
coverage.nPresentAll = nPresentAll;
coverage.missingFracAll = missingFracAll;
coverage.nPresent = anchorTable.n_scales_present;
coverage.missingFrac = anchorTable.missing_scale_fraction;
coverage.perScaleAnchorFracAll = mean(presence, 1);
coverage.perScaleAnchorFrac = mean(presenceRet, 1);
coverage.perAnchorPresentQuantilesAll = local_quantiles_safe(nPresentAll);
coverage.perAnchorMissingQuantilesAll = local_quantiles_safe(missingFracAll);
coverage.perAnchorPresentQuantiles = local_quantiles_safe(anchorTable.n_scales_present);
coverage.perAnchorMissingQuantiles = local_quantiles_safe(anchorTable.missing_scale_fraction);
coverage.frameDeltaMedian = cellfun(@(x) median(abs(x(isfinite(x))), 'omitnan'), frameDelta);
coverage.frameDeltaP90 = cellfun(@(x) local_prctile_safe(abs(x(isfinite(x))), 90), frameDelta);
coverage.colBlock = colBlock;

source = struct();
source.scaleData = scaleData;
source.rowMap = rowMap;
source.frameDelta = frameDelta;

Data = struct();
Data.X = X;
Data.anchorTable = anchorTable;
Data.presence = presenceRet;
Data.coverage = coverage;
Data.selectedScaleSec = selectedScaleSec(:)';
Data.selectedScaleIdx = selectedScaleIdx(:)';
Data.retainedMask = retainedMask;
Data.retention = retention;
Data.source = source;
Data.integrationMode = sprintf('anchor_concat_%s_%s', char(anchorMode), char(alignMode));

if P.Verbose
    fprintf('build_anchor_multiscale_matrix | selected scales = %d | anchors retained = %d / %d | dim = %d\n', ...
        nSel, nnz(retainedMask), nAnchorAll, size(X,2));
    fprintf('  requested scales (s): %s\n', sprintf('%.4g ', coverage.scaleSecRequested));
    fprintf('  matched scales   (s): %s\n', sprintf('%.4g ', coverage.scaleSecMatched));
    fprintf('  scale match |abs error| (s) = %s\n', sprintf('%.4g ', coverage.scaleMatchAbsError));
    fprintf('  reference scale (s): %.4g\n', coverage.referenceScaleSec);
    fprintf('  anchor tolerance (frames): %s\n', sprintf('%.1f ', coverage.anchorToleranceFrames));
    fprintf('  ALL anchors | nScalesPresent quantiles [min p10 p25 p50 p75 p90 max] = %s\n', ...
        sprintf('%.2f ', coverage.perAnchorPresentQuantilesAll));
    fprintf('  ALL anchors | missingFrac quantiles  [min p10 p25 p50 p75 p90 max] = %s\n', ...
        sprintf('%.3f ', coverage.perAnchorMissingQuantilesAll));
    fprintf('  ALL anchors | per-scale coverage = %s\n', sprintf('%.3f ', coverage.perScaleAnchorFracAll));
    fprintf('  ALL anchors | median |frame delta| by scale = %s\n', sprintf('%.2f ', coverage.frameDeltaMedian));
    fprintf('  retained anchors | nScalesPresent quantiles [min p10 p25 p50 p75 p90 max] = %s\n', ...
        sprintf('%.2f ', coverage.perAnchorPresentQuantiles));
    fprintf('  retained anchors | missingFrac quantiles  [min p10 p25 p50 p75 p90 max] = %s\n', ...
        sprintf('%.3f ', coverage.perAnchorMissingQuantiles));
    fprintf('  retained anchors | per-scale coverage = %s\n', sprintf('%.3f ', coverage.perScaleAnchorFrac));
    fprintf('  retention rule: %s\n', retention.description);
end
end

function scaleSec = local_get_scale_seconds(ChunkSet)
assert(isfield(ChunkSet, 'scale'), 'ChunkSet.scale is required.');
n = numel(ChunkSet.scale);
scaleSec = nan(1,n);
for i = 1:n
    sc = ChunkSet.scale(i);
    candidates = {'chunkSec','sec','scaleSec','windowSec','historySec'};
    for j = 1:numel(candidates)
        if isfield(sc, candidates{j}) && ~isempty(sc.(candidates{j}))
            scaleSec(i) = sc.(candidates{j});
            break
        end
    end
end
assert(all(isfinite(scaleSec)), 'Could not determine scale durations from ChunkSet.scale.');
end

function selectedScaleSec = local_parse_selected_scales(selectedScales)
if istable(selectedScales)
    vars = string(selectedScales.Properties.VariableNames);
    candidates = ["chunkSec","scaleSec","selectedScaleSec","sec","windowSec","historySec","scale_s","chunk_sec"];
    hit = intersect(candidates, vars, 'stable');
    if isempty(hit)
        hit = strings(0,1);
        for i = 1:numel(vars)
            vn = vars(i);
            x = selectedScales.(vn);
            if isnumeric(x) && isvector(x)
                nameOk = contains(lower(vn), 'scale') || contains(lower(vn), 'sec') || ...
                    contains(lower(vn), 'window') || contains(lower(vn), 'history') || contains(lower(vn), 'chunk');
                if nameOk
                    hit = vn;
                    break
                end
            end
        end
    end
    if isempty(hit)
        error('build_anchor_multiscale_matrix:BadSelectedScales', ...
            'selectedScales table does not contain a recognizable scale column.');
    end
    selectedScaleSec = selectedScales.(char(hit(1)));
else
    selectedScaleSec = selectedScales;
end
selectedScaleSec = double(selectedScaleSec(:)');
selectedScaleSec = selectedScaleSec(isfinite(selectedScaleSec));
end

function [idx, matchedScaleSec, info] = local_match_selected_scales(allScaleSec, requestedScaleSec, absTolUser, relTol)
allScaleSec = double(allScaleSec(:)');
requestedScaleSec = double(requestedScaleSec(:)');
if isempty(requestedScaleSec)
    idx = zeros(1,0); matchedScaleSec = zeros(1,0);
    info = struct('requestedScaleSec', zeros(1,0), 'matchedScaleSec', zeros(1,0), 'absError', zeros(1,0));
    return
end
if numel(allScaleSec) > 1
    minGap = min(diff(sort(allScaleSec)));
else
    minGap = inf;
end
defaultAbsTol = max([1e-3, 0.02, 0.5 * minGap]);
if isempty(absTolUser), absTol = defaultAbsTol; else, absTol = absTolUser; end
idx = zeros(1, numel(requestedScaleSec));
matchedScaleSec = nan(1, numel(requestedScaleSec));
absErr = nan(1, numel(requestedScaleSec));
for i = 1:numel(requestedScaleSec)
    [d, j] = min(abs(allScaleSec - requestedScaleSec(i)));
    localTol = max(absTol, relTol * max(abs(requestedScaleSec(i)), eps));
    if d > localTol
        error('build_anchor_multiscale_matrix:ScaleMatchFailed', ...
            ['Could not match requested scale %.4g s to ChunkSet. ' ...
             'Nearest available scale is %.4g s (|error| = %.4g s, tolerance = %.4g s).'], ...
            requestedScaleSec(i), allScaleSec(j), d, localTol);
    end
    idx(i) = j; matchedScaleSec(i) = allScaleSec(j); absErr(i) = d;
end
[idxUnique, ia] = unique(idx, 'stable');
idx = idxUnique;
matchedScaleSec = matchedScaleSec(ia);
requestedScaleSec = requestedScaleSec(ia);
absErr = absErr(ia);
info = struct('requestedScaleSec', requestedScaleSec, 'matchedScaleSec', matchedScaleSec, 'absError', absErr);
end

function meta = local_get_scale_meta(ChunkSet, sIdx)
sc = ChunkSet.scale(sIdx);
if isfield(sc, 'meta') && istable(sc.meta)
    meta = sc.meta;
elseif isfield(sc, 'meta') && isstruct(sc.meta)
    meta = struct2table(sc.meta);
else
    error('build_anchor_multiscale_matrix:MissingMeta', 'ChunkSet.scale(%d).meta table is required.', sIdx);
end
meta = local_standardize_meta(meta);
end

function meta = local_standardize_meta(meta)
vars = string(meta.Properties.VariableNames);
meta = local_rename_if_present(meta, vars, ["anchor_frame","anchorFrame","center_frame","frame"], "anchor_frame");
vars = string(meta.Properties.VariableNames);
meta = local_rename_if_present(meta, vars, ["session_index","sessionIdx","session_index_local","session"], "session_index");
assert(ismember('anchor_frame', meta.Properties.VariableNames), 'Scale meta missing anchor_frame.');
assert(ismember('session_index', meta.Properties.VariableNames), 'Scale meta missing session_index.');
meta.anchor_frame = double(meta.anchor_frame(:));
meta.session_index = double(meta.session_index(:));
end

function T = local_rename_if_present(T, vars, candidates, newName)
if ismember(newName, vars), return; end
hit = intersect(candidates, vars, 'stable');
if ~isempty(hit)
    T.Properties.VariableNames{strcmp(T.Properties.VariableNames, hit(1))} = newName;
end
end

function Y = local_get_scale_embedding(EmbedModel, ChunkSet, sIdx)
Y = [];
if isfield(EmbedModel, 'scale') && numel(EmbedModel.scale) >= sIdx
    st = EmbedModel.scale(sIdx);
    candidates = {'score','Z','Y','embedding','pcs'};
    for i = 1:numel(candidates)
        if isfield(st, candidates{i}) && ~isempty(st.(candidates{i}))
            Y = st.(candidates{i});
            break
        end
    end
end
if isempty(Y)
    sc = ChunkSet.scale(sIdx);
    if isfield(sc, 'Xraw') && ~isempty(sc.Xraw)
        Xraw = sc.Xraw;
        if ndims(Xraw) == 3
            Y = reshape(Xraw, size(Xraw,1), []);
        else
            Y = Xraw;
        end
    else
        error('build_anchor_multiscale_matrix:MissingEmbedding', ...
            'Could not find scale embedding for scale %d.', sIdx);
    end
end
Y = double(Y);
end

function refLocalIdx = local_choose_reference_scale(scaleData, modeIn, refScaleSec)
mode = lower(string(modeIn));
selSec = cellfun(@(s)s.scaleSec, scaleData);
switch mode
    case "median_selected"
        target = median(selSec);
        [~, refLocalIdx] = min(abs(selSec - target));
    case "largest"
        [~, refLocalIdx] = max(selSec);
    case "smallest"
        [~, refLocalIdx] = min(selSec);
    case "user"
        assert(~isempty(refScaleSec), 'ReferenceScaleSec is required when ReferenceScaleMode=''user''.');
        [~, refLocalIdx] = min(abs(selSec - refScaleSec));
    otherwise
        error('build_anchor_multiscale_matrix:BadReferenceScaleMode', ...
            'Unknown ReferenceScaleMode: %s', mode);
end
end

function tol = local_choose_anchor_tolerance(scaleData, refLocalIdx, tolUser)
nSel = numel(scaleData);
tol = nan(nSel,1);
if ~isempty(tolUser)
    tol(:) = tolUser;
    return
end
refFrames = scaleData{refLocalIdx}.meta.anchor_frame(:);
refStep = local_median_positive_diff(refFrames);
if ~isfinite(refStep), refStep = 1; end
for i = 1:nSel
    thisFrames = scaleData{i}.meta.anchor_frame(:);
    thisStep = local_median_positive_diff(thisFrames);
    if ~isfinite(thisStep), thisStep = refStep; end
    tol(i) = max(1, ceil(0.5 * max(refStep, thisStep)));
end
end

function d = local_median_positive_diff(x)
x = unique(double(x(:)));
dx = diff(x);
dx = dx(dx > 0);
if isempty(dx), d = NaN; else, d = median(dx); end
end

function [sessAll, frameAll] = local_build_union_anchor_set(scaleData)
sessAll = [];
frameAll = [];
for i = 1:numel(scaleData)
    sessAll = [sessAll; scaleData{i}.meta.session_index(:)]; %#ok<AGROW>
    frameAll = [frameAll; scaleData{i}.meta.anchor_frame(:)]; %#ok<AGROW>
end
X = [sessAll, frameAll];
[~, ia] = unique(X, 'rows', 'stable');
sessAll = sessAll(ia);
frameAll = frameAll(ia);
end

function [loc, dFrame] = local_match_rows_to_reference(refSess, refFrame, srcSess, srcFrame, alignMode, tolFrames)
nRef = numel(refFrame);
loc = zeros(nRef,1);
dFrame = nan(nRef,1);
if lower(string(alignMode)) == "exact"
    keyRef = string(refSess) + "|" + string(refFrame);
    keySrc = string(srcSess) + "|" + string(srcFrame);
    [tf, idx] = ismember(keyRef, keySrc);
    loc(tf) = idx(tf);
    dFrame(tf) = 0;
    return
end

uSess = unique(refSess(:));
for s = 1:numel(uSess)
    sid = uSess(s);
    ridx = find(refSess == sid);
    sidx = find(srcSess == sid);
    if isempty(sidx), continue; end
    srcF = srcFrame(sidx);
    for j = 1:numel(ridx)
        rf = refFrame(ridx(j));
        [d, k] = min(abs(srcF - rf));
        if d <= tolFrames
            loc(ridx(j)) = sidx(k);
            dFrame(ridx(j)) = srcF(k) - rf;
        end
    end
end
end

function retention = local_choose_retention_rule(nPresent, missingFrac, P, nSel)
rule = lower(string(P.RetentionRule));
if isempty(nPresent)
    minScales = NaN; maxMissing = NaN; retainedMask = false(0,1);
    desc = 'empty';
elseif P.RequireAllSelectedScales
    minScales = nSel; maxMissing = 0;
    retainedMask = (nPresent >= minScales);
    desc = sprintf('require_all_selected_scales (%d/%d)', nSel, nSel);
elseif rule == "manual"
    if isempty(P.MinScalesPresent), minScales = 1; else, minScales = P.MinScalesPresent; end
    if isempty(P.MaxMissingScaleFraction), maxMissing = 1 - minScales ./ max(nSel,1); else, maxMissing = P.MaxMissingScaleFraction; end
    retainedMask = (nPresent >= minScales) & (missingFrac <= maxMissing);
    desc = sprintf('manual | min scales present >= %d | max missing frac <= %.3f', minScales, maxMissing);
else
    q = quantile(nPresent, [0.25 0.5]);
    minScales = max(2, floor(q(1)));
    minScales = min(minScales, nSel);
    maxMissing = min(0.5, 1 - minScales ./ max(nSel,1));
    retainedMask = (nPresent >= minScales) & (missingFrac <= maxMissing);
    desc = sprintf('auto | min scales present >= %d | max missing frac <= %.3f', minScales, maxMissing);
end
retention = struct('minScalesPresent', minScales, ...
    'maxMissingScaleFraction', maxMissing, ...
    'retainedMask', retainedMask, ...
    'description', desc);
end

function q = local_quantiles_safe(x)
if isempty(x)
    q = nan(1,7);
else
    q = quantile(double(x(:)), [0 0.1 0.25 0.5 0.75 0.9 1]);
end
end

function v = local_prctile_safe(x, p)
if isempty(x)
    v = NaN;
else
    v = prctile(double(x(:)), p);
end
end
