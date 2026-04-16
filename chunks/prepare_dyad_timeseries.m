function Seq = prepare_dyad_timeseries(dyad, varargin)
%PREPARE_DYAD_TIMESERIES Convert one dyad struct into a model-ready raw time-series matrix.
%
% Supports both:
%   A) simulator-style dyad with:
%        dyad.time
%        dyad.<featureName>
%
%   B) compute_dyad_features-style dyad with:
%        dyad.time_s
%        dyad.X
%        dyad.featureNames
%        dyad.featureMeta
%        dyad.raw
%        dyad.frameMask

p = inputParser;
p.addParameter('anchorMask', [], @(x)islogical(x) || isnumeric(x));
p.addParameter('transformConfig', struct(), @isstruct);
p.addParameter('stats', struct(), @isstruct);
p.parse(varargin{:});
P = p.Results;

assert(isfield(dyad, 'featureNames'), 'dyad.featureNames is required.');
assert(isfield(dyad, 'featureMeta'), 'dyad.featureMeta is required.');

% -------------------------------------------------------------------------
% Time handling: support dyad.time or dyad.time_s
% -------------------------------------------------------------------------
if isfield(dyad, 'time')
    time = dyad.time(:);
elseif isfield(dyad, 'time_s')
    time = dyad.time_s(:);
else
    error('prepare_dyad_timeseries:MissingTime', ...
        'dyad.time or dyad.time_s is required.');
end

featureNames = dyad.featureNames(:);
featureMeta = dyad.featureMeta;
T = numel(time);

% -------------------------------------------------------------------------
% FPS handling
% -------------------------------------------------------------------------
if isfield(dyad, 'fps') && ~isempty(dyad.fps)
    fps = dyad.fps;
else
    dt = median(diff(time), 'omitnan');
    fps = 1 / max(dt, eps);
end

% -------------------------------------------------------------------------
% Valid mask handling
% anchorMask overrides frameMask if provided
% -------------------------------------------------------------------------
if isempty(P.anchorMask)
    if isfield(dyad, 'frameMask') && ~isempty(dyad.frameMask)
        validMask = logical(dyad.frameMask(:));
    else
        validMask = true(T,1);
    end
else
    validMask = logical(P.anchorMask(:));
end
assert(numel(validMask) == T, 'valid/anchor mask length mismatch.');

% -------------------------------------------------------------------------
% Feature extraction helper:
% 1) top-level field dyad.(fname)
% 2) dyad.raw.(fname)
% 3) dyad.X + featureNames lookup
% -------------------------------------------------------------------------
X = [];
obsNames = strings(0,1);
baseFeature = strings(0,1);
channelType = strings(0,1);
transformHint = strings(0,1);

for f = 1:numel(featureNames)
    fname = featureNames{f};
    x = local_get_feature_vector(dyad, fname, f, T);
    meta = featureMeta(f,:);
    hint = string(meta.TransformHint);

    if meta.IsCircular
        th = deg2rad(x);
        X = [X, cos(th), sin(th)]; %#ok<AGROW>
        obsNames = [obsNames; string(fname) + "_cos"; string(fname) + "_sin"]; %#ok<AGROW>
        baseFeature = [baseFeature; string(fname); string(fname)]; %#ok<AGROW>
        channelType = [channelType; "circular_cos"; "circular_sin"]; %#ok<AGROW>
        transformHint = [transformHint; hint; hint]; %#ok<AGROW>
    else
        if hint == "log1p"
            x(x < 0) = NaN;
            x = log1p(x);
        end
        X = [X, x]; %#ok<AGROW>
        obsNames = [obsNames; string(fname)]; %#ok<AGROW>
        baseFeature = [baseFeature; string(fname)]; %#ok<AGROW>

        if meta.IsBoolean
            channelType = [channelType; "boolean"]; %#ok<AGROW>
        else
            channelType = [channelType; "continuous"]; %#ok<AGROW>
        end
        transformHint = [transformHint; hint]; %#ok<AGROW>
    end
end

channelMeta = table(obsNames, baseFeature, channelType, transformHint, ...
    'VariableNames', {'ObsName','BaseFeature','ChannelType','TransformHint'});

if isempty(fieldnames(P.stats))
    stats = local_fit_stats(X, channelMeta);
else
    stats = P.stats;
end

Xscaled = local_apply_stats(X, channelMeta, stats);

Seq = struct();
Seq.X = X;
Seq.Xscaled = Xscaled;
Seq.obsNames = cellstr(obsNames);
Seq.channelMeta = channelMeta;
Seq.validMask = validMask;
Seq.time = time;
Seq.fps = fps;
Seq.featureNames = dyad.featureNames;
Seq.featureMeta = dyad.featureMeta;
Seq.stats = stats;
end

function x = local_get_feature_vector(dyad, fname, fIdx, T)
if isfield(dyad, fname)
    x = dyad.(fname);
elseif isfield(dyad, 'raw') && isstruct(dyad.raw) && isfield(dyad.raw, fname)
    x = dyad.raw.(fname);
elseif isfield(dyad, 'X') && ~isempty(dyad.X)
    assert(size(dyad.X,2) >= fIdx, ...
        'dyad.X does not contain expected feature column %d for %s.', fIdx, fname);
    x = dyad.X(:, fIdx);
else
    error('prepare_dyad_timeseries:MissingFeature', ...
        'Could not find feature "%s" in dyad struct.', fname);
end

x = x(:);
assert(numel(x) == T, 'Feature %s length mismatch.', fname);
end

function stats = local_fit_stats(X, channelMeta)
D = size(X,2);
stats = struct();
stats.center = zeros(1,D);
stats.scale = ones(1,D);
stats.impute = zeros(1,D);

for d = 1:D
    x = X(:,d);
    ok = isfinite(x);
    if ~any(ok)
        continue
    end

    if channelMeta.ChannelType(d) == "boolean"
        stats.center(d) = 0;
        stats.scale(d) = 1;
        stats.impute(d) = median(x(ok));
    else
        stats.center(d) = median(x(ok));
        s = iqr(x(ok));
        if ~(isfinite(s) && s > 0)
            s = std(x(ok), 0);
        end
        if ~(isfinite(s) && s > 0)
            s = 1;
        end
        stats.scale(d) = s;
        stats.impute(d) = stats.center(d);
    end
end
end

function Xscaled = local_apply_stats(X, channelMeta, stats)
Xscaled = X;
for d = 1:size(X,2)
    x = X(:,d);
    x(~isfinite(x)) = stats.impute(d);

    if channelMeta.ChannelType(d) == "boolean"
        Xscaled(:,d) = x;
    else
        Xscaled(:,d) = (x - stats.center(d)) ./ stats.scale(d);
    end
end
end