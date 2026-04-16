function windows = summarize_windows(dyad, fps, windowSec, opts)
%SUMMARIZE_WINDOWS Summarize frame-level dyadic features into multi-scale windows.
%
%   windows = summarize_windows(dyad, fps, [0.2 0.5 1 2])
%
% Output
%   windows.table      one row per valid center frame per window scale
%   windows.meta       metadata for summary columns
%   windows.windowSec  window sizes used
%
% Summary rules
%   numeric features  -> mean, std, min, max, slope, delta
%   circular features -> circmean, resultant_length
%   boolean features  -> occupancy, onset_rate

arguments
    dyad struct
    fps (1,1) double {mustBePositive}
    windowSec (1,:) double {mustBePositive} = [0.2 0.5 1 2]
    opts.requireFullWindow (1,1) logical = true
    opts.maxMissingFrac (1,1) double {mustBeGreaterThanOrEqual(opts.maxMissingFrac,0),mustBeLessThanOrEqual(opts.maxMissingFrac,1)} = 0.10
    opts.sessionIdx (1,1) double = 1
end

X = dyad.X;
T = size(X,1);
featureNames = dyad.featureNames;
meta = dyad.featureMeta;
frameMask = dyad.frameMask(:);
time_s = dyad.time_s(:);

rows = struct([]);
rowCount = 0;
summaryMetaNames = {};
summaryMetaSource = {};
summaryMetaStat = {};
summaryMetaFamily = {};

for w = 1:numel(windowSec)
    winFrames = max(1, round(windowSec(w) * fps));
    halfLeft = floor((winFrames - 1)/2);
    halfRight = ceil((winFrames - 1)/2);

    startIdx = 1 + halfLeft;
    endIdx = T - halfRight;
    if endIdx < startIdx
        continue;
    end

    for t = startIdx:endIdx
        idx = (t-halfLeft):(t+halfRight);
        validFrameFrac = mean(frameMask(idx));
        if opts.requireFullWindow && validFrameFrac < (1 - opts.maxMissingFrac)
            continue;
        end

        rowCount = rowCount + 1;
        S = struct();
        S.session_idx = opts.sessionIdx;
        S.center_frame = t;
        S.center_time_s = time_s(t);
        S.window_sec = windowSec(w);
        S.valid_frame_fraction = validFrameFrac;

        for f = 1:numel(featureNames)
            x = X(idx, f);
            xValid = x(~isnan(x));

            if isempty(xValid)
                stats = nan_stats(meta, f);
            elseif meta.IsBoolean(f)
                stats = bool_stats(x);
            elseif meta.IsCircular(f)
                stats = circ_stats_deg(x);
            else
                stats = scalar_stats(x, fps);
            end

            fn = featureNames{f};
            statNames = fieldnames(stats);
            for k = 1:numel(statNames)
                colName = matlab.lang.makeValidName(sprintf('%s__%s__%0.1fs', fn, statNames{k}, windowSec(w)));
                S.(colName) = stats.(statNames{k});
                if rowCount == 1
                    summaryMetaNames{end+1,1} = colName; %#ok<AGROW>
                    summaryMetaSource{end+1,1} = fn; %#ok<AGROW>
                    summaryMetaStat{end+1,1} = statNames{k}; %#ok<AGROW>
                    summaryMetaFamily{end+1,1} = string(meta.Family{f}); %#ok<AGROW>
                end
            end
        end

        rows = append_row(rows, S); %#ok<AGROW>
    end
end

if isempty(rows)
    windows.table = table();
    windows.meta = table();
    windows.windowSec = windowSec;
    return;
end

windows.table = struct2table(rows);
windows.meta = table(summaryMetaNames, summaryMetaSource, summaryMetaStat, summaryMetaFamily, ...
    'VariableNames', {'SummaryName','SourceFeature','Statistic','Family'});
windows.windowSec = windowSec;
end

function rows = append_row(rows, S)
%APPEND_ROW Append struct row while harmonizing fields across rows.

if isempty(rows)
    rows = S;
    return;
end

rowsFields = fieldnames(rows);
sFields = fieldnames(S);

allFields = union(rowsFields, sFields, 'stable');

% Add missing fields to existing rows
for i = 1:numel(allFields)
    fn = allFields{i};
    if ~isfield(rows, fn)
        [rows.(fn)] = deal(NaN);
    end
end

% Add missing fields to new struct
for i = 1:numel(allFields)
    fn = allFields{i};
    if ~isfield(S, fn)
        S.(fn) = NaN;
    end
end

% Reorder S fields to match existing struct array order
S = orderfields(S, rows);

rows(end+1,1) = S;
end

function stats = scalar_stats(x, fps)
x = x(:);
idx = ~isnan(x);
xi = x(idx);
if isempty(xi)
    stats = nan_stats();
    return;
end
stats.mean = mean(xi);
stats.std = std(xi, 0);
stats.min = min(xi);
stats.max = max(xi);
stats.delta = xi(end) - xi(1);
if numel(xi) >= 2
    tt = (0:numel(xi)-1)' ./ fps;
    p = polyfit(tt, xi, 1);
    stats.slope = p(1);
else
    stats.slope = NaN;
end
end

function stats = bool_stats(x)
x = x(:);
idx = ~isnan(x);
xi = x(idx);
if isempty(xi)
    stats.occupancy = NaN;
    stats.onset_rate = NaN;
    return;
end
xi = xi > 0.5;
stats.occupancy = mean(xi);
onsets = sum(diff([false; xi]) == 1);
stats.onset_rate = onsets / max(numel(xi),1);
end

function stats = circ_stats_deg(x)
x = x(:);
xi = x(~isnan(x));
if isempty(xi)
    stats.circmean = NaN;
    stats.resultant = NaN;
    return;
end
r = deg2rad(xi);
C = mean(cos(r));
S = mean(sin(r));
stats.circmean = rad2deg(atan2(S, C));
stats.resultant = hypot(C, S);
end

function stats = nan_stats(meta, f)
if nargin < 2
    stats.mean = NaN; stats.std = NaN; stats.min = NaN; stats.max = NaN; stats.delta = NaN; stats.slope = NaN;
    return;
end
if meta.IsBoolean(f)
    stats.occupancy = NaN; stats.onset_rate = NaN;
elseif meta.IsCircular(f)
    stats.circmean = NaN; stats.resultant = NaN;
else
    stats.mean = NaN; stats.std = NaN; stats.min = NaN; stats.max = NaN; stats.delta = NaN; stats.slope = NaN;
end
end
