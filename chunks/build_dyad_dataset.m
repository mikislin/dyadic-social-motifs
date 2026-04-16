function DyadSet = build_dyad_dataset(dyadCell, varargin)
%BUILD_DYAD_DATASET Build a consistent dataset from one or more dyad structs.
%
% Input
%   dyadCell can be:
%     - cell array of dyad structs
%     - struct array of dyad structs
%
% Output
%   DyadSet.seqs       : cell array of prepared sequences
%   DyadSet.stats      : global scaling stats fit on all sequences
%   DyadSet.obsNames   : observation names after circular expansion
%   DyadSet.channelMeta: observation metadata

p = inputParser;
p.addParameter('anchorMasks', {}, @(x)iscell(x) || islogical(x) || isnumeric(x));
p.addParameter('sessionIds', [], @(x)isnumeric(x) || iscell(x));
p.parse(varargin{:});
P = p.Results;

if ~iscell(dyadCell)
    dyadCell = num2cell(dyadCell);
end
nSeq = numel(dyadCell);

anchorMasks = P.anchorMasks;
if isempty(anchorMasks)
    anchorMasks = cell(nSeq,1);
    for i = 1:nSeq
        if isfield(dyadCell{i}, 'frameMask') && ~isempty(dyadCell{i}.frameMask)
            anchorMasks{i} = logical(dyadCell{i}.frameMask(:));
        else
            anchorMasks{i} = [];
        end
    end
elseif ~iscell(anchorMasks)
    tmp = anchorMasks;
    anchorMasks = cell(nSeq,1);
    for i = 1:nSeq
        anchorMasks{i} = tmp;
    end
end


if isempty(P.sessionIds)
    sessionIds = num2cell(1:nSeq);
elseif iscell(P.sessionIds)
    sessionIds = P.sessionIds;
else
    sessionIds = num2cell(P.sessionIds);
end

% First pass: collect raw transformed channels to fit global stats.
seqsRaw = cell(nSeq,1);
Xall = [];
for i = 1:nSeq
    seqsRaw{i} = prepare_dyad_timeseries(dyadCell{i}, 'anchorMask', anchorMasks{i});
    Xall = [Xall; seqsRaw{i}.X(seqsRaw{i}.validMask,:)]; %#ok<AGROW>
end

channelMeta = seqsRaw{1}.channelMeta;
stats = local_fit_stats(Xall, channelMeta);

% Second pass: apply global stats consistently.
seqs = cell(nSeq,1);
for i = 1:nSeq
    seqs{i} = prepare_dyad_timeseries(dyadCell{i}, ...
        'anchorMask', anchorMasks{i}, ...
        'stats', stats);
    seqs{i}.sessionId = sessionIds{i};
end

DyadSet = struct();
DyadSet.seqs = seqs;
DyadSet.stats = stats;
DyadSet.obsNames = seqs{1}.obsNames;
DyadSet.channelMeta = channelMeta;
DyadSet.featureNames = dyadCell{1}.featureNames;
DyadSet.featureMeta = dyadCell{1}.featureMeta;
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
