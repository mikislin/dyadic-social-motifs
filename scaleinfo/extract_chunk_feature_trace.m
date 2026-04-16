function x = extract_chunk_feature_trace(ChunkSet, scaleIdx, chunkIdx, featureName, varargin)
%EXTRACT_CHUNK_FEATURE_TRACE Return one feature trace for one chunk.
%
% x is returned in a human-readable feature space:
% - continuous / boolean features: direct channel values
% - circular features: reconstructed angle in degrees
%
% Optional:
%   'invertLog1p' : true/false, default false

p = inputParser;
p.addParameter('invertLog1p', false, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});
P = p.Results;

featureName = string(featureName);

cm = ChunkSet.channelMeta;
base = string(cm.BaseFeature);
ctype = string(cm.ChannelType);
obsNames = string(ChunkSet.obsNames);

fMeta = ChunkSet.featureMeta;
vn = string(fMeta.Properties.VariableNames);

% Robustly detect the feature-name column
if any(vn == "Name")
    metaNames = string(fMeta.Name);
elseif any(vn == "featureNames")
    metaNames = string(fMeta.featureNames);
elseif any(vn == "FeatureName")
    metaNames = string(fMeta.FeatureName);
else
    error('extract_chunk_feature_trace:MissingNameColumn', ...
        'ChunkSet.featureMeta has no recognized feature-name column.');
end

fRow = find(metaNames == featureName, 1);
if isempty(fRow)
    error('extract_chunk_feature_trace:FeatureNotFound', ...
        'Feature not found in ChunkSet.featureMeta: %s', featureName);
end

isCirc = logical(fMeta.IsCircular(fRow));
isBool = logical(fMeta.IsBoolean(fRow));
hint = string(fMeta.TransformHint(fRow));

X = ChunkSet.scale(scaleIdx).Xraw;

if isCirc
    idxCos = find(base == featureName & ctype == "circular_cos", 1);
    idxSin = find(base == featureName & ctype == "circular_sin", 1);

    assert(~isempty(idxCos) && ~isempty(idxSin), ...
        'Missing circular cos/sin channels for %s.', featureName);

    c = squeeze(X(chunkIdx,:,idxCos));
    s = squeeze(X(chunkIdx,:,idxSin));
    x = atan2d(s, c);
else
    idx = find(base == featureName, 1);

    if isempty(idx)
        idx = find(obsNames == featureName, 1);
    end

    assert(~isempty(idx), 'Missing channel for feature %s.', featureName);

    x = squeeze(X(chunkIdx,:,idx));

    if P.invertLog1p && hint == "log1p"
        x = expm1(x);
    end

    if isBool
        x = double(x > 0.5);
    end
end

x = x(:);
end
