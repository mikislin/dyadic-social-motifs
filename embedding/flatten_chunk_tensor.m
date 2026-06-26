function [Xflat, dimMeta] = flatten_chunk_tensor(Sc, ChunkSet)
%FLATTEN_CHUNK_TENSOR Flatten a multiscale chunk tensor for embedding models.
%   [Xflat, dimMeta] = flatten_chunk_tensor(Sc, ChunkSet) converts
%   Sc.X [N x L x D] to [N x (L*D)], masking invalid time points with NaN.

X = double(Sc.X);
V = logical(Sc.valid);
[nChunks, L, D] = size(X);

for t = 1:L
    bad = ~V(:,t);
    if any(bad)
        X(bad,t,:) = NaN;
    end
end

Xflat = reshape(X, nChunks, L * D);
obsNames = string(ChunkSet.obsNames(:));
channelMeta = ChunkSet.channelMeta;
baseFeature = strings(L*D,1);
channelType = strings(L*D,1);
timeIndex = zeros(L*D,1);
obsName = strings(L*D,1);
family = strings(L*D,1);
idx = 0;

for t = 1:L
    for d = 1:D
        idx = idx + 1;
        obsName(idx) = obsNames(d);
        baseFeature(idx) = string(channelMeta.BaseFeature(d));
        channelType(idx) = string(channelMeta.ChannelType(d));
        row = strcmp(string(ChunkSet.featureNames), string(channelMeta.BaseFeature(d)));
        if any(row)
            family(idx) = string(ChunkSet.featureMeta.Family{find(row,1,'first')});
        else
            family(idx) = "unknown";
        end
        timeIndex(idx) = t;
    end
end

dimMeta = table(obsName, baseFeature, channelType, family, timeIndex, ...
    'VariableNames', {'ObsName','BaseFeature','ChannelType','Family','TimeIndex'});
end
