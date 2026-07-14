function [Xsummary, dimMeta, audit] = summarize_multiresolution_chunks(Sc, ChunkSet, varargin)
%SUMMARIZE_MULTIRESOLUTION_CHUNKS Build compact condition-blind chunk summaries.
%
% The summary is designed for run_06 scale scoring. It replaces raw
% frame-by-frame flattening with interpretable per-channel history features:
% distribution summaries, temporal-bin means, slopes, endpoint change,
% optional low-frequency DCT coefficients, and boolean transition rates.
% It uses only chunk tensors and frame-valid masks; provenance labels never
% enter the representation.

p = inputParser;
p.addParameter('featureNames', strings(0,1), @(x)isstring(x) || iscell(x));
p.addParameter('nTemporalBins', 6, @(x)isscalar(x) && x >= 1);
p.addParameter('nDctCoeffs', 4, @(x)isscalar(x) && x >= 0);
p.addParameter('includeDct', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('useScaledFeatures', true, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});
P = p.Results;

assert(isfield(Sc, 'valid') && isfield(Sc, 'X') && isfield(Sc, 'Xraw'), ...
    'summarize_multiresolution_chunks:BadScale', ...
    'Scale struct must contain X, Xraw, and valid fields.');
assert(isfield(ChunkSet, 'channelMeta') && istable(ChunkSet.channelMeta), ...
    'summarize_multiresolution_chunks:MissingChannelMeta', ...
    'ChunkSet.channelMeta table is required.');

featureNames = string(P.featureNames);
if isempty(featureNames)
    featureNames = string(ChunkSet.featureNames(:));
end

channelMetaAll = ChunkSet.channelMeta;
channelBase = string(channelMetaAll.BaseFeature);
channelIdx = find(ismember(channelBase, featureNames));
assert(~isempty(channelIdx), 'summarize_multiresolution_chunks:NoChannels', ...
    'No transformed channels matched the requested featureNames.');

V = logical(Sc.valid);
[nChunks, L, nAllChannels] = size(Sc.X);
assert(size(V, 1) == nChunks && size(V, 2) == L, ...
    'summarize_multiresolution_chunks:MaskSizeMismatch', ...
    'Scale valid mask size does not match chunk tensor size.');
assert(max(channelIdx) <= nAllChannels, ...
    'summarize_multiresolution_chunks:ChannelIndexOutOfBounds', ...
    'ChunkSet.channelMeta has more channels than the scale tensor.');

nBins = max(1, floor(P.nTemporalBins));
nDct = max(0, floor(P.nDctCoeffs));
includeDct = logical(P.includeDct) && nDct > 0;

nChannels = numel(channelIdx);
baseKinds = ["mean","std","median","q10","q90","slope","delta_late_minus_early"];
nBase = numel(baseKinds);
nBool = nnz(string(channelMetaAll.ChannelType(channelIdx)) == "boolean");
nDims = nChannels * (nBase + nBins + includeDct * nDct) + nBool;

Xsummary = nan(nChunks, nDims);
dimRows = table();
col = 0;

binEdges = unique(round(linspace(1, L + 1, nBins + 1)));
if numel(binEdges) < nBins + 1
    binEdges = 1:(L + 1);
    nBins = numel(binEdges) - 1;
end
relT = linspace(-1, 1, L);
relT = relT - mean(relT);
slopeDen = sum(relT .^ 2);
if slopeDen <= 0
    slopeDen = 1;
end
dctBasis = local_dct_basis(L, nDct);

for c = 1:nChannels
    ch = channelIdx(c);
    if logical(P.useScaledFeatures)
        Xc = double(Sc.X(:, :, ch));
        sourceTensor = "condition_blind_scaled";
    else
        Xc = double(Sc.Xraw(:, :, ch));
        sourceTensor = "transformed_unscaled";
    end
    Xc(~V) = NaN;
    validFraction = mean(isfinite(Xc), 2);

    rowMean = mean(Xc, 2, 'omitnan');
    rowStd = std(Xc, 0, 2, 'omitnan');
    rowMedian = median(Xc, 2, 'omitnan');
    rowMean(~isfinite(rowMean)) = 0;
    rowStd(~isfinite(rowStd)) = 0;
    rowMedian(~isfinite(rowMedian)) = rowMean(~isfinite(rowMedian));

    Xfill = Xc;
    for i = 1:nChunks
        bad = ~isfinite(Xfill(i, :));
        if any(bad)
            Xfill(i, bad) = rowMedian(i);
        end
    end

    q10 = prctile(Xfill, 10, 2);
    q90 = prctile(Xfill, 90, 2);
    slope = (Xfill * relT(:)) ./ slopeDen;
    earlyIdx = 1:max(1, floor(L / 3));
    lateIdx = max(1, L - floor(L / 3) + 1):L;
    deltaLateEarly = mean(Xfill(:, lateIdx), 2) - mean(Xfill(:, earlyIdx), 2);

    baseVals = {rowMean, rowStd, rowMedian, q10, q90, slope, deltaLateEarly};
    for k = 1:nBase
        col = col + 1;
        Xsummary(:, col) = baseVals{k};
        dimRows = [dimRows; local_dim_row(col, Sc, channelMetaAll, ch, baseKinds(k), ...
            NaN, NaN, sourceTensor, validFraction)]; %#ok<AGROW>
    end

    for b = 1:nBins
        idx = binEdges(b):(binEdges(b + 1) - 1);
        idx = idx(idx >= 1 & idx <= L);
        col = col + 1;
        Xsummary(:, col) = mean(Xc(:, idx), 2, 'omitnan');
        miss = ~isfinite(Xsummary(:, col));
        if any(miss)
            Xsummary(miss, col) = rowMedian(miss);
        end
        dimRows = [dimRows; local_dim_row(col, Sc, channelMetaAll, ch, "temporal_bin_mean", ...
            b, NaN, sourceTensor, validFraction)]; %#ok<AGROW>
    end

    if includeDct
        Xcentered = Xfill - mean(Xfill, 2);
        coeffs = Xcentered * dctBasis';
        for k = 1:nDct
            col = col + 1;
            Xsummary(:, col) = coeffs(:, k);
            dimRows = [dimRows; local_dim_row(col, Sc, channelMetaAll, ch, "low_frequency_dct", ...
                NaN, k, sourceTensor, validFraction)]; %#ok<AGROW>
        end
    end

    if string(channelMetaAll.ChannelType(ch)) == "boolean"
        col = col + 1;
        xb = Xfill > 0.5;
        Xsummary(:, col) = sum(abs(diff(double(xb), 1, 2)), 2) ./ max(L - 1, 1);
        dimRows = [dimRows; local_dim_row(col, Sc, channelMetaAll, ch, "boolean_transition_rate", ...
            NaN, NaN, sourceTensor, validFraction)]; %#ok<AGROW>
    end
end

Xsummary = Xsummary(:, 1:col);
dimMeta = dimRows;

audit = table();
audit.scale_index = local_scale_index(Sc);
audit.chunk_sec = Sc.chunkSec;
audit.chunk_frames = L;
audit.n_chunks = nChunks;
audit.n_input_channels = nChannels;
audit.n_raw_flattened_dims = L * nChannels;
audit.n_summary_dims = size(Xsummary, 2);
audit.compression_ratio_raw_to_summary = audit.n_raw_flattened_dims ./ max(audit.n_summary_dims, 1);
audit.n_temporal_bins = nBins;
audit.n_dct_coeffs = includeDct * nDct;
audit.n_boolean_transition_channels = nBool;
audit.summary_method = "multiresolution_channel_statistics_temporal_bins_low_frequency_dct";
audit.source_tensor = sourceTensor;
audit.uses_frame_mask = true;
audit.labels_used_for_summary = "none";
end

function row = local_dim_row(col, Sc, channelMeta, ch, kind, temporalBin, dctCoeff, sourceTensor, validFraction)
row = table();
row.summary_dim_index = col;
row.scale_index = local_scale_index(Sc);
row.chunk_sec = Sc.chunkSec;
row.obs_name = string(channelMeta.ObsName(ch));
row.base_feature = string(channelMeta.BaseFeature(ch));
row.channel_type = string(channelMeta.ChannelType(ch));
row.summary_kind = string(kind);
row.temporal_bin = temporalBin;
row.dct_coefficient = dctCoeff;
row.source_tensor = sourceTensor;
row.mean_frame_valid_fraction_for_channel = mean(validFraction, 'omitnan');
row.uses_frame_mask = true;
row.labels_used_for_summary = "none";
end

function scaleIndex = local_scale_index(Sc)
scaleIndex = NaN;
if isfield(Sc, 'meta') && istable(Sc.meta) && ismember('scale_index', Sc.meta.Properties.VariableNames) && ~isempty(Sc.meta)
    scaleIndex = Sc.meta.scale_index(1);
end
end

function B = local_dct_basis(L, nCoeff)
if nCoeff <= 0
    B = zeros(0, L);
    return
end
n = 0:(L - 1);
B = zeros(nCoeff, L);
for k = 1:nCoeff
    B(k, :) = sqrt(2 / L) .* cos(pi .* (n + 0.5) .* k ./ L);
end
end
