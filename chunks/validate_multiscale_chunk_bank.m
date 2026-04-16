function Report = validate_multiscale_chunk_bank(ChunkSet, varargin)
%VALIDATE_MULTISCALE_CHUNK_BANK Summarize dense scale-bank extraction quality.

p = inputParser;
p.addParameter('makePlots', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('verbose', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('maxScalesForExamples', 9, @(x)isscalar(x) && x >= 1);
p.parse(varargin{:});
P = p.Results;

nScale = numel(ChunkSet.scale);
rows = table();
issueList = strings(0,1);
for s = 1:nScale
    Sc = ChunkSet.scale(s);
    meta = Sc.meta;
    nChunks = size(Sc.X,1);
    minVF = NaN; meanVF = NaN; medStride = NaN;
    if ~isempty(meta)
        minVF = min(meta.valid_frac); meanVF = mean(meta.valid_frac);
        if height(meta) > 1
            medStride = median(diff(meta.anchor_frame));
        end
        durFrames = unique(meta.stop_frame - meta.start_frame + 1);
        if numel(durFrames) ~= 1 || durFrames ~= Sc.nFrames
            issueList(end+1,1) = sprintf('Scale %.4gs duration mismatch.', Sc.chunkSec); %#ok<AGROW>
        end
    end
    rows = [rows; table(s, Sc.chunkSec, Sc.nFrames, nChunks, meanVF, minVF, medStride, string(ChunkSet.scaleBankMeta.temporal_band(s)), ...
        'VariableNames', {'scale_index','chunk_sec','chunk_frames','n_chunks','mean_valid_frac','min_valid_frac','median_anchor_stride_frames','temporal_band'})]; %#ok<AGROW>
end

% Smoothness / redundancy diagnostics based on adjacent scales.
logSec = log10(rows.chunk_sec);
logGap = [NaN; diff(logSec)];
rows.log10_scale_gap = logGap;
rows.chunk_density_per_sec = rows.n_chunks ./ max(rows.chunk_sec, eps);

% Band summaries
bands = unique(rows.temporal_band, 'stable');
bandRows = table();
for b = 1:numel(bands)
    idx = rows.temporal_band == bands(b);
    bandRows = [bandRows; table(bands(b), nnz(idx), mean(rows.chunk_sec(idx)), mean(rows.n_chunks(idx)), mean(rows.mean_valid_frac(idx)), ...
        'VariableNames', {'temporal_band','n_scales','mean_chunk_sec','mean_n_chunks','mean_valid_frac'})]; %#ok<AGROW>
end

Report = struct();
Report.scaleSummary = rows;
Report.bandSummary = bandRows;
Report.issues = issueList;
Report.isValid = isempty(issueList);

if P.verbose
    disp('=== Dense scale-bank validation: scale summary ==='); disp(rows)
    disp('=== Dense scale-bank validation: temporal-band summary ==='); disp(bandRows)
    if isempty(issueList), disp('No validation issues detected.'); else, disp(issueList); end
end
if P.makePlots
    plot_multiscale_chunk_bank_diagnostics(ChunkSet, Report, 'maxScalesForExamples', P.maxScalesForExamples);
end
end
