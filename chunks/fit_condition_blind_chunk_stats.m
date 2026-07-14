function stats = fit_condition_blind_chunk_stats(X, channelMeta)
%FIT_CONDITION_BLIND_CHUNK_STATS Fit robust transform stats for chunks.
%
% X contains transformed frame-level channels from prepare_dyad_timeseries.
% Rows may be sampled for tractability, but they must be selected without
% condition, cohort, arena, drug, genotype, or outcome labels.

assert(ismatrix(X) && isnumeric(X), ...
    'fit_condition_blind_chunk_stats:BadX', ...
    'X must be a numeric matrix.');
assert(istable(channelMeta) && height(channelMeta) == size(X, 2), ...
    'fit_condition_blind_chunk_stats:BadChannelMeta', ...
    'channelMeta must have one row per X column.');
assert(ismember('ChannelType', channelMeta.Properties.VariableNames), ...
    'fit_condition_blind_chunk_stats:MissingChannelType', ...
    'channelMeta must contain ChannelType.');

D = size(X, 2);
stats = struct();
stats.center = zeros(1, D);
stats.scale = ones(1, D);
stats.impute = zeros(1, D);
stats.n_fit_rows = size(X, 1);
stats.fit_scope = "condition_blind_valid_frame_sample";

for d = 1:D
    x = X(:, d);
    ok = isfinite(x);
    if ~any(ok)
        continue
    end

    if string(channelMeta.ChannelType(d)) == "boolean"
        stats.center(d) = 0;
        stats.scale(d) = 1;
        stats.impute(d) = median(x(ok));
    else
        stats.center(d) = median(x(ok));
        sc = iqr(x(ok));
        if ~(isfinite(sc) && sc > 0)
            sc = std(x(ok), 0);
        end
        if ~(isfinite(sc) && sc > 0)
            sc = 1;
        end
        stats.scale(d) = sc;
        stats.impute(d) = stats.center(d);
    end
end
end
