function Data = prepare_run10_independent_feature_matrix( ...
        Measurement, Registry)
%PREPARE_RUN10_INDEPENDENT_FEATURE_MATRIX Select the frozen feature panel.

featureNames = string(Registry.feature_name( ...
    Registry.included_in_primary_validation));
T = Measurement(Measurement.eligible_for_automated_analysis, :);
T = sortrows(T, 'graph_node_id');
X = double(T{:, cellstr(featureNames)});
med = median(X, 1, 'omitnan');
assert(all(isfinite(med)), ...
    'prepare_run10_independent_feature_matrix:AllMissingFeature', ...
    'An included validation feature has no finite observations.');
Ximputed = X;
for j = 1:size(X,2)
    Ximputed(~isfinite(Ximputed(:,j)), j) = med(j);
end
madValue = median(abs(Ximputed - med), 1, 'omitnan');
scale = 1.4826 .* madValue;
fallback = std(Ximputed, 0, 1, 'omitnan');
bad = ~isfinite(scale) | scale <= 1e-12;
scale(bad) = fallback(bad);
scale(~isfinite(scale) | scale <= 1e-12) = 1;
Z = (Ximputed - med) ./ scale;

Data = struct();
Data.table = T;
Data.feature_names = featureNames;
Data.raw = X;
Data.global_imputed = Ximputed;
Data.global_center = med;
Data.global_scale = scale;
Data.global_standardized = Z;
Data.missing_fraction = mean(~isfinite(X), 1);
end
