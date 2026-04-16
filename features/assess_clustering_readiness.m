function readiness = assess_clustering_readiness(X, sampleTable, featureSelection, opts)
%ASSESS_CLUSTERING_READINESS Assess whether the window feature matrix is ready for clustering.

arguments
    X double
    sampleTable table
    featureSelection struct
    opts.minWindows (1,1) double {mustBeInteger,mustBePositive} = 2000
    opts.maxMedianMissing (1,1) double {mustBeGreaterThanOrEqual(opts.maxMedianMissing,0),mustBeLessThanOrEqual(opts.maxMedianMissing,1)} = 0.05
    opts.minFeatures (1,1) double {mustBeInteger,mustBePositive} = 20
    opts.maxFeatures (1,1) double {mustBeInteger,mustBePositive} = 200
end

if isempty(X)
    readiness = struct('isReady', false, 'reasons', {{'Empty feature matrix.'}}, 'table', table());
    return;
end

rowMissing = mean(isnan(X), 2);
medianMissing = median(rowMissing, 'omitnan');
numWindows = size(X,1);
numFeatures = sum(featureSelection.keepMask);

% Effective dimensionality via PCA on median-imputed selected features.
Xs = featureSelection.Xselected;
Xi = Xs;
for j = 1:size(Xi,2)
    medj = median(Xi(:,j), 'omitnan');
    Xi(isnan(Xi(:,j)), j) = medj;
end

pcNum = min([size(Xi,1), size(Xi,2), 20]);
if pcNum >= 2
    [~,~,latent,~,explained] = pca(Xi, 'NumComponents', pcNum);
    effectiveRank = sum(explained > 1);
    pc1Frac = explained(1);
else
    latent = []; %#ok<NASGU>
    explained = NaN;
    effectiveRank = NaN;
    pc1Frac = NaN;
end

% Session dominance.
if any(strcmp(sampleTable.Properties.VariableNames, 'session_idx'))
    sessionCounts = groupcounts(sampleTable(:, 'session_idx'));
    largestSessionFrac = max(sessionCounts.GroupCount) / sum(sessionCounts.GroupCount);
else
    largestSessionFrac = NaN;
end

checks = {
    'Enough windows', numWindows >= opts.minWindows;
    'Median missingness acceptable', medianMissing <= opts.maxMedianMissing;
    'Feature count in range', numFeatures >= opts.minFeatures && numFeatures <= opts.maxFeatures;
    'PC1 not dominant', isnan(pc1Frac) || pc1Frac < 60;
    'Session dominance acceptable', isnan(largestSessionFrac) || largestSessionFrac < 0.7};

pass = cell2mat(checks(:,2));
reasons = checks(~pass,1);

summary = table(numWindows, medianMissing, numFeatures, effectiveRank, pc1Frac, largestSessionFrac, ...
    'VariableNames', {'NumWindows','MedianRowMissing','NumSelectedFeatures','EffectiveRank','PC1ExplainedPct','LargestSessionFraction'});

readiness = struct();
readiness.isReady = all(pass);
readiness.reasons = reasons;
readiness.table = summary;
readiness.checks = cell2table(checks, 'VariableNames', {'Criterion','Pass'});
end
