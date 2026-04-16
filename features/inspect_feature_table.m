function report = inspect_feature_table(X, featureNames, opts)
%INSPECT_FEATURE_TABLE Feature inspection stats before clustering.

arguments
    X double
    featureNames cell
    opts.corrType char {mustBeMember(opts.corrType, {'Pearson','Spearman'})} = 'Spearman'
end

assert(size(X,2) == numel(featureNames), 'featureNames must match columns in X.');

n = size(X,1);
p = size(X,2);
missingFrac = mean(isnan(X), 1)';
validN = sum(~isnan(X), 1)';
medianVal = median(X, 1, 'omitnan')';
madVal = mad(X, 1, 1)';
stdVal = std(X, 0, 1, 'omitnan')';
IQRval = iqr(X)';

isBinary = false(p,1);
for j = 1:p
    u = unique(X(~isnan(X(:,j)), j));
    isBinary(j) = ~isempty(u) && all(ismember(u, [0 1]));
end

prevalence = nan(p,1);
for j = 1:p
    if isBinary(j)
        prevalence(j) = mean(X(:,j) > 0.5, 'omitnan');
    end
end

% Pairwise correlations on pairwise-complete rows.
R = corr(X, 'Type', opts.corrType, 'Rows', 'pairwise');
absR = abs(R);
absR(1:p+1:end) = NaN;
maxAbsCorr = max(absR, [], 2, 'omitnan');

report = table(featureNames(:), validN, missingFrac, medianVal, madVal, IQRval, stdVal, isBinary, prevalence, maxAbsCorr, ...
    'VariableNames', {'Feature','ValidN','MissingFrac','Median','MAD','IQR','Std','IsBinary','Prevalence','MaxAbsCorr'});
report.Properties.UserData.corr = R;
report.Properties.UserData.nRows = n;
end
