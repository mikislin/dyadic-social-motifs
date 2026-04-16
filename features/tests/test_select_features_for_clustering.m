function tests = test_select_features_for_clustering
tests = functiontests(localfunctions);
end

function testFeatureSelection(testCase)
rng(1);
X = randn(500, 10);
X(:,10) = X(:,1) + 0.001*randn(500,1); % redundant with feature 1
names = arrayfun(@(k) sprintf('f%d', k), 1:10, 'UniformOutput', false);
S = select_features_for_clustering(X, names, 'maxKeep', 8, 'minKeep', 5, 'maxPairCorr', 0.95);

verifyGreaterThanOrEqual(testCase, numel(S.keepIdx), 5);
verifyLessThanOrEqual(testCase, numel(S.keepIdx), 8);
verifyFalse(testCase, S.keepMask(1) && S.keepMask(10));
end
