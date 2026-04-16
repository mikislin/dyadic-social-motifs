function tests = test_adaptive_confidence_thresholds
tests = functiontests(localfunctions);
end

function testAdaptiveThresholdPerNode(~)
scores = [linspace(0,1,100)' linspace(0.2,0.8,100)'];
P = default_preprocessing_params();
[thr, info] = derive_confidence_thresholds(scores, P.confidence);
assert(numel(thr) == 2);
assert(all(thr >= P.confidence.min_threshold));
assert(all(isfinite(info.maskFractionFinal)));
end

function testConfidenceMaskRelaxation(~)
scores = 0.001 * ones(100,3);
P = default_preprocessing_params();
[thr, info] = derive_confidence_thresholds(scores, P.confidence);
assert(all(abs(thr - P.confidence.relaxed_threshold) < 1e-12));
assert(all(info.wasRelaxed));
end
