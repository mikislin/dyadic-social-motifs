function tests = test_dense_chunk_bank_basic
tests = functiontests(localfunctions);
end

function testMakeScaleBank(testCase)
sc = make_logscale_chunk_bank('minSec',0.2,'maxSec',8.0,'nScales',25);
verifyGreaterThanOrEqual(testCase, numel(sc), 20);
verifyEqual(testCase, sc(1), 0.2, 'AbsTol', 1e-4);
verifyEqual(testCase, sc(end), 8.0, 'AbsTol', 1e-4);
verifyTrue(testCase, all(diff(sc) > 0));
end

function testDenseChunkDatasetBuilds(testCase)
dyad = i_make_mock_dyad(2000, 80);
ChunkSet = build_multiscale_chunk_dataset({dyad}, 'useLogScaleBank', true, 'minScaleSec', 0.2, 'maxScaleSec', 8.0, 'nScales', 25, 'strideSec', 0.1);
verifyEqual(testCase, numel(ChunkSet.scale), numel(ChunkSet.chunkSec));
verifyGreaterThanOrEqual(testCase, numel(ChunkSet.chunkSec), 20);
verifyTrue(testCase, height(ChunkSet.chunkTable) > 0);
end

function testDenseValidationPasses(testCase)
dyad = i_make_mock_dyad(1500, 80);
ChunkSet = build_multiscale_chunk_dataset({dyad}, 'useLogScaleBank', true, 'strideSec', 0.1);
R = validate_multiscale_chunk_bank(ChunkSet, 'makePlots', false, 'verbose', false);
verifyTrue(testCase, R.isValid);
verifyEqual(testCase, height(R.scaleSummary), numel(ChunkSet.scale));
end

function dyad = i_make_mock_dyad(T, fps)
if nargin < 2, fps = 80; end
[featureNames, featureMeta] = default_dyad_feature_metadata();
time_s = (0:T-1)' ./ fps;
X = zeros(T, numel(featureNames));
for f = 1:numel(featureNames)
    if featureMeta.IsBoolean(f)
        x = double(rand(T,1) > 0.8);
    elseif featureMeta.IsCircular(f)
        x = mod(cumsum(randn(T,1) * 4), 360) - 180;
    elseif featureMeta.TransformHint(f) == "log1p"
        x = abs(50 + cumsum(randn(T,1)));
    else
        x = cumsum(randn(T,1) * 0.25);
    end
    X(:,f) = x;
end
raw = struct();
for f = 1:numel(featureNames), raw.(featureNames{f}) = X(:,f); end
frameMask = true(T,1); frameMask(randperm(T, round(0.03*T))) = false;
dyad = struct('time_s',time_s,'fps',fps,'X',X,'featureNames',{featureNames},'featureMeta',featureMeta,'raw',raw,'frameMask',frameMask);
end
