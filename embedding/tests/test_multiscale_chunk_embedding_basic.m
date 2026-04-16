function tests = test_multiscale_chunk_embedding_basic
%TEST_MULTISCALE_CHUNK_EMBEDDING_BASIC Unit tests for chunk embedding layer.
tests = functiontests(localfunctions);
end

function testEmbeddingFitsOnMockChunkSet(testCase)
dyad = i_make_mock_dyad(600, 80);
ChunkSet = build_multiscale_chunk_dataset({dyad}, ...
    'chunkSec', [0.4 0.8 1.5], 'strideSec', 0.1, 'anchorMode', 'center');
EmbedModel = fit_multiscale_chunk_embedding(ChunkSet, 'nPCsPerScale', 6, 'globalNPCs', 8, 'verbose', false);
verifyTrue(testCase, isfield(EmbedModel, 'embedding'));
verifyGreaterThan(testCase, size(EmbedModel.embedding,2), 0);
verifyEqual(testCase, size(EmbedModel.embedding,1), height(EmbedModel.chunkTable));
end

function testValidationRuns(testCase)
dyad = i_make_mock_dyad(600, 80);
ChunkSet = build_multiscale_chunk_dataset({dyad}, ...
    'chunkSec', [0.4 0.8 1.5], 'strideSec', 0.1);
EmbedModel = fit_multiscale_chunk_embedding(ChunkSet, 'nPCsPerScale', 5, 'globalNPCs', 6, 'verbose', false);
R = validate_multiscale_chunk_embedding(EmbedModel, 'makePlots', false, 'verbose', false);
verifyTrue(testCase, R.isReadyForClustering);
verifyEqual(testCase, height(R.scaleSummary), 3);
end

function testApplyEmbedding(testCase)
dyad1 = i_make_mock_dyad(600, 80);
dyad2 = i_make_mock_dyad(500, 80);
ChunkSet1 = build_multiscale_chunk_dataset({dyad1}, ...
    'chunkSec', [0.4 0.8], 'strideSec', 0.1);
ChunkSet2 = build_multiscale_chunk_dataset({dyad2}, ...
    'chunkSec', [0.4 0.8], 'strideSec', 0.1);
EmbedModel = fit_multiscale_chunk_embedding(ChunkSet1, 'nPCsPerScale', 4, 'globalNPCs', 5, 'verbose', false);
Out = apply_multiscale_chunk_embedding(ChunkSet2, EmbedModel);
verifyEqual(testCase, size(Out.embedding,1), height(Out.chunkTable));
verifyGreaterThan(testCase, size(Out.embedding,2), 0);
end

function dyad = i_make_mock_dyad(T, fps)
if nargin < 2
    fps = 80;
end
[featureNames, featureMeta] = default_dyad_feature_metadata();
time_s = (0:T-1)' ./ fps;
X = zeros(T, numel(featureNames));
for f = 1:numel(featureNames)
    if featureMeta.IsBoolean(f)
        x = double(rand(T,1) > 0.8);
    elseif featureMeta.IsCircular(f)
        x = mod(cumsum(randn(T,1) * 5), 360) - 180;
    elseif string(featureMeta.TransformHint{f}) == "log1p"
        x = abs(50 + cumsum(randn(T,1) * 1.2));
    else
        x = cumsum(randn(T,1) * 0.15);
    end
    X(:,f) = x;
end
raw = struct();
for f = 1:numel(featureNames)
    raw.(featureNames{f}) = X(:,f);
end
frameMask = true(T,1);
frameMask(randperm(T, round(0.03*T))) = false;

dyad = struct();
dyad.time_s = time_s;
dyad.fps = fps;
dyad.X = X;
dyad.featureNames = featureNames;
dyad.featureMeta = featureMeta;
dyad.raw = raw;
dyad.frameMask = frameMask;
end
