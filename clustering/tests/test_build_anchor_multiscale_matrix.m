function tests = test_build_anchor_multiscale_matrix
tests = functiontests(localfunctions);
end

function testNearestReferenceAlignment(testCase)
scaleVals = [0.93 1.26 1.72 2.34];
baseFrames = (100:10:200)';
for s = 1:numel(scaleVals)
    ChunkSet.scale(s).chunkSec = scaleVals(s);
    T = table();
    T.session_index = ones(numel(baseFrames),1);
    T.anchor_frame = baseFrames + (s-1); % slightly shifted across scales
    ChunkSet.scale(s).meta = T;
    ChunkSet.scale(s).Xraw = randn(numel(baseFrames), 5, 3);
    EmbedModel.scale(s).score = randn(numel(baseFrames), 4);
end

selectedScales = scaleVals;
Data = build_anchor_multiscale_matrix(ChunkSet, EmbedModel, selectedScales, ...
    'AnchorSetMode', 'reference', ...
    'AnchorAlignment', 'nearest', ...
    'RequireAllSelectedScales', true, ...
    'AnchorToleranceFrames', 3, ...
    'Verbose', false);

verifyEqual(testCase, height(Data.anchorTable), numel(baseFrames));
verifyTrue(testCase, all(Data.anchorTable.n_scales_present == numel(scaleVals)));
end
