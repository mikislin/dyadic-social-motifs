function tests = test_04_dyad_feature_metadata
tests = functiontests(localfunctions);
end

function testFeatureDictionaryHasReviewerFields(testCase)
[featureNames, featureMeta] = default_dyad_feature_metadata();

required = ["FeatureIndex","Name","Family","FamilyLabel","Unit","UnitLabel", ...
    "FeatureFamilyRole","FamilyDefinition","IsDirected","Directionality", ...
    "FocalAnimal","ReferenceFrame", ...
    "IsCircular","IsBoolean","TransformHint","ExpectedRange", ...
    "BiologicalInterpretation","ClusteringCandidate","FeatureLayerRole", ...
    "WindowSummaryStats"];
verifyTrue(testCase, all(ismember(required, string(featureMeta.Properties.VariableNames))));
verifyEqual(testCase, string(featureNames(:)), string(featureMeta.Name));
verifyEqual(testCase, numel(unique(featureMeta.Name)), height(featureMeta));
verifyTrue(testCase, all(strlength(featureMeta.Unit) > 0));
verifyTrue(testCase, all(strlength(featureMeta.TransformHint) > 0));
verifyTrue(testCase, all(strlength(featureMeta.BiologicalInterpretation) > 0));
verifyTrue(testCase, all(featureMeta.FeatureFamilyRole == ...
    "computational_measurement_group_not_biological_motif_family"));
verifyTrue(testCase, all(strlength(featureMeta.FamilyDefinition) > 0));
verifyTrue(testCase, all(featureMeta.ClusteringCandidate));
verifyTrue(testCase, all(featureMeta.FeatureLayerRole == "condition_blind_frame_feature_for_chunking"));
verifyTrue(testCase, all(strlength(featureMeta.WindowSummaryStats) > 0));
verifyTrue(testCase, all(featureMeta.Unit(featureMeta.IsCircular) == "degrees"));
verifyTrue(testCase, all(featureMeta.Unit(featureMeta.IsBoolean) == "boolean"));
verifyTrue(testCase, all(featureMeta.WindowSummaryStats(featureMeta.IsCircular) == "circmean;resultant"));
verifyTrue(testCase, all(featureMeta.WindowSummaryStats(featureMeta.IsBoolean) == "occupancy;onset_rate"));
end

function testCanonicalNodeMapComesFromSleapMetadata(testCase)
[nodeMap, partNames, bodypartMeta] = default_sleap_node_map();

verifyEqual(testCase, partNames, cellstr(bodypartMeta.bodypart_name)');
verifyEqual(testCase, nodeMap.nose, 1);
verifyEqual(testCase, nodeMap.neck, 2);
verifyEqual(testCase, nodeMap.body, 9);
verifyEqual(testCase, nodeMap.tailBase, 10);
verifyEqual(testCase, nodeMap.midBody, 13);
verifyEqual(testCase, bodypartMeta.node_index', 1:13);
end

function testFeatureExtractionConfigDeclaresUnitScale(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
params = load_motif_dyad_feature_extraction_config( ...
    fullfile(repoRoot, 'config', 'motif_dyad_feature_extraction_config.csv'));

verifyEqual(testCase, params.fps, 80);
verifyEqual(testCase, params.pixel_size_mm, 1/1.97, 'AbsTol', 1e-12);
verifyGreaterThanOrEqual(testCase, params.unit_parameter_tolerance, 0);
end

function testFeatureUnitsUsePixelSizeMM(testCase)
T = 20;
tracks = local_mock_tracks(T, 20);
nodeMap = default_sleap_node_map();

D = compute_dyad_features(tracks, 80, nodeMap, ...
    'pixelSizeMM', 0.5, ...
    'smoothSpanFrames', 1);

centroidIdx = strcmp(D.featureNames, 'centroid_dist');
verifyTrue(testCase, any(centroidIdx));
verifyEqual(testCase, D.X(:, centroidIdx), repmat(10, T, 1), 'AbsTol', 1e-12);
verifyEqual(testCase, string(D.featureMeta.Unit(centroidIdx)), "mm");
verifyEqual(testCase, D.params.pixelSizeMM, 0.5);
verifyTrue(testCase, all(D.frameMask));
end

function testBadframesBecomeInvalidAndMissing(testCase)
T = 24;
tracks = local_mock_tracks(T, 20);
nodeMap = default_sleap_node_map();
badframes = false(T, 2);
badframes(5:7, 2) = true;

D = compute_dyad_features(tracks, 80, nodeMap, ...
    'pixelSizeMM', 0.5, ...
    'smoothSpanFrames', 1, ...
    'badframes', badframes);

verifyFalse(testCase, any(D.frameMask(5:7)));
verifyTrue(testCase, all(all(isnan(D.X(5:7,:)))));
verifyEqual(testCase, D.maskAudit.nBadframeRows, 3);
verifyEqual(testCase, D.maskAudit.nCoreNanOnlyFrames, 0);
verifyEqual(testCase, D.badframeMask, any(badframes, 2));
end

function testNoncanonicalNodeMapRejected(testCase)
T = 12;
tracks = local_mock_tracks(T, 20);
nodeMap = default_sleap_node_map();
nodeMap.nose = 2;

verifyError(testCase, @() compute_dyad_features(tracks, 80, nodeMap), ...
    'compute_dyad_features:NoncanonicalNodeMap');
end

function tracks = local_mock_tracks(T, offsetPx)
N = 13;
tracks = nan(T, N, 2, 2);

tracks(:,1,:,1) = repmat([10 0], T, 1);  % nose
tracks(:,2,:,1) = repmat([8 0], T, 1);   % neck
tracks(:,9,:,1) = repmat([5 0], T, 1);   % body
tracks(:,10,:,1) = repmat([0 0], T, 1);  % tailBase
tracks(:,11,:,1) = repmat([-2 0], T, 1); % tailMid
tracks(:,13,:,1) = repmat([4 0], T, 1);  % bodyMid

tracks(:,1,:,2) = repmat([10 + offsetPx 0], T, 1);
tracks(:,2,:,2) = repmat([8 + offsetPx 0], T, 1);
tracks(:,9,:,2) = repmat([5 + offsetPx 0], T, 1);
tracks(:,10,:,2) = repmat([0 + offsetPx 0], T, 1);
tracks(:,11,:,2) = repmat([-2 + offsetPx 0], T, 1);
tracks(:,13,:,2) = repmat([4 + offsetPx 0], T, 1);
end
