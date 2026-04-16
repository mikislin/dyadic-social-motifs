function tests = test_compute_dyad_features
tests = functiontests(localfunctions);
end

function testBasicFeatureComputation(testCase)
T = 40;
N = 13;
tracks = nan(T, N, 2, 2);

% Animal 1 moves right.
tracks(:,1,:,1) = [linspace(0,39,T)' zeros(T,1)];  % nose
tracks(:,9,:,1) = [linspace(-5,34,T)' zeros(T,1)];
tracks(:,10,:,1) = [linspace(-10,29,T)' zeros(T,1)];
tracks(:,13,:,1) = [linspace(-7,32,T)' zeros(T,1)];

% Animal 2 is ahead, stationary.
tracks(:,1,:,2) = [ones(T,1)*50 zeros(T,1)];
tracks(:,9,:,2) = [ones(T,1)*45 zeros(T,1)];
tracks(:,10,:,2) = [ones(T,1)*40 zeros(T,1)];
tracks(:,13,:,2) = [ones(T,1)*43 zeros(T,1)];

nodeMap = struct('nose',1,'body',9,'tailBase',10,'midBody',13);
D = compute_dyad_features(tracks, 80, nodeMap);

verifyEqual(testCase, size(D.X,1), T);
verifyTrue(testCase, any(strcmp(D.featureNames, 'centroid_dist')));
verifyTrue(testCase, all(D.frameMask));
end
