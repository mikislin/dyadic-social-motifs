function tests = test_cluster_multiscale_chunks_fields
tests = functiontests(localfunctions);
end

function testOutputFields(testCase)
N = 30;
K = 3;

Cluster = struct();
Cluster.labels = repelem((1:K)', 10);
Cluster.NumClusters = K;
Cluster.maxPosterior = rand(N,1);
Cluster.occupancyFrac = accumarray(Cluster.labels, 1, [K 1]) / N;
Cluster.stability = struct('meanARI', 0.5);

Eval = evaluate_multiscale_clustering(Cluster);

verifyTrue(testCase, isfield(Eval, 'medianPosteriorByCluster'));
verifyEqual(testCase, numel(Eval.medianPosteriorByCluster), K);
verifyEqual(testCase, Eval.nAnchors, N);
end