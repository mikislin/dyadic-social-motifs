function Eval = evaluate_multiscale_clustering(Cluster)
%EVALUATE_MULTISCALE_CLUSTERING Evaluate clustering output from cluster_multiscale_chunks.

assert(isstruct(Cluster), 'Cluster must be a struct.');
assert(isfield(Cluster, 'labels'), 'Cluster.labels is required.');
assert(isfield(Cluster, 'NumClusters'), 'Cluster.NumClusters is required.');
assert(isfield(Cluster, 'maxPosterior'), 'Cluster.maxPosterior is required.');
assert(isfield(Cluster, 'occupancyFrac'), 'Cluster.occupancyFrac is required.');

labels = Cluster.labels(:);
K = Cluster.NumClusters;

confByCluster = accumarray(labels, Cluster.maxPosterior(:), [K 1], @median, NaN);
countByCluster = accumarray(labels, 1, [K 1], @sum, 0);

runBreaks = [1; find(diff(labels) ~= 0) + 1; numel(labels) + 1];
runLengths = diff(runBreaks);

Eval = struct();
Eval.NumClusters = K;
Eval.nAnchors = numel(labels);
Eval.clusterCounts = countByCluster;
Eval.clusterFrac = Cluster.occupancyFrac(:);
Eval.medianPosteriorByCluster = confByCluster;
Eval.medianMaxPosterior = median(Cluster.maxPosterior, 'omitnan');
Eval.runLengthQuantiles = prctile(runLengths, [10 25 50 75 90]);
Eval.meanARI = NaN;
if isfield(Cluster, 'stability') && isfield(Cluster.stability, 'meanARI')
    Eval.meanARI = Cluster.stability.meanARI;
end
end
