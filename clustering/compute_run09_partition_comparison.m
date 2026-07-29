function Metrics = compute_run09_partition_comparison(a, b)
%COMPUTE_RUN09_PARTITION_COMPARISON Exact ARI and VI for two partitions.

a = double(canonicalize_run09_membership(a));
b = double(canonicalize_run09_membership(b));
assert(numel(a) == numel(b) && ~isempty(a), ...
    'compute_run09_partition_comparison:DimensionMismatch', ...
    'Partitions must be nonempty and have equal node counts.');
n = numel(a);
joint = sparse(a, b, 1, max(a), max(b));
jointCounts = nonzeros(joint);
rowCounts = full(sum(joint, 2));
columnCounts = full(sum(joint, 1))';
comb2 = @(x) x .* (x - 1) / 2;
index = sum(comb2(jointCounts));
rowPairs = sum(comb2(rowCounts));
columnPairs = sum(comb2(columnCounts));
totalPairs = comb2(n);
expected = rowPairs * columnPairs / max(totalPairs, eps);
maximum = 0.5 * (rowPairs + columnPairs);
denominator = maximum - expected;
if abs(denominator) <= eps(maximum)
    ari = double(isequal(a, b));
else
    ari = (index - expected) / denominator;
end

[jointRow, jointColumn, counts] = find(joint);
pJoint = counts / n;
pRow = rowCounts / n;
pColumn = columnCounts / n;
hA = -sum(pRow(pRow > 0) .* log(pRow(pRow > 0)));
hB = -sum(pColumn(pColumn > 0) .* log(pColumn(pColumn > 0)));
mutualInformation = sum(pJoint .* log(pJoint ./ ...
    (pRow(jointRow) .* pColumn(jointColumn))));
vi = max(0, hA + hB - 2 * mutualInformation);

Metrics = struct();
Metrics.adjusted_rand_index = ari;
Metrics.variation_of_information = vi;
Metrics.normalized_variation_of_information = ...
    vi / max(log(n), eps);
end
