%DEMO_RUN08_CONDITION_BLIND_GRAPH_PRINCIPLES Synthetic run-08 principles demo.
%
% This demo is illustrative only. It builds a condition-blind kNN graph from
% numeric synthetic embedding scores, then audits rare/social event flags only
% after graph construction.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
cd(repoRoot);
addpath(genpath(repoRoot));

outRoot = fullfile(repoRoot, 'derived', 'demo_run08_condition_blind_graph_principles');
figDir = fullfile(outRoot, 'figures');
if ~exist(outRoot, 'dir')
    mkdir(outRoot);
end
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

rng(808);
nPerScale = 220;
scaleIndex = repelem((1:3)', nPerScale, 1);
chunkSecLookup = [0.2, 2.1, 13.2];
chunkSec = chunkSecLookup(scaleIndex)';
n = numel(scaleIndex);

theta = rand(n, 1) .* 2 .* pi + 0.42 .* (scaleIndex - 2);
radius = 1 + 0.22 .* randn(n, 1) + 0.08 .* scaleIndex;
contactAxis = 1.15 .* cos(theta) + 0.28 .* randn(n, 1);
approachAxis = 1.05 .* sin(theta) + 0.25 .* randn(n, 1);
withdrawAxis = -0.75 .* sin(theta + 0.7) + 0.32 .* randn(n, 1);
roleAxis = 0.85 .* cos(2 .* theta) + 0.32 .* randn(n, 1);
speedAxis = 0.55 .* radius + 0.40 .* randn(n, 1) + 0.10 .* scaleIndex;
turnAxis = abs(sin(3 .* theta)) + 0.25 .* randn(n, 1);

Xraw = [radius .* cos(theta), radius .* sin(theta), contactAxis, approachAxis, ...
    withdrawAxis, roleAxis, speedAxis, turnAxis];
X = zeros(size(Xraw));
for j = 1:size(Xraw, 2)
    med = median(Xraw(:, j), 'omitnan');
    sc = iqr(Xraw(:, j));
    if ~(isfinite(sc) && sc > 1e-10)
        sc = std(Xraw(:, j), 0, 'omitnan');
    end
    X(:, j) = (Xraw(:, j) - med) ./ sc;
end

rowMeta = table();
rowMeta.graph_node_id = (1:n)';
rowMeta.embedding_row_id = (1:n)';
rowMeta.scale_index = scaleIndex;
rowMeta.chunk_sec = chunkSec;
rowMeta.labels_used_for_demo_graph = repmat("none", n, 1);
rowMeta.arena_used_for_demo_graph = false(n, 1);
rowMeta.condition_used_for_demo_graph = false(n, 1);

params = struct();
params.k_neighbors = 18;
params.knn_block_size = 128;
Graph = build_condition_blind_knn_graph(X, rowMeta, params);

highSpeedEvent = speedAxis >= quantile(speedAxis, 0.92);
contactTransitionEvent = contactAxis > 0.85 & abs(approachAxis) < 0.35;
roleAsymmetryEvent = abs(roleAxis) >= quantile(abs(roleAxis), 0.90);
rareSocialEvent = highSpeedEvent | contactTransitionEvent | roleAsymmetryEvent;

targetRare = double(rareSocialEvent(Graph.Edges.target_node_id));
neighborRareFraction = accumarray(Graph.Edges.source_node_id, targetRare, [n 1], @mean, NaN);
sourceRareFraction = accumarray(Graph.Edges.source_node_id, double(rareSocialEvent(Graph.Edges.source_node_id)), [n 1], @mean, NaN);

nodeAudit = rowMeta;
nodeAudit.demo_axis_1 = Xraw(:, 1);
nodeAudit.demo_axis_2 = Xraw(:, 2);
nodeAudit.high_speed_event_postfit = highSpeedEvent;
nodeAudit.contact_transition_event_postfit = contactTransitionEvent;
nodeAudit.role_asymmetry_event_postfit = roleAsymmetryEvent;
nodeAudit.any_rare_social_event_postfit = rareSocialEvent;
nodeAudit.neighbor_rare_event_fraction = neighborRareFraction;
nodeAudit.labels_used_for_demo_graph = repmat("none", n, 1);
nodeAudit.event_flags_used_for_demo_graph = false(n, 1);

scales = unique(scaleIndex, 'stable')';
coverage = table();
for s = scales
    idx = scaleIndex == s;
    one = table();
    one.scale_index = double(s);
    one.chunk_sec = chunkSecLookup(s);
    one.n_nodes = nnz(idx);
    one.high_speed_event_fraction = mean(highSpeedEvent(idx));
    one.contact_transition_event_fraction = mean(contactTransitionEvent(idx));
    one.role_asymmetry_event_fraction = mean(roleAsymmetryEvent(idx));
    one.any_rare_social_event_fraction = mean(rareSocialEvent(idx));
    one.mean_neighbor_rare_event_fraction = mean(neighborRareFraction(idx), 'omitnan');
    one.event_flags_used_for_graph = false;
    one.labels_used_for_graph = "none";
    one.condition_used_for_graph = false;
    coverage = [coverage; one]; %#ok<AGROW>
end

summary = table();
summary.n_nodes = n;
summary.n_numeric_graph_dimensions = size(X, 2);
summary.k_neighbors = Graph.k;
summary.n_directed_edges = height(Graph.Edges);
summary.overall_rare_social_event_fraction = mean(rareSocialEvent);
summary.mean_neighbor_rare_event_fraction = mean(neighborRareFraction, 'omitnan');
summary.mean_source_rare_event_fraction = mean(sourceRareFraction, 'omitnan');
summary.graph_inputs = "synthetic numeric embedding scores only";
summary.postfit_audits = "synthetic rare/social event flags and scale labels";
summary.labels_used_for_graph = "none";
summary.event_flags_used_for_graph = false;

writetable(nodeAudit, fullfile(outRoot, 'demo_run08_node_event_audit.csv'));
writetable(Graph.Edges, fullfile(outRoot, 'demo_run08_edge_list.csv'));
writetable(coverage, fullfile(outRoot, 'demo_run08_event_coverage_by_scale.csv'));
writetable(summary, fullfile(outRoot, 'demo_run08_principles_summary.csv'));

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1500 760]);
tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
scatter(Xraw(:, 1), Xraw(:, 2), 14, chunkSec, 'filled');
xlabel('Synthetic embedding axis 1', 'Interpreter', 'none');
ylabel('Synthetic embedding axis 2', 'Interpreter', 'none');
title('Graph inputs are numeric only', 'Interpreter', 'none');
colorbar;
box off; grid on;

nexttile;
scatter(Xraw(:, 1), Xraw(:, 2), 14, rareSocialEvent, 'filled');
xlabel('Synthetic embedding axis 1', 'Interpreter', 'none');
ylabel('Synthetic embedding axis 2', 'Interpreter', 'none');
title('Rare/social events are post-fit audits', 'Interpreter', 'none');
colormap(gca, [0.62 0.62 0.62; 0.80 0.20 0.20]);
colorbar;
box off; grid on;

nexttile;
plot(coverage.chunk_sec, coverage.any_rare_social_event_fraction, 'o-', ...
    'Color', [0.80 0.20 0.20], 'MarkerFaceColor', [0.80 0.20 0.20], 'LineWidth', 1.5);
hold on;
plot(coverage.chunk_sec, coverage.mean_neighbor_rare_event_fraction, 's-', ...
    'Color', [0.00 0.62 0.45], 'MarkerFaceColor', [0.00 0.62 0.45], 'LineWidth', 1.5);
set(gca, 'XScale', 'log');
ylim([0 1]);
xlabel('Synthetic chunk scale (s)', 'Interpreter', 'none');
ylabel('Fraction', 'Interpreter', 'none');
title('Coverage is audited after graph build', 'Interpreter', 'none');
legend({'Node event fraction','Neighbor event fraction'}, 'Location', 'best', 'Box', 'off');
box off; grid on;

exportgraphics(fig, fullfile(figDir, 'demo_run08_condition_blind_graph_principles.png'), ...
    'Resolution', 300);
exportgraphics(fig, fullfile(figDir, 'demo_run08_condition_blind_graph_principles.pdf'), ...
    'ContentType', 'vector');
close(fig);

fprintf('Demo outputs written to %s\n', outRoot);
fprintf('Synthetic graph nodes: %d | k: %d | directed edges: %d\n', ...
    n, Graph.k, height(Graph.Edges));
