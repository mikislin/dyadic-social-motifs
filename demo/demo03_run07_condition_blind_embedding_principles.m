%DEMO_RUN07_CONDITION_BLIND_EMBEDDING_PRINCIPLES
% Illustrates the key design principles behind run_07 on synthetic data.
%
% This demo is intentionally separate from the paper pipeline. It does not
% read real condition labels or write run_07 artifacts. It shows why the
% embedding step uses scale-specific dimensions, label-free scale weights,
% post-fit arena/event audits, and a separate rare-event coverage audit.

clear;
clc;

thisFile = mfilename('fullpath');
if isempty(thisFile)
    repoRoot = pwd;
else
    repoRoot = fileparts(fileparts(thisFile));
end
cd(repoRoot);
addpath(genpath(repoRoot));

rng(707);
outRoot = fullfile(repoRoot, 'derived', 'demo_run07_condition_blind_embedding_principles');
figDir = fullfile(outRoot, 'figures');
if ~exist(outRoot, 'dir')
    mkdir(outRoot);
end
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

scales = [0.25; 0.75; 2.5; 8; 20];
roles = ["micro"; "micro"; "motif"; "context"; "context"];
truthDims = [14; 18; 32; 42; 20];
roleCaps = [24; 24; 48; 64; 64];
nPerScale = 240;
nFeatures = 34;
nSummaryDims = 6 * nFeatures;
nRareRequested = 18;
globalTargetPCs = 12;
fixedPCs = 5;

scoreCells = cell(numel(scales), 1);
rowCells = cell(numel(scales), 1);
pcaAudit = table();
coverageAudit = table();
scaleWeightAudit = table();

fprintf('demo_run07_condition_blind_embedding_principles\n');
fprintf('Output root: %s\n', outRoot);

for s = 1:numel(scales)
    n = nPerScale;
    d = nSummaryDims;
    kTrue = truthDims(s);

    sessionIndex = (1:n)';
    anchorTime = linspace(5, 900, n)' + 0.01 .* randn(n, 1);
    arenaLabel = repmat(["big"; "small"], ceil(n / 2), 1);
    arenaLabel = arenaLabel(1:n);
    conditionLabel = repmat(["held_out_A"; "held_out_B"], ceil(n / 2), 1);
    conditionLabel = conditionLabel(randperm(n));

    rareEvent = false(n, 1);
    rareEvent(round(linspace(12, n - 8, nRareRequested))) = true;

    loading = randn(d, kTrue);
    loading = loading ./ max(vecnorm(loading, 2, 1), eps);
    latent = randn(n, kTrue);
    slowContext = sin(anchorTime ./ (45 + 3 * scales(s)));
    latent(:, 1) = latent(:, 1) + 1.5 .* slowContext;
    latent(:, 2) = latent(:, 2) + double(arenaLabel == "small") .* 0.35;
    latent(rareEvent, min(3, kTrue)) = latent(rareEvent, min(3, kTrue)) + 4.0;

    X = latent * loading' + 0.35 .* randn(n, d);
    rareFeatureBlock = 1:min(18, d);
    X(rareEvent, rareFeatureBlock) = X(rareEvent, rareFeatureBlock) + 2.5;

    Xproc = zeros(size(X));
    keep = false(1, d);
    for j = 1:d
        x = X(:, j);
        med = median(x, 'omitnan');
        sc = iqr(x);
        if ~(isfinite(sc) && sc > 1e-8)
            sc = std(x, 0, 'omitnan');
        end
        if ~(isfinite(sc) && sc > 1e-8)
            Xproc(:, j) = 0;
        else
            Xproc(:, j) = (x - med) ./ sc;
            keep(j) = any(abs(Xproc(:, j)) > 0);
        end
    end
    Xproc = Xproc(:, keep);

    [~, scoreAll, latentVar] = pca(Xproc);
    explained = 100 .* latentVar ./ sum(latentVar, 'omitnan');
    cexpl = cumsum(explained);
    targetPCs = find(cexpl >= 85, 1, 'first');
    if isempty(targetPCs)
        targetPCs = min(size(scoreAll, 2), roleCaps(s));
    end
    selectedPCs = min([targetPCs, roleCaps(s), size(scoreAll, 2), n - 1]);
    selectedPCs = max(1, selectedPCs);
    fixedCum = sum(explained(1:min(fixedPCs, numel(explained))));
    selectedCum = sum(explained(1:selectedPCs));

    score = scoreAll(:, 1:selectedPCs);
    for j = 1:size(score, 2)
        med = median(score(:, j), 'omitnan');
        sc = iqr(score(:, j));
        if ~(isfinite(sc) && sc > 1e-8)
            sc = std(score(:, j), 0, 'omitnan');
        end
        if ~(isfinite(sc) && sc > 1e-8)
            score(:, j) = 0;
        else
            score(:, j) = (score(:, j) - med) ./ sc;
            score(:, j) = min(max(score(:, j), -8), 8);
        end
    end
    scaleWeight = 1 ./ sqrt(size(score, 2));
    score = score .* scaleWeight;

    rows = table();
    rows.demo_row_id = ((s - 1) * n + (1:n))';
    rows.scale_index = repmat(s, n, 1);
    rows.chunk_sec = repmat(scales(s), n, 1);
    rows.hierarchical_role = repmat(roles(s), n, 1);
    rows.session_index = sessionIndex;
    rows.anchor_time_s = anchorTime;
    rows.rare_event_audit_only = rareEvent;
    rows.arena_label_audit_only = arenaLabel;
    rows.labels_used_for_fit = repmat("none", n, 1);
    rows.arena_used_for_fit = false(n, 1);
    rows.condition_used_for_fit = false(n, 1);
    scoreCells{s} = score;
    rowCells{s} = rows;

    one = table();
    one.scale_index = s;
    one.chunk_sec = scales(s);
    one.hierarchical_role = roles(s);
    one.synthetic_true_latent_dims = kTrue;
    one.fixed_5_pc_cum_explained = fixedCum;
    one.scale_specific_target_pcs_85pct = targetPCs;
    one.role_pc_cap = roleCaps(s);
    one.selected_pcs = selectedPCs;
    one.selected_pc_cum_explained = selectedCum;
    one.labels_used_for_pca = "none";
    one.arena_used_for_pca = false;
    one.condition_used_for_pca = false;
    pcaAudit = [pcaAudit; one]; %#ok<AGROW>

    two = table();
    two.scale_index = s;
    two.chunk_sec = scales(s);
    two.n_uniform_anchors = n;
    two.n_rare_events_available = nnz(rareEvent);
    two.rare_event_anchor_fraction = mean(rareEvent);
    two.condition_labels_used_for_rare_event_coverage = false;
    two.recommended_audit = "track rare event coverage separately; augment anchors only with condition-blind event definitions";
    coverageAudit = [coverageAudit; two]; %#ok<AGROW>

    three = table();
    three.scale_index = s;
    three.chunk_sec = scales(s);
    three.selected_pcs = selectedPCs;
    three.scale_weight = scaleWeight;
    three.scale_weight_rule = "equal_total_weight";
    three.labels_used_for_scale_weight = "none";
    three.arena_used_for_scale_weight = false;
    three.condition_used_for_scale_weight = false;
    scaleWeightAudit = [scaleWeightAudit; three]; %#ok<AGROW>
end

maxPC = max(cellfun(@(x) size(x, 2), scoreCells));
Xglobal = zeros(numel(scales) * nPerScale, maxPC);
rowManifest = table();
cursor = 1;
for s = 1:numel(scales)
    score = scoreCells{s};
    nr = size(score, 1);
    Xglobal(cursor:(cursor + nr - 1), 1:size(score, 2)) = score;
    rowManifest = [rowManifest; rowCells{s}]; %#ok<AGROW>
    cursor = cursor + nr;
end

[~, globalScore, globalLatent] = pca(Xglobal, 'NumComponents', globalTargetPCs);
globalExplained = 100 .* globalLatent ./ sum(var(Xglobal, 0, 1), 'omitnan');

globalScores = rowManifest(:, {'demo_row_id', 'scale_index', 'chunk_sec', ...
    'hierarchical_role', 'session_index', 'anchor_time_s', ...
    'labels_used_for_fit', 'arena_used_for_fit', 'condition_used_for_fit'});
for j = 1:size(globalScore, 2)
    globalScores.(sprintf('global_pc%02d', j)) = globalScore(:, j);
end

arenaAudit = rowManifest(:, {'demo_row_id', 'scale_index', 'chunk_sec', ...
    'arena_label_audit_only', 'rare_event_audit_only'});
arenaAudit.audit_only_not_used_for_fit = true(height(arenaAudit), 1);

globalAudit = table();
globalAudit.n_rows = size(Xglobal, 1);
globalAudit.n_columns = size(Xglobal, 2);
globalAudit.n_global_pcs = size(globalScore, 2);
globalAudit.global_pc1_explained = globalExplained(1);
globalAudit.global_cum5_explained = sum(globalExplained(1:min(5, numel(globalExplained))));
globalAudit.global_cum_selected_explained = sum(globalExplained);
globalAudit.labels_used_for_global_pca = "none";
globalAudit.arena_used_for_global_pca = false;
globalAudit.condition_used_for_global_pca = false;

writetable(pcaAudit, fullfile(outRoot, 'demo_pca_dimension_audit.csv'));
writetable(coverageAudit, fullfile(outRoot, 'demo_rare_event_coverage_audit.csv'));
writetable(scaleWeightAudit, fullfile(outRoot, 'demo_scale_weight_audit.csv'));
writetable(globalScores, fullfile(outRoot, 'demo_global_scores_label_free.csv'));
writetable(arenaAudit, fullfile(outRoot, 'demo_postfit_arena_and_event_audit.csv'));
writetable(globalAudit, fullfile(outRoot, 'demo_global_embedding_audit.csv'));

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1500 980]);
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
plot(pcaAudit.chunk_sec, pcaAudit.fixed_5_pc_cum_explained, 'o-', 'LineWidth', 1.5);
plot(pcaAudit.chunk_sec, pcaAudit.selected_pc_cum_explained, 's-', 'LineWidth', 1.5);
set(gca, 'XScale', 'log');
xlabel('Synthetic chunk scale (s)', 'Interpreter', 'none');
ylabel('Cumulative variance (%)', 'Interpreter', 'none');
title('Fixed 5 PCs can underdescribe some scales', 'Interpreter', 'none');
legend({'Fixed 5 PCs', 'Scale-specific PCs'}, 'Location', 'best', 'Box', 'off');
grid on;

nexttile;
hold on;
plot(pcaAudit.chunk_sec, pcaAudit.scale_specific_target_pcs_85pct, 'o-', 'LineWidth', 1.5);
plot(pcaAudit.chunk_sec, pcaAudit.selected_pcs, 's-', 'LineWidth', 1.5);
set(gca, 'XScale', 'log');
xlabel('Synthetic chunk scale (s)', 'Interpreter', 'none');
ylabel('PC count', 'Interpreter', 'none');
title('Dimensionality is scale-specific and capped', 'Interpreter', 'none');
legend({'85% target', 'Selected after cap'}, 'Location', 'best', 'Box', 'off');
grid on;

nexttile;
scatter(globalScores.global_pc01, globalScores.global_pc02, 14, globalScores.chunk_sec, 'filled');
xlabel('Global PC1', 'Interpreter', 'none');
ylabel('Global PC2', 'Interpreter', 'none');
title('Label-free global embedding colored by scale', 'Interpreter', 'none');
colorbar;
grid on;

nexttile;
hold on;
idxRare = arenaAudit.rare_event_audit_only;
scatter(globalScores.global_pc01(~idxRare), globalScores.global_pc02(~idxRare), ...
    10, [0.45 0.45 0.45], 'filled');
scatter(globalScores.global_pc01(idxRare), globalScores.global_pc02(idxRare), ...
    22, [0.80 0.20 0.20], 'filled');
xlabel('Global PC1', 'Interpreter', 'none');
ylabel('Global PC2', 'Interpreter', 'none');
title('Rare events are audited after fitting', 'Interpreter', 'none');
legend({'background anchors', 'rare event audit'}, 'Location', 'best', 'Box', 'off');
grid on;

exportgraphics(fig, fullfile(figDir, 'demo_run07_embedding_principles.png'), 'Resolution', 320);
exportgraphics(fig, fullfile(figDir, 'demo_run07_embedding_principles.pdf'), 'ContentType', 'vector');
close(fig);

fprintf('Wrote demo audits and figure to: %s\n', outRoot);
