function [Embedding, Audit] = fit_condition_blind_embedding_pca(Input, dimensionAudit, params)
%FIT_CONDITION_BLIND_EMBEDDING_PCA Fit run-07 per-scale and global PCA.
%
% The fitter consumes only numeric summary matrices and scale guidance from
% run_06 embedding_dimension_audit.csv. Provenance labels remain in metadata
% and are used only after fitting for arena sensitivity diagnostics.

rng(params.rng_seed);

nScale = numel(Input.scale);
ScaleModel = repmat(struct( ...
    'scale_index', [], ...
    'chunk_sec', [], ...
    'hierarchical_role', "", ...
    'n_input_rows', [], ...
    'n_input_dims', [], ...
    'n_pcs_selected', [], ...
    'coeff', [], ...
    'explained', [], ...
    'mu', [], ...
    'score', [], ...
    'score_for_global', [], ...
    'score_center', [], ...
    'score_scale', [], ...
    'preprocess', struct(), ...
    'dimMeta', table(), ...
    'rowMeta', table()), nScale, 1);

pcaByScale = table();
scaleWeights = table();
stabilityAudit = table();
arenaAudit = table();
globalCells = cell(nScale, 1);
rowCells = cell(nScale, 1);
maxPcs = 0;

for s = 1:nScale
    S = Input.scale(s);
    guidance = local_dimension_guidance(dimensionAudit, S.scale_index);
    role = string(S.hierarchical_role);
    if role == ""
        role = string(guidance.initial_band(1));
    end
    targetPcs = double(guidance.recommended_pcs_for_next_embedding(1));
    capPcs = local_role_pc_cap(role, params);

    [Xproc, prep] = local_preprocess_matrix(S.X, params);
    nPC = min([ceil(targetPcs), capPcs, size(Xproc, 1) - 1, size(Xproc, 2)]);
    nPC = max(1, floor(nPC));
    if nPC < params.min_pcs_per_scale && size(Xproc, 1) > params.min_pcs_per_scale && size(Xproc, 2) > params.min_pcs_per_scale
        nPC = min([params.min_pcs_per_scale, size(Xproc, 1) - 1, size(Xproc, 2)]);
    end
    assert(nPC >= 1, 'fit_condition_blind_embedding_pca:NoPCs', ...
        'Scale %.4g s has no usable PCA dimensions.', S.chunk_sec);

    [coeff, score, latent] = pca(Xproc, 'NumComponents', nPC);
    mu = mean(Xproc, 1);
    explained = local_explained_from_latent(latent, Xproc);
    [scoreForGlobal, scoreCenter, scoreScale] = local_prepare_scores_for_global(score, params);
    scaleWeight = local_scale_weight(size(scoreForGlobal, 2), params);
    scoreForGlobal = scoreForGlobal .* scaleWeight;

    dimKept = S.dimMeta(prep.keepMask(:), :);
    ScaleModel(s).scale_index = S.scale_index;
    ScaleModel(s).chunk_sec = S.chunk_sec;
    ScaleModel(s).hierarchical_role = role;
    ScaleModel(s).n_input_rows = size(S.X, 1);
    ScaleModel(s).n_input_dims = size(S.X, 2);
    ScaleModel(s).n_pcs_selected = nPC;
    ScaleModel(s).coeff = coeff;
    ScaleModel(s).explained = explained(:);
    ScaleModel(s).mu = mu;
    ScaleModel(s).score = score;
    ScaleModel(s).score_for_global = scoreForGlobal;
    ScaleModel(s).score_center = scoreCenter;
    ScaleModel(s).score_scale = scoreScale;
    ScaleModel(s).preprocess = prep;
    ScaleModel(s).dimMeta = dimKept;
    ScaleModel(s).rowMeta = local_sanitize_row_manifest(S.rowMeta, true);

    globalCells{s} = scoreForGlobal;
    rowCells{s} = ScaleModel(s).rowMeta;
    maxPcs = max(maxPcs, size(scoreForGlobal, 2));

    pcaByScale = [pcaByScale; local_pca_audit_row(S, guidance, role, targetPcs, capPcs, nPC, explained, prep)]; %#ok<AGROW>
    scaleWeights = [scaleWeights; local_scale_weight_row(S, role, nPC, scaleWeight, params)]; %#ok<AGROW>
    stabilityAudit = [stabilityAudit; local_split_half_stability(Xproc, S, role, nPC, params)]; %#ok<AGROW>
    arenaAudit = [arenaAudit; local_arena_sensitivity(S, role, score, coeff, dimKept, params)]; %#ok<AGROW>
end

[Xglobal, rowManifest] = local_stack_global_matrix(globalCells, rowCells, maxPcs);
nGlobal = min([params.global_n_pcs, size(Xglobal, 1) - 1, size(Xglobal, 2)]);
nGlobal = max(1, floor(nGlobal));
[globalCoeff, globalScore, globalLatent] = pca(Xglobal, 'NumComponents', nGlobal);
globalMu = mean(Xglobal, 1);
globalExplained = local_explained_from_latent(globalLatent, Xglobal);

scoreTable = local_sanitize_row_manifest(rowManifest, false);
for j = 1:size(globalScore, 2)
    scoreTable.(sprintf('global_pc%02d', j)) = globalScore(:, j);
end

Embedding = struct();
Embedding.scale = ScaleModel;
Embedding.global_matrix = Xglobal;
Embedding.global_score = globalScore;
Embedding.global_coeff = globalCoeff;
Embedding.global_mu = globalMu;
Embedding.global_explained = globalExplained(:);
Embedding.global_score_names = compose('global_pc%02d', 1:size(globalScore, 2))';
Embedding.rowManifest = rowManifest;
Embedding.scoreTable = scoreTable;
Embedding.params = params;
Embedding.labels_used_for_fit = "none";
Embedding.arena_used_for_fit = false;
Embedding.condition_used_for_fit = false;

globalAudit = table();
globalAudit.matrix_name = "global_ordinal_pc_stack";
globalAudit.n_rows = size(Xglobal, 1);
globalAudit.n_columns = size(Xglobal, 2);
globalAudit.n_global_pcs = size(globalScore, 2);
globalAudit.finite_value_fraction = mean(isfinite(Xglobal), 'all');
globalAudit.global_pc1_explained = globalExplained(1);
globalAudit.global_cum5_explained = sum(globalExplained(1:min(5, numel(globalExplained))));
globalAudit.global_cum_selected_explained = sum(globalExplained);
globalAudit.labels_used_for_global_pca = "none";
globalAudit.arena_used_for_global_pca = false;
globalAudit.condition_used_for_global_pca = false;
globalAudit.global_matrix_mode = string(params.global_matrix_mode);

Audit = struct();
Audit.pcaByScale = pcaByScale;
Audit.scaleWeights = scaleWeights;
Audit.stability = stabilityAudit;
Audit.arenaSensitivity = arenaAudit;
Audit.global = globalAudit;
end

function guidance = local_dimension_guidance(dimensionAudit, scaleIndex)
required = ["scale_index", "recommended_pcs_for_next_embedding"];
missing = setdiff(required, string(dimensionAudit.Properties.VariableNames));
assert(isempty(missing), 'fit_condition_blind_embedding_pca:BadDimensionAudit', ...
    'embedding_dimension_audit missing required columns: %s', strjoin(missing, ', '));
row = find(double(dimensionAudit.scale_index) == double(scaleIndex), 1);
assert(~isempty(row), 'fit_condition_blind_embedding_pca:MissingDimensionGuidance', ...
    'Missing embedding dimension guidance for scale_index=%g.', scaleIndex);
guidance = dimensionAudit(row, :);
end

function cap = local_role_pc_cap(role, params)
switch string(role)
    case "micro"
        cap = params.micro_max_pcs;
    case "motif"
        cap = params.motif_max_pcs;
    otherwise
        cap = params.context_max_pcs;
end
end

function explained = local_explained_from_latent(latent, X)
latent = double(latent(:));
total = sum(var(double(X), 0, 1), 'omitnan');
if ~(isfinite(total) && total > 0)
    explained = zeros(size(latent));
else
    explained = 100 .* latent ./ total;
end
end

function [Xproc, prep] = local_preprocess_matrix(X, params)
X = double(X);
[n, d] = size(X);
Xproc = zeros(n, d);
center = zeros(1, d);
scale = ones(1, d);
winsorLow = nan(1, d);
winsorHigh = nan(1, d);

for j = 1:d
    x = X(:, j);
    ok = isfinite(x);
    if nnz(ok) < 5
        Xproc(:, j) = 0;
        continue
    end
    q = quantile(x(ok), [params.preprocess_winsor_low, params.preprocess_winsor_high]);
    x(ok) = min(max(x(ok), q(1)), q(2));
    med = median(x(ok), 'omitnan');
    sc = iqr(x(ok));
    if ~(isfinite(sc) && sc > params.preprocess_min_robust_scale)
        sc = std(x(ok), 0, 'omitnan');
    end
    if ~(isfinite(sc) && sc > params.preprocess_min_robust_scale)
        Xproc(:, j) = 0;
        center(j) = med;
        scale(j) = 1;
        winsorLow(j) = q(1);
        winsorHigh(j) = q(2);
        continue
    end
    x(~ok) = med;
    Xproc(:, j) = (x - med) ./ sc;
    center(j) = med;
    scale(j) = sc;
    winsorLow(j) = q(1);
    winsorHigh(j) = q(2);
end

keepMask = any(abs(Xproc) > 0, 1);
assert(any(keepMask), 'fit_condition_blind_embedding_pca:NoFiniteColumns', ...
    'No nonconstant finite columns remained after preprocessing.');
Xproc = Xproc(:, keepMask);

prep = struct();
prep.center = center(keepMask);
prep.scale = scale(keepMask);
prep.winsorLow = winsorLow(keepMask);
prep.winsorHigh = winsorHigh(keepMask);
prep.keepMask = keepMask;
prep.n_input_dims = d;
prep.n_kept_dims = nnz(keepMask);
prep.finite_input_fraction = mean(isfinite(X), 'all');
prep.labels_used_for_preprocessing = "none";
end

function [scoreOut, center, scale] = local_prepare_scores_for_global(score, params)
scoreOut = score;
center = zeros(1, size(score, 2));
scale = ones(1, size(score, 2));
if ~logical(params.standardize_scale_scores_before_global)
    return
end
for j = 1:size(scoreOut, 2)
    x = scoreOut(:, j);
    med = median(x, 'omitnan');
    sc = iqr(x);
    if ~(isfinite(sc) && sc > params.preprocess_min_robust_scale)
        sc = std(x, 0, 'omitnan');
    end
    if ~(isfinite(sc) && sc > params.preprocess_min_robust_scale) || ...
            max(abs(x - med), [], 'omitnan') <= params.preprocess_min_robust_scale
        scoreOut(:, j) = 0;
        center(j) = med;
        scale(j) = 1;
        continue
    end
    if ~(isfinite(sc) && sc > 0)
        sc = 1;
    end
    scoreOut(:, j) = (x - med) ./ sc;
    scoreOut(:, j) = min(max(scoreOut(:, j), -params.scale_score_winsor_abs), ...
        params.scale_score_winsor_abs);
    center(j) = med;
    scale(j) = sc;
end
end

function w = local_scale_weight(nPC, params)
switch string(params.scale_weight_mode)
    case "equal_total_weight"
        w = 1 ./ sqrt(max(nPC, 1));
    otherwise
        w = 1;
end
end

function row = local_pca_audit_row(S, guidance, role, targetPcs, capPcs, nPC, explained, prep)
row = table();
row.scale_index = S.scale_index;
row.chunk_sec = S.chunk_sec;
row.hierarchical_role = role;
row.n_input_rows = size(S.X, 1);
row.n_input_dims = size(S.X, 2);
row.n_kept_dims_after_preprocess = prep.n_kept_dims;
row.input_finite_fraction = prep.finite_input_fraction;
row.run06_recommended_pcs_for_next_embedding = targetPcs;
row.run07_role_pc_cap = capPcs;
row.n_pcs_selected = nPC;
row.pc1_explained = explained(1);
row.cum5_explained = sum(explained(1:min(5, numel(explained))));
row.cum_selected_explained = sum(explained);
if ismember('n_pcs_90pct', guidance.Properties.VariableNames)
    row.run06_n_pcs_90pct = double(guidance.n_pcs_90pct(1));
else
    row.run06_n_pcs_90pct = NaN;
end
row.reaches_run06_recommended_pc_count = nPC >= targetPcs;
row.embedding_dimension_rule = "run06_embedding_dimension_audit_recommendation_with_predefined_run07_role_caps";
row.labels_used_for_pca = "none";
row.arena_used_for_pca = false;
row.condition_used_for_pca = false;
end

function row = local_scale_weight_row(S, role, nPC, scaleWeight, params)
row = table();
row.scale_index = S.scale_index;
row.chunk_sec = S.chunk_sec;
row.hierarchical_role = role;
row.n_pcs_weighted = nPC;
row.scale_weight = scaleWeight;
row.scale_weight_mode = string(params.scale_weight_mode);
row.standardized_before_weighting = logical(params.standardize_scale_scores_before_global);
row.scale_score_winsor_abs_before_weighting = params.scale_score_winsor_abs;
row.scale_weight_rule = string(params.scale_weight_rule);
row.labels_used_for_scale_weight = "none";
row.arena_used_for_scale_weight = false;
row.condition_used_for_scale_weight = false;
end

function row = local_split_half_stability(Xproc, S, role, nPC, params)
if params.run_mode == "smoke"
    nSplits = params.smoke_stability_splits;
else
    nSplits = params.stability_splits;
end
n = size(Xproc, 1);
similarity = nan(nSplits, 1);
pc1Corr = nan(nSplits, 1);
for b = 1:nSplits
    ord = randperm(n);
    nA = floor(n / 2);
    idxA = ord(1:nA);
    idxB = ord((nA + 1):end);
    k = min([params.stability_n_pcs_compared, nPC, numel(idxA) - 1, numel(idxB) - 1, size(Xproc, 2)]);
    if k < 1
        continue
    end
    coeffA = pca(Xproc(idxA, :), 'NumComponents', k);
    coeffB = pca(Xproc(idxB, :), 'NumComponents', k);
    singularVals = svd(coeffA(:, 1:k)' * coeffB(:, 1:k));
    similarity(b) = mean(singularVals, 'omitnan');
    pc1Corr(b) = abs(coeffA(:, 1)' * coeffB(:, 1));
end
row = table();
row.scale_index = S.scale_index;
row.chunk_sec = S.chunk_sec;
row.hierarchical_role = role;
row.n_chunks = n;
row.n_split_repeats = nSplits;
row.n_pcs_compared = min([params.stability_n_pcs_compared, nPC]);
row.median_subspace_similarity = median(similarity, 'omitnan');
row.subspace_similarity_p10 = prctile(similarity, 10);
row.subspace_similarity_p90 = prctile(similarity, 90);
row.median_pc1_abs_correlation = median(pc1Corr, 'omitnan');
row.pc1_abs_correlation_p10 = prctile(pc1Corr, 10);
row.pc1_abs_correlation_p90 = prctile(pc1Corr, 90);
row.bootstrap_unit = "chunks_split_half_within_scale";
row.labels_used_for_stability = "none";
row.arena_used_for_stability = false;
row.condition_used_for_stability = false;
end

function row = local_arena_sensitivity(S, role, score, coeff, dimMeta, params)
k = min([params.arena_sensitivity_n_pcs, size(score, 2)]);
arena = string(S.rowMeta.arena_label);
levels = unique(arena, 'stable');
arenaShift = NaN;
nA = NaN;
nB = NaN;
arenaLevels = "";
if numel(levels) >= 2 && k >= 1
    idxA = arena == levels(1);
    idxB = arena == levels(2);
    nA = nnz(idxA);
    nB = nnz(idxB);
    medA = median(score(idxA, 1:k), 1, 'omitnan');
    medB = median(score(idxB, 1:k), 1, 'omitnan');
    pooled = iqr(score(:, 1:k), 1);
    pooled(~isfinite(pooled) | pooled <= 0) = 1;
    arenaShift = norm((medB - medA) ./ pooled) ./ sqrt(k);
    arenaLevels = strjoin(levels(1:2), ";");
end

[familyText, topFeature, topFeatureFrac, topFeatureText] = local_loading_fractions(coeff, dimMeta, k);

row = table();
row.scale_index = S.scale_index;
row.chunk_sec = S.chunk_sec;
row.hierarchical_role = role;
row.n_embedding_pcs_audited = k;
row.embedding_arena_shift_iqr_units = arenaShift;
row.arena_levels = arenaLevels;
row.n_arena_level_a = nA;
row.n_arena_level_b = nB;
row.loading_fraction_by_feature_family = familyText;
row.top_base_feature_by_loading = topFeature;
row.top_base_feature_loading_fraction = topFeatureFrac;
row.top_feature_loading_fractions = topFeatureText;
row.audit_only_not_selection = true;
row.labels_used_for_embedding_fit = "none";
row.arena_used_for_embedding_fit = false;
row.arena_used_for_posthoc_audit = true;
row.condition_used_for_embedding_fit = false;
end

function [familyText, topFeature, topFeatureFrac, topFeatureText] = local_loading_fractions(coeff, dimMeta, k)
if k < 1 || isempty(coeff)
    familyText = "";
    topFeature = "";
    topFeatureFrac = NaN;
    topFeatureText = "";
    return
end
energy = sum(coeff(:, 1:k).^2, 2);
total = sum(energy);
if ~(isfinite(total) && total > 0)
    total = 1;
end
families = string(dimMeta.feature_family);
uFam = unique(families, 'stable');
parts = strings(numel(uFam), 1);
for i = 1:numel(uFam)
    frac = sum(energy(families == uFam(i))) ./ total;
    parts(i) = uFam(i) + ":" + string(sprintf('%.4f', frac));
end
familyText = strjoin(parts, ";");

features = string(dimMeta.base_feature);
uFeat = unique(features, 'stable');
featFrac = zeros(numel(uFeat), 1);
for i = 1:numel(uFeat)
    featFrac(i) = sum(energy(features == uFeat(i))) ./ total;
end
[featFrac, ord] = sort(featFrac, 'descend');
uFeat = uFeat(ord);
topFeature = uFeat(1);
topFeatureFrac = featFrac(1);
nShow = min(8, numel(uFeat));
topParts = strings(nShow, 1);
for i = 1:nShow
    topParts(i) = uFeat(i) + ":" + string(sprintf('%.4f', featFrac(i)));
end
topFeatureText = strjoin(topParts, ";");
end

function [Xglobal, rowManifest] = local_stack_global_matrix(globalCells, rowCells, maxPcs)
Xglobal = zeros(0, maxPcs);
rowManifest = table();
for s = 1:numel(globalCells)
    X = globalCells{s};
    if isempty(X)
        continue
    end
    P = zeros(size(X, 1), maxPcs);
    P(:, 1:size(X, 2)) = X;
    Xglobal = [Xglobal; P]; %#ok<AGROW>
    rowManifest = [rowManifest; rowCells{s}]; %#ok<AGROW>
end
[rowManifest, ord] = sortrows(rowManifest, 'embedding_row_id');
Xglobal = Xglobal(ord, :);
end

function T = local_sanitize_row_manifest(T, keepArenaAudit)
names = string(T.Properties.VariableNames);
lowerNames = lower(names);
drop = false(size(names));
blockedTokens = ["condition", "cohort", "group", "drug", "genotype", "outcome"];
for i = 1:numel(blockedTokens)
    drop = drop | contains(lowerNames, blockedTokens(i));
end

auditFlag = startsWith(lowerNames, "labels_used") | contains(lowerNames, "_used_for_");
drop(auditFlag) = false;
if ~logical(keepArenaAudit)
    drop = drop | lowerNames == "arena" | lowerNames == "arena_label";
end
T(:, names(drop)) = [];
end
