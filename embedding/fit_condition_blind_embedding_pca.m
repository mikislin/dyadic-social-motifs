function [Embedding, Audit] = fit_condition_blind_embedding_pca(Input, dimensionAudit, params)
%FIT_CONDITION_BLIND_EMBEDDING_PCA Fit run-07 per-scale and global PCA.
%
% The fitter consumes only numeric summary matrices and scale guidance from
% run_06 embedding_dimension_audit.csv. Provenance labels remain in metadata
% and are used only after fitting for arena sensitivity diagnostics.

rng(params.rng_seed);
assert(~isfield(params, 'use_anchor_weights_for_pca') || ~logical(params.use_anchor_weights_for_pca), ...
    'fit_condition_blind_embedding_pca:WeightedPcaNotEnabled', ...
    'Initial enriched embedding must remain an unweighted PCA fit.');

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
preprocessDimensionAudit = table();
anchorStageSensitivityAudit = table();
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

    [Xproc, prep, prepDimension] = local_preprocess_matrix(S.X, params);
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

    pcaByScale = [pcaByScale; local_pca_audit_row(S, guidance, role, targetPcs, capPcs, nPC, explained, prep, coeff, dimKept)]; %#ok<AGROW>
    scaleWeights = [scaleWeights; local_scale_weight_row(S, role, nPC, scaleWeight, params)]; %#ok<AGROW>
    stabilityAudit = [stabilityAudit; local_split_half_stability(Xproc, S, role, nPC, params)]; %#ok<AGROW>
    arenaAudit = [arenaAudit; local_arena_sensitivity(S, role, score, coeff, dimKept, params)]; %#ok<AGROW>
    preprocessDimensionAudit = [preprocessDimensionAudit; ...
        local_preprocess_dimension_audit(S, prepDimension, params)]; %#ok<AGROW>
    anchorStageSensitivityAudit = [anchorStageSensitivityAudit; ...
        local_anchor_stage_pca_sensitivity(Xproc, S, role, nPC, coeff, explained, params)]; %#ok<AGROW>
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
globalAudit.global_cum_selected_explained = ...
    sum(globalExplained(1:min(nGlobal, numel(globalExplained))));
globalAudit.labels_used_for_global_pca = "none";
globalAudit.arena_used_for_global_pca = false;
globalAudit.condition_used_for_global_pca = false;
globalAudit.global_matrix_mode = string(params.global_matrix_mode);
globalAudit.weights_used_for_global_pca = false;
if ismember('anchor_manifest_mode', rowManifest.Properties.VariableNames)
    globalAudit.anchor_manifest_mode = strjoin(unique(string(rowManifest.anchor_manifest_mode), 'stable'), ";");
else
    globalAudit.anchor_manifest_mode = "primary";
end

Audit = struct();
Audit.pcaByScale = pcaByScale;
Audit.scaleWeights = scaleWeights;
Audit.stability = stabilityAudit;
Audit.arenaSensitivity = arenaAudit;
Audit.preprocessDimensions = preprocessDimensionAudit;
Audit.anchorStageSensitivity = anchorStageSensitivityAudit;
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

function [Xproc, prep, diagnostic] = local_preprocess_matrix(X, params)
X = double(X);
[n, d] = size(X);
Xproc = zeros(n, d);
center = zeros(1, d);
scale = ones(1, d);
winsorLow = nan(1, d);
winsorHigh = nan(1, d);
finiteCount = zeros(1, d);
winsorIqr = nan(1, d);
winsorStd = nan(1, d);
iqrToStdRatio = nan(1, d);
scaleMethod = repmat("dropped_insufficient_finite", 1, d);
safeguardTriggered = false(1, d);
standardizedStd = nan(1, d);
standardizedAbsP99 = nan(1, d);
standardizedAbsMax = nan(1, d);
tailCount = zeros(1, d);
severeTailCount = zeros(1, d);

for j = 1:d
    x = X(:, j);
    ok = isfinite(x);
    finiteCount(j) = nnz(ok);
    if nnz(ok) < 5
        Xproc(:, j) = 0;
        continue
    end
    q = quantile(x(ok), [params.preprocess_winsor_low, params.preprocess_winsor_high]);
    x(ok) = min(max(x(ok), q(1)), q(2));
    med = median(x(ok), 'omitnan');
    robustSc = iqr(x(ok));
    stdSc = std(x(ok), 0, 'omitnan');
    winsorIqr(j) = robustSc;
    winsorStd(j) = stdSc;
    if isfinite(stdSc) && stdSc > params.preprocess_min_robust_scale
        iqrToStdRatio(j) = robustSc ./ stdSc;
    end

    useSparseGuard = logical(params.preprocess_sparse_scale_safeguard_enabled) && ...
        isfinite(robustSc) && robustSc > params.preprocess_min_robust_scale && ...
        isfinite(stdSc) && stdSc > params.preprocess_min_robust_scale && ...
        isfinite(iqrToStdRatio(j)) && ...
        iqrToStdRatio(j) < params.preprocess_min_iqr_to_std_ratio;

    if useSparseGuard
        sc = stdSc;
        scaleMethod(j) = "winsorized_std_sparse_guard";
        safeguardTriggered(j) = true;
    elseif isfinite(robustSc) && robustSc > params.preprocess_min_robust_scale
        sc = robustSc;
        scaleMethod(j) = "winsorized_iqr";
    else
        sc = stdSc;
        scaleMethod(j) = "winsorized_std_low_iqr";
    end
    if ~(isfinite(sc) && sc > params.preprocess_min_robust_scale)
        Xproc(:, j) = 0;
        center(j) = med;
        scale(j) = 1;
        winsorLow(j) = q(1);
        winsorHigh(j) = q(2);
        scaleMethod(j) = "dropped_no_usable_spread";
        continue
    end
    x(~ok) = med;
    z = (x - med) ./ sc;
    Xproc(:, j) = z;
    center(j) = med;
    scale(j) = sc;
    winsorLow(j) = q(1);
    winsorHigh(j) = q(2);
    standardizedStd(j) = std(z, 0, 'omitnan');
    standardizedAbsP99(j) = prctile(abs(z), 99);
    standardizedAbsMax(j) = max(abs(z), [], 'omitnan');
    tailCount(j) = nnz(abs(z) > params.preprocess_tail_audit_abs_threshold);
    severeTailCount(j) = nnz(abs(z) > params.preprocess_severe_tail_audit_abs_threshold);
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

diagnostic = table();
diagnostic.input_dimension_index = (1:d)';
diagnostic.kept_after_preprocess = keepMask(:);
diagnostic.kept_dimension_index = nan(d, 1);
diagnostic.kept_dimension_index(keepMask(:)) = (1:nnz(keepMask))';
diagnostic.n_rows = repmat(n, d, 1);
diagnostic.n_finite_input = finiteCount(:);
diagnostic.finite_input_fraction = finiteCount(:) ./ max(n, 1);
diagnostic.winsor_low = winsorLow(:);
diagnostic.winsor_high = winsorHigh(:);
diagnostic.winsorized_median = center(:);
diagnostic.winsorized_iqr = winsorIqr(:);
diagnostic.winsorized_std = winsorStd(:);
diagnostic.iqr_to_std_ratio = iqrToStdRatio(:);
diagnostic.selected_scale = scale(:);
diagnostic.scale_method = scaleMethod(:);
diagnostic.sparse_scale_safeguard_triggered = safeguardTriggered(:);
diagnostic.standardized_std = standardizedStd(:);
diagnostic.standardized_abs_p99 = standardizedAbsP99(:);
diagnostic.standardized_abs_max = standardizedAbsMax(:);
diagnostic.n_abs_gt_tail_threshold = tailCount(:);
diagnostic.fraction_abs_gt_tail_threshold = tailCount(:) ./ max(n, 1);
diagnostic.n_abs_gt_severe_tail_threshold = severeTailCount(:);
diagnostic.fraction_abs_gt_severe_tail_threshold = severeTailCount(:) ./ max(n, 1);
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

function row = local_pca_audit_row(S, guidance, role, targetPcs, capPcs, nPC, explained, prep, coeff, dimKept)
row = table();
row.scale_index = S.scale_index;
row.chunk_sec = S.chunk_sec;
row.hierarchical_role = role;
row.n_input_rows = size(S.X, 1);
if ismember('anchor_manifest_mode', S.rowMeta.Properties.VariableNames)
    row.anchor_manifest_mode = strjoin(unique(string(S.rowMeta.anchor_manifest_mode), 'stable'), ";");
else
    row.anchor_manifest_mode = "primary";
end
if ismember('anchor_stage', S.rowMeta.Properties.VariableNames)
    stage = string(S.rowMeta.anchor_stage);
    row.n_base_rows = nnz(stage == "base_time_even");
    row.n_rare_enriched_rows = nnz(stage == "rare_strata_enriched");
else
    row.n_base_rows = size(S.X, 1);
    row.n_rare_enriched_rows = 0;
end
if ismember('audit_inverse_probability_weight', S.rowMeta.Properties.VariableNames)
    row.audit_inverse_probability_effective_sample_size = ...
        local_effective_sample_size(double(S.rowMeta.audit_inverse_probability_weight));
else
    row.audit_inverse_probability_effective_sample_size = size(S.X, 1);
end
row.weights_used_for_pca = false;
row.n_input_dims = size(S.X, 2);
row.n_kept_dims_after_preprocess = prep.n_kept_dims;
row.input_finite_fraction = prep.finite_input_fraction;
row.run06_recommended_pcs_for_next_embedding = targetPcs;
row.run07_role_pc_cap = capPcs;
row.n_pcs_selected = nPC;
row.pc1_explained = explained(1);
row.cum5_explained = sum(explained(1:min(5, numel(explained))));
row.cum_selected_explained = sum(explained(1:min(nPC, numel(explained))));
[topLoading, topFeature, topKind, topDct, topDim] = ...
    local_pc1_loading_dominance(coeff, dimKept);
row.pc1_top_dimension_loading_fraction = topLoading;
row.pc1_top_base_feature = topFeature;
row.pc1_top_summary_kind = topKind;
row.pc1_top_dct_coefficient = topDct;
row.pc1_top_embedding_feature_index = topDim;
if ismember('n_pcs_90pct', guidance.Properties.VariableNames)
    row.run06_n_pcs_90pct = double(guidance.n_pcs_90pct(1));
else
    row.run06_n_pcs_90pct = NaN;
end
if ismember('summary_profile', guidance.Properties.VariableNames)
    row.run06_summary_profile = string(guidance.summary_profile(1));
else
    row.run06_summary_profile = "";
end
if ismember('n_temporal_bins', guidance.Properties.VariableNames)
    row.run06_n_temporal_bins = double(guidance.n_temporal_bins(1));
else
    row.run06_n_temporal_bins = NaN;
end
if ismember('n_dct_coeffs', guidance.Properties.VariableNames)
    row.run06_n_dct_coeffs = double(guidance.n_dct_coeffs(1));
else
    row.run06_n_dct_coeffs = NaN;
end
if ismember('n_summary_dims', guidance.Properties.VariableNames)
    row.run06_summary_dims = double(guidance.n_summary_dims(1));
else
    row.run06_summary_dims = NaN;
end
row.reaches_run06_recommended_pc_count = nPC >= targetPcs;
row.embedding_dimension_rule = "run06_embedding_dimension_audit_recommendation_with_predefined_run07_role_caps";
row.labels_used_for_pca = "none";
row.arena_used_for_pca = false;
row.condition_used_for_pca = false;
end

function [topLoading, topFeature, topKind, topDct, topDim] = ...
        local_pc1_loading_dominance(coeff, dimMeta)
topLoading = NaN;
topFeature = "";
topKind = "";
topDct = NaN;
topDim = NaN;
if isempty(coeff) || isempty(dimMeta)
    return
end
energy = double(coeff(:, 1)).^2;
denom = sum(energy, 'omitnan');
if ~(isfinite(denom) && denom > 0)
    return
end
[topLoading, idx] = max(energy ./ denom);
topFeature = local_table_string_value(dimMeta, 'base_feature', idx);
topKind = local_table_string_value(dimMeta, 'summary_kind', idx);
topDct = local_table_numeric_value(dimMeta, 'dct_coefficient', idx);
topDim = local_table_numeric_value(dimMeta, 'embedding_feature_index', idx);
end

function T = local_preprocess_dimension_audit(S, diagnostic, params)
n = height(diagnostic);
assert(height(S.dimMeta) == n, ...
    'fit_condition_blind_embedding_pca:PreprocessAuditDimensionMismatch', ...
    'Scale %.4g s preprocessing audit has %d dimensions but dimMeta has %d rows.', ...
    S.chunk_sec, n, height(S.dimMeta));

T = table();
T.scale_index = repmat(S.scale_index, n, 1);
T.chunk_sec = repmat(S.chunk_sec, n, 1);
T.hierarchical_role = repmat(string(S.hierarchical_role), n, 1);
T.embedding_feature_index = local_table_numeric_column(S.dimMeta, 'embedding_feature_index', n);
T.summary_dim_index = local_table_numeric_column(S.dimMeta, 'summary_dim_index', n);
T.obs_name = local_table_string_column(S.dimMeta, 'obs_name', n);
T.base_feature = local_table_string_column(S.dimMeta, 'base_feature', n);
T.channel_type = local_table_string_column(S.dimMeta, 'channel_type', n);
T.summary_kind = local_table_string_column(S.dimMeta, 'summary_kind', n);
T.temporal_bin = local_table_numeric_column(S.dimMeta, 'temporal_bin', n);
T.dct_coefficient = local_table_numeric_column(S.dimMeta, 'dct_coefficient', n);
T.feature_family = local_table_string_column(S.dimMeta, 'feature_family', n);
T.unit = local_table_string_column(S.dimMeta, 'unit', n);
T.summary_profile = local_table_string_column(S.dimMeta, 'summary_profile', n);
T = [T, diagnostic];
T.sparse_scale_safeguard_enabled = repmat( ...
    logical(params.preprocess_sparse_scale_safeguard_enabled), n, 1);
T.min_iqr_to_std_ratio = repmat(params.preprocess_min_iqr_to_std_ratio, n, 1);
T.tail_abs_threshold = repmat(params.preprocess_tail_audit_abs_threshold, n, 1);
T.severe_tail_abs_threshold = repmat(params.preprocess_severe_tail_audit_abs_threshold, n, 1);
T.scale_selection_rule = repmat( ...
    "winsorized_iqr_unless_low_absolute_spread_or_condition_blind_iqr_to_std_sparse_guard", n, 1);
T.labels_used_for_preprocessing = repmat("none", n, 1);
T.arena_used_for_preprocessing = false(n, 1);
T.condition_used_for_preprocessing = false(n, 1);
end

function row = local_anchor_stage_pca_sensitivity(Xproc, S, role, nPC, combinedCoeff, combinedExplained, params)
row = table();
row.scale_index = S.scale_index;
row.chunk_sec = S.chunk_sec;
row.hierarchical_role = role;
row.n_combined_rows = size(Xproc, 1);
row.n_base_rows = 0;
row.n_rare_enriched_rows = 0;
row.n_pcs_compared = 0;
row.combined_pc1_explained = combinedExplained(1);
row.combined_cum_compared_explained = NaN;
row.base_pc1_explained = NaN;
row.base_cum_compared_explained = NaN;
row.rare_enriched_pc1_explained = NaN;
row.rare_enriched_cum_compared_explained = NaN;
row.combined_vs_base_subspace_similarity = NaN;
row.combined_vs_rare_enriched_subspace_similarity = NaN;
row.base_vs_rare_enriched_subspace_similarity = NaN;
row.combined_vs_base_pc1_abs_alignment = NaN;
row.combined_vs_rare_enriched_pc1_abs_alignment = NaN;
row.base_vs_rare_enriched_pc1_abs_alignment = NaN;
row.audit_status = "not_applicable_no_anchor_stage";

if ~ismember('anchor_stage', S.rowMeta.Properties.VariableNames)
    row = local_finish_anchor_stage_audit(row);
    return
end
stage = string(S.rowMeta.anchor_stage);
base = stage == "base_time_even";
rare = stage == "rare_strata_enriched";
row.n_base_rows = nnz(base);
row.n_rare_enriched_rows = nnz(rare);
if nnz(base) < 5 || nnz(rare) < 5
    row.audit_status = "not_applicable_requires_base_and_rare_enriched_rows";
    row = local_finish_anchor_stage_audit(row);
    return
end

k = min([params.enrichment_sensitivity_n_pcs_compared, nPC, ...
    nnz(base) - 1, nnz(rare) - 1, size(Xproc, 2), size(combinedCoeff, 2)]);
if k < 1
    row.audit_status = "not_applicable_no_comparable_pcs";
    row = local_finish_anchor_stage_audit(row);
    return
end

[baseCoeff, ~, baseLatent] = pca(Xproc(base, :), 'NumComponents', k);
[rareCoeff, ~, rareLatent] = pca(Xproc(rare, :), 'NumComponents', k);
baseExplained = local_explained_from_latent(baseLatent, Xproc(base, :));
rareExplained = local_explained_from_latent(rareLatent, Xproc(rare, :));

row.n_pcs_compared = k;
row.combined_cum_compared_explained = sum(combinedExplained(1:min(k, numel(combinedExplained))));
row.base_pc1_explained = baseExplained(1);
row.base_cum_compared_explained = sum(baseExplained(1:min(k, numel(baseExplained))));
row.rare_enriched_pc1_explained = rareExplained(1);
row.rare_enriched_cum_compared_explained = sum(rareExplained(1:min(k, numel(rareExplained))));
row.combined_vs_base_subspace_similarity = local_subspace_similarity(combinedCoeff, baseCoeff, k);
row.combined_vs_rare_enriched_subspace_similarity = local_subspace_similarity(combinedCoeff, rareCoeff, k);
row.base_vs_rare_enriched_subspace_similarity = local_subspace_similarity(baseCoeff, rareCoeff, k);
row.combined_vs_base_pc1_abs_alignment = abs(combinedCoeff(:, 1)' * baseCoeff(:, 1));
row.combined_vs_rare_enriched_pc1_abs_alignment = abs(combinedCoeff(:, 1)' * rareCoeff(:, 1));
row.base_vs_rare_enriched_pc1_abs_alignment = abs(baseCoeff(:, 1)' * rareCoeff(:, 1));
row.audit_status = "complete";
row = local_finish_anchor_stage_audit(row);
end

function row = local_finish_anchor_stage_audit(row)
row.preprocessing_fit_scope = "combined_condition_blind_anchor_bank";
row.comparison_fit_scope = "audit_only_base_and_rare_enriched_unweighted_pca";
row.anchor_stage_used_for_primary_pca = false;
row.audit_only_not_selection = true;
row.labels_used_for_stage_sensitivity = "none";
row.arena_used_for_stage_sensitivity = false;
row.condition_used_for_stage_sensitivity = false;
end

function similarity = local_subspace_similarity(coeffA, coeffB, k)
singularVals = svd(coeffA(:, 1:k)' * coeffB(:, 1:k));
similarity = mean(singularVals, 'omitnan');
end

function values = local_table_string_column(T, name, n)
if ismember(name, T.Properties.VariableNames)
    values = string(T.(name));
else
    values = repmat("", n, 1);
end
values = values(:);
end

function values = local_table_numeric_column(T, name, n)
if ismember(name, T.Properties.VariableNames)
    values = double(T.(name));
else
    values = nan(n, 1);
end
values = values(:);
end

function value = local_table_string_value(T, name, idx)
values = local_table_string_column(T, name, height(T));
value = values(idx);
end

function value = local_table_numeric_value(T, name, idx)
values = local_table_numeric_column(T, name, height(T));
value = values(idx);
end

function ess = local_effective_sample_size(w)
w = double(w(:));
w = w(isfinite(w) & w > 0);
if isempty(w)
    ess = NaN;
else
    ess = sum(w).^2 ./ sum(w.^2);
end
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
