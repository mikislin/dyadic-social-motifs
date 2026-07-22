function params = load_multiscale_embedding_config(configPath)
%LOAD_MULTISCALE_EMBEDDING_CONFIG Load run-07 embedding parameters.

if nargin < 1 || strlength(string(configPath)) == 0
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    configPath = fullfile(repoRoot, 'config', 'multiscale_embedding_config.csv');
end

importOpts = detectImportOptions(configPath, 'FileType', 'text', ...
    'Delimiter', ',', 'TextType', 'string');
importOpts.DataLines = [2 Inf];
importOpts.VariableNamesLine = 1;
importOpts.VariableUnitsLine = 0;
importOpts = setvartype(importOpts, importOpts.VariableNames, 'string');
T = readtable(configPath, importOpts);

required = ["parameter", "value", "type", "env_override", "description"];
missing = setdiff(required, string(T.Properties.VariableNames));
assert(isempty(missing), 'load_multiscale_embedding_config:MissingColumn', ...
    'Config %s is missing required columns: %s', configPath, strjoin(missing, ', '));
assert(height(T) >= 1, 'load_multiscale_embedding_config:EmptyConfig', ...
    'Config %s must contain at least one parameter row.', configPath);

params = struct();
effectiveValue = strings(height(T), 1);
envOverrideUsed = false(height(T), 1);
for i = 1:height(T)
    name = matlab.lang.makeValidName(T.parameter(i));
    rawValue = string(T.value(i));
    envName = strtrim(string(T.env_override(i)));
    if ismissing(envName)
        envName = "";
    end
    if envName ~= ""
        envValue = strtrim(string(getenv(envName)));
        if envValue ~= ""
            rawValue = envValue;
            envOverrideUsed(i) = true;
        end
    end
    params.(name) = local_parse_value(rawValue, lower(strtrim(string(T.type(i)))), T.parameter(i));
    effectiveValue(i) = local_value_to_string(params.(name));
end

T.env_override_used = envOverrideUsed;
[params, T] = local_apply_mode_defaults(params, T);
T.effective_value = effectiveValue;
T = local_refresh_effective_values(T, params);
T.config_path = repmat(string(configPath), height(T), 1);
params.config_path = string(configPath);
params.config_table = T;

local_validate_params(params);
end

function [params, T] = local_apply_mode_defaults(params, T)
canonicalSmokeOutput = "derived/embedding_motif_discovery_smoke";
canonicalFullOutput = "derived/embedding_motif_discovery";
canonicalRareSmokeOutput = "derived/embedding_motif_discovery_rare_enriched_smoke";
canonicalRareFullOutput = "derived/embedding_motif_discovery_rare_enriched";
canonicalChunkInput = "derived/chunks_motif_discovery";
canonicalRareChunkSmokeInput = "derived/chunks_motif_discovery_rare_enriched_smoke";
canonicalRareChunkFullInput = "derived/chunks_motif_discovery_rare_enriched";

outputExplicit = local_env_override_used(T, "output_dir");
inputExplicit = local_env_override_used(T, "chunk_input_dir");

if params.anchor_manifest_mode == "rare_enriched"
    params.allow_reviewed_snapshot_fallback = false;
    if params.run_mode == "full"
        if ~outputExplicit && any(params.output_dir == [canonicalSmokeOutput, canonicalFullOutput, canonicalRareSmokeOutput, canonicalRareFullOutput])
            params.output_dir = canonicalRareFullOutput;
        end
        if ~inputExplicit && any(params.chunk_input_dir == [canonicalChunkInput, canonicalRareChunkSmokeInput, canonicalRareChunkFullInput])
            params.chunk_input_dir = canonicalRareChunkFullInput;
        end
    else
        if ~outputExplicit && any(params.output_dir == [canonicalSmokeOutput, canonicalFullOutput, canonicalRareSmokeOutput, canonicalRareFullOutput])
            params.output_dir = canonicalRareSmokeOutput;
        end
        if ~inputExplicit && any(params.chunk_input_dir == [canonicalChunkInput, canonicalRareChunkSmokeInput, canonicalRareChunkFullInput])
            params.chunk_input_dir = canonicalRareChunkSmokeInput;
        end
    end
    return
end

if params.run_mode == "full"
    if ~outputExplicit && any(params.output_dir == [canonicalSmokeOutput, canonicalFullOutput])
        params.output_dir = canonicalFullOutput;
    end
else
    if ~outputExplicit && any(params.output_dir == [canonicalSmokeOutput, canonicalFullOutput])
        params.output_dir = canonicalSmokeOutput;
    end
end
end

function tf = local_env_override_used(T, parameterName)
idx = string(T.parameter) == string(parameterName);
tf = any(T.env_override_used(idx));
end

function T = local_refresh_effective_values(T, params)
effectiveValue = strings(height(T), 1);
for i = 1:height(T)
    name = matlab.lang.makeValidName(T.parameter(i));
    if isfield(params, name)
        effectiveValue(i) = local_value_to_string(params.(name));
    end
end
T.effective_value = effectiveValue;
end

function local_validate_params(params)
assert(any(params.run_mode == ["smoke", "full"]), ...
    'load_multiscale_embedding_config:BadRunMode', ...
    'run_mode must be smoke or full.');
assert(any(params.anchor_manifest_mode == ["primary", "rare_enriched"]), ...
    'load_multiscale_embedding_config:BadAnchorManifestMode', ...
    'anchor_manifest_mode must be primary or rare_enriched.');
assert(~logical(params.use_anchor_weights_for_pca), ...
    'load_multiscale_embedding_config:WeightedPcaNotYetEnabled', ...
    'Initial enriched run_07 must fit unweighted PCA and use anchor weights for audit only.');
assert(any(lower(params.anchor_mode) == ["center", "past"]), ...
    'load_multiscale_embedding_config:BadAnchorMode', ...
    'anchor_mode must be center or past.');
assert(params.fps > 0, 'load_multiscale_embedding_config:BadFps', ...
    'fps must be positive.');
assert(params.smoke_max_anchors_per_scale >= 1, ...
    'load_multiscale_embedding_config:BadSmokeAnchorCap', ...
    'smoke_max_anchors_per_scale must be positive.');
assert(isnan(params.max_anchors_per_scale) || params.max_anchors_per_scale >= 1, ...
    'load_multiscale_embedding_config:BadFullAnchorCap', ...
    'max_anchors_per_scale must be NaN or positive.');
assert(params.summary_temporal_bins >= 1 && params.summary_dct_coeffs >= 0, ...
    'load_multiscale_embedding_config:BadSummaryParams', ...
    'summary_temporal_bins must be positive and summary_dct_coeffs nonnegative.');
assert(params.preprocess_winsor_low >= 0 && params.preprocess_winsor_high <= 1 && ...
    params.preprocess_winsor_low < params.preprocess_winsor_high, ...
    'load_multiscale_embedding_config:BadWinsor', ...
    'winsor quantiles must lie in [0, 1] and be ordered.');
assert(params.preprocess_min_robust_scale > 0, ...
    'load_multiscale_embedding_config:BadMinRobustScale', ...
    'preprocess_min_robust_scale must be positive.');
assert(isfinite(params.preprocess_min_iqr_to_std_ratio) && ...
    params.preprocess_min_iqr_to_std_ratio >= 0, ...
    'load_multiscale_embedding_config:BadIqrStdRatio', ...
    'preprocess_min_iqr_to_std_ratio must be finite and nonnegative.');
assert(isfinite(params.preprocess_tail_audit_abs_threshold) && ...
    params.preprocess_tail_audit_abs_threshold > 0 && ...
    isfinite(params.preprocess_severe_tail_audit_abs_threshold) && ...
    params.preprocess_severe_tail_audit_abs_threshold > params.preprocess_tail_audit_abs_threshold, ...
    'load_multiscale_embedding_config:BadTailAuditThresholds', ...
    'Preprocessing tail thresholds must be finite positive and severe must exceed the primary threshold.');
assert(params.min_pcs_per_scale >= 1 && params.micro_max_pcs >= params.min_pcs_per_scale && ...
    params.motif_max_pcs >= params.min_pcs_per_scale && params.context_max_pcs >= params.min_pcs_per_scale, ...
    'load_multiscale_embedding_config:BadPcCaps', ...
    'PC caps must be positive and not below min_pcs_per_scale.');
assert(params.global_n_pcs >= 1, 'load_multiscale_embedding_config:BadGlobalPCs', ...
    'global_n_pcs must be positive.');
assert(params.global_matrix_mode == "ordinal_pc_stack", ...
    'load_multiscale_embedding_config:BadGlobalMatrixMode', ...
    'global_matrix_mode currently must be ordinal_pc_stack.');
assert(isfinite(params.scale_score_winsor_abs) && params.scale_score_winsor_abs > 0, ...
    'load_multiscale_embedding_config:BadScaleScoreWinsor', ...
    'scale_score_winsor_abs must be positive and finite.');
assert(any(params.scale_weight_mode == ["equal_total_weight", "none"]), ...
    'load_multiscale_embedding_config:BadScaleWeightMode', ...
    'scale_weight_mode must be equal_total_weight or none.');
assert(params.stability_splits >= 1 && params.smoke_stability_splits >= 1 && ...
    params.stability_n_pcs_compared >= 1, ...
    'load_multiscale_embedding_config:BadStabilityParams', ...
    'Stability split counts and compared PC count must be positive.');
assert(params.enrichment_sensitivity_n_pcs_compared >= 1, ...
    'load_multiscale_embedding_config:BadEnrichmentSensitivityPcs', ...
    'enrichment_sensitivity_n_pcs_compared must be positive.');
assert(params.arena_sensitivity_n_pcs >= 1, ...
    'load_multiscale_embedding_config:BadArenaSensitivityPCs', ...
    'arena_sensitivity_n_pcs must be positive.');
local_validate_run_mode_paths(params);
end

function local_validate_run_mode_paths(params)
outputDir = lower(replace(string(params.output_dir), "\", "/"));
chunkInputDir = lower(replace(string(params.chunk_input_dir), "\", "/"));
fallbackDir = lower(replace(string(params.fallback_chunk_input_dir), "\", "/"));
if params.run_mode == "full"
    assert(~contains(outputDir, "smoke"), ...
        'load_multiscale_embedding_config:BadFullOutputDir', ...
        'Full run_07 output_dir must not point to a smoke directory.');
    assert(~contains(chunkInputDir, "smoke"), ...
        'load_multiscale_embedding_config:BadFullChunkInputDir', ...
        'Full run_07 chunk_input_dir must not point to a smoke run_06 directory.');
    assert(~contains(fallbackDir, "smoke"), ...
        'load_multiscale_embedding_config:BadFullFallbackInputDir', ...
        'Full run_07 fallback_chunk_input_dir must not point to a smoke run_06 directory.');
else
    assert(contains(outputDir, "smoke"), ...
        'load_multiscale_embedding_config:BadSmokeOutputDir', ...
        'Smoke run_07 output_dir must point to a smoke directory.');
end
if params.anchor_manifest_mode == "rare_enriched"
    assert(contains(outputDir, "rare_enriched") && contains(chunkInputDir, "rare_enriched"), ...
        'load_multiscale_embedding_config:RareModePathMismatch', ...
        'rare_enriched mode must use separate rare_enriched input and output roots.');
    assert(~logical(params.allow_reviewed_snapshot_fallback), ...
        'load_multiscale_embedding_config:RareFallbackBlocked', ...
        'rare_enriched mode cannot fall back to a reviewed primary snapshot.');
end
end

function value = local_parse_value(rawValue, valueType, parameterName)
rawValue = string(rawValue);
switch valueType
    case "number"
        if lower(strtrim(rawValue)) == "nan"
            value = NaN;
        else
            value = str2double(rawValue);
        end
        assert(isfinite(value) || isnan(value), ...
            'load_multiscale_embedding_config:BadNumber', ...
            '%s must be numeric.', parameterName);
    case "logical"
        raw = lower(strtrim(rawValue));
        assert(any(raw == ["true", "false", "1", "0", "yes", "no"]), ...
            'load_multiscale_embedding_config:BadLogical', ...
            '%s must be logical.', parameterName);
        value = any(raw == ["true", "1", "yes"]);
    case "string"
        value = rawValue;
    case "string_list"
        if lower(strtrim(rawValue)) == "all"
            value = "all";
        else
            value = string(split(rawValue, {';', '|', ','}))';
            value = strtrim(value);
            value = value(value ~= "");
        end
    otherwise
        error('load_multiscale_embedding_config:BadType', ...
            'Unknown config type "%s" for parameter "%s".', valueType, parameterName);
end
end

function value = local_value_to_string(x)
if islogical(x)
    value = string(x);
elseif isnumeric(x)
    if isscalar(x) && isnan(x)
        value = "NaN";
    elseif isscalar(x)
        value = string(sprintf('%.15g', x));
    else
        value = strjoin(string(compose('%.15g', x(:)')), ";");
    end
elseif isstring(x) || ischar(x)
    value = strjoin(string(x), ";");
else
    value = string(evalc('disp(x)'));
    value = strtrim(value);
end
end
