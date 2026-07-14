function params = load_multiscale_graph_config(configPath)
%LOAD_MULTISCALE_GRAPH_CONFIG Load run-08 graph/topology parameters.

if nargin < 1 || strlength(string(configPath)) == 0
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    configPath = fullfile(repoRoot, 'config', 'multiscale_graph_config.csv');
end

opts = detectImportOptions(configPath, 'FileType', 'text', ...
    'Delimiter', ',', 'TextType', 'string');
opts.DataLines = [2 Inf];
opts.VariableNamesLine = 1;
opts.VariableUnitsLine = 0;
opts = setvartype(opts, opts.VariableNames, 'string');
T = readtable(configPath, opts);

required = ["parameter", "value", "type", "env_override", "description"];
missing = setdiff(required, string(T.Properties.VariableNames));
assert(isempty(missing), 'load_multiscale_graph_config:MissingColumn', ...
    'Config %s is missing required columns: %s', configPath, strjoin(missing, ', '));

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
canonicalSmokeOutput = "derived/graph_motif_discovery_smoke";
canonicalFullOutput = "derived/graph_motif_discovery";
canonicalSmokeEmbedding = "derived/embedding_motif_discovery_smoke";
canonicalFullEmbedding = "derived/embedding_motif_discovery";

outputExplicit = local_env_override_used(T, "output_dir");
embeddingExplicit = local_env_override_used(T, "embedding_input_dir");

if params.run_mode == "full"
    if ~outputExplicit && any(params.output_dir == [canonicalSmokeOutput, canonicalFullOutput])
        params.output_dir = canonicalFullOutput;
    end
    if ~embeddingExplicit && any(params.embedding_input_dir == [canonicalSmokeEmbedding, canonicalFullEmbedding])
        params.embedding_input_dir = canonicalFullEmbedding;
    end
else
    if ~outputExplicit && any(params.output_dir == [canonicalSmokeOutput, canonicalFullOutput])
        params.output_dir = canonicalSmokeOutput;
    end
    if ~embeddingExplicit && any(params.embedding_input_dir == [canonicalSmokeEmbedding, canonicalFullEmbedding])
        params.embedding_input_dir = canonicalSmokeEmbedding;
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
    'load_multiscale_graph_config:BadRunMode', ...
    'run_mode must be smoke or full.');
assert(params.graph_n_pcs >= 2, ...
    'load_multiscale_graph_config:BadGraphPCs', ...
    'graph_n_pcs must be at least 2.');
assert(params.k_neighbors >= 2, ...
    'load_multiscale_graph_config:BadK', ...
    'k_neighbors must be at least 2.');
assert(all(params.k_sensitivity_values >= 2), ...
    'load_multiscale_graph_config:BadKValues', ...
    'k_sensitivity_values must all be at least 2.');
assert(params.knn_block_size >= 16, ...
    'load_multiscale_graph_config:BadBlockSize', ...
    'knn_block_size must be at least 16.');
assert(isfinite(params.graph_score_winsor_abs) && params.graph_score_winsor_abs > 0, ...
    'load_multiscale_graph_config:BadWinsorAbs', ...
    'graph_score_winsor_abs must be positive and finite.');
assert(params.min_event_valid_fraction >= 0 && params.min_event_valid_fraction <= 1, ...
    'load_multiscale_graph_config:BadEventValidFraction', ...
    'min_event_valid_fraction must be in [0, 1].');
assert(params.high_quantile_threshold > 0 && params.high_quantile_threshold < 1, ...
    'load_multiscale_graph_config:BadHighQuantile', ...
    'high_quantile_threshold must be in (0, 1).');
local_validate_run_mode_paths(params);
end

function local_validate_run_mode_paths(params)
outputDir = lower(replace(string(params.output_dir), "\", "/"));
embeddingDir = lower(replace(string(params.embedding_input_dir), "\", "/"));
if params.run_mode == "full"
    assert(~contains(outputDir, "smoke"), ...
        'load_multiscale_graph_config:BadFullOutputDir', ...
        'Full run_08 output_dir must not point to a smoke directory.');
    assert(~contains(embeddingDir, "smoke"), ...
        'load_multiscale_graph_config:BadFullEmbeddingInputDir', ...
        'Full run_08 embedding_input_dir must not point to a smoke run_07 directory.');
else
    assert(contains(outputDir, "smoke"), ...
        'load_multiscale_graph_config:BadSmokeOutputDir', ...
        'Smoke run_08 output_dir must point to a smoke directory.');
    assert(contains(embeddingDir, "smoke"), ...
        'load_multiscale_graph_config:BadSmokeEmbeddingInputDir', ...
        'Smoke run_08 embedding_input_dir must point to a smoke run_07 directory.');
end
end

function value = local_parse_value(rawValue, valueType, parameterName)
rawValue = string(rawValue);
switch valueType
    case "number"
        value = str2double(rawValue);
        assert(isfinite(value) || isnan(value), ...
            'load_multiscale_graph_config:BadNumber', ...
            '%s must be numeric.', parameterName);
    case "number_list"
        parts = string(split(rawValue, {';', '|', ','}))';
        parts = strtrim(parts);
        parts = parts(parts ~= "");
        value = str2double(parts);
        assert(~isempty(value) && all(isfinite(value)), ...
            'load_multiscale_graph_config:BadNumberList', ...
            '%s must be a finite numeric list.', parameterName);
    case "logical"
        raw = lower(strtrim(rawValue));
        assert(any(raw == ["true", "false", "1", "0", "yes", "no"]), ...
            'load_multiscale_graph_config:BadLogical', ...
            '%s must be logical.', parameterName);
        value = any(raw == ["true", "1", "yes"]);
    case "string"
        value = rawValue;
    otherwise
        error('load_multiscale_graph_config:BadType', ...
            'Unknown config type "%s" for parameter "%s".', valueType, parameterName);
end
end

function value = local_value_to_string(x)
if islogical(x)
    value = string(x);
elseif isnumeric(x)
    if isscalar(x)
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
