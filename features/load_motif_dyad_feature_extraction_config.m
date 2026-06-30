function params = load_motif_dyad_feature_extraction_config(configPath)
%LOAD_MOTIF_DYAD_FEATURE_EXTRACTION_CONFIG Load run-05 feature parameters.
%
% Defaults live in config/motif_dyad_feature_extraction_config.csv. Optional
% environment variables listed in that file may override values for a run.

if nargin < 1 || strlength(string(configPath)) == 0
    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    configPath = fullfile(repoRoot, 'config', 'motif_dyad_feature_extraction_config.csv');
end

importOpts = detectImportOptions(configPath, 'FileType','text', ...
    'Delimiter', ',', 'TextType','string');
importOpts.DataLines = [2 Inf];
importOpts.VariableNamesLine = 1;
importOpts.VariableUnitsLine = 0;
importOpts = setvartype(importOpts, importOpts.VariableNames, 'string');
T = readtable(configPath, importOpts);
assert(height(T) >= 1, 'load_motif_dyad_feature_extraction_config:EmptyConfig', ...
    'Config %s must contain at least one parameter row.', configPath);
required = ["parameter","value","type","env_override","description"];
missing = setdiff(required, string(T.Properties.VariableNames));
assert(isempty(missing), 'load_motif_dyad_feature_extraction_config:MissingColumn', ...
    'Config %s is missing required columns: %s', configPath, strjoin(missing, ', '));

params = struct();
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
        end
    end
    params.(name) = local_parse_value(rawValue, lower(strtrim(string(T.type(i)))), T.parameter(i));
end

params.config_path = string(configPath);
params.config_table = T;

assert(any(params.run_mode == ["smoke","full"]), ...
    'load_motif_dyad_feature_extraction_config:BadRunMode', ...
    'run_mode must be smoke or full.');
assert(params.smoke_sessions_per_arena >= 1, ...
    'load_motif_dyad_feature_extraction_config:BadSmokeCount', ...
    'smoke_sessions_per_arena must be positive.');
assert(params.fps > 0 && params.pixel_size_mm > 0, ...
    'load_motif_dyad_feature_extraction_config:BadUnitScale', ...
    'fps and pixel_size_mm must be positive.');
assert(params.unit_parameter_tolerance >= 0, ...
    'load_motif_dyad_feature_extraction_config:BadTolerance', ...
    'unit_parameter_tolerance must be nonnegative.');
assert(params.contact_threshold_mm > 0 && params.close_threshold_mm > 0, ...
    'load_motif_dyad_feature_extraction_config:BadThreshold', ...
    'contact_threshold_mm and close_threshold_mm must be positive.');
assert(params.close_threshold_mm >= params.contact_threshold_mm, ...
    'load_motif_dyad_feature_extraction_config:ThresholdOrder', ...
    'close_threshold_mm must be >= contact_threshold_mm.');
end

function value = local_parse_value(rawValue, valueType, parameterName)
switch valueType
    case "number"
        if lower(rawValue) == "nan"
            value = NaN;
        else
            value = str2double(rawValue);
        end
        assert(isfinite(value) || isnan(value), ...
            'load_motif_dyad_feature_extraction_config:BadNumber', ...
            '%s must be numeric.', parameterName);
    case "logical"
        raw = lower(strtrim(rawValue));
        assert(any(raw == ["true","false","1","0","yes","no"]), ...
            'load_motif_dyad_feature_extraction_config:BadLogical', ...
            '%s must be logical.', parameterName);
        value = any(raw == ["true","1","yes"]);
    case "string"
        value = string(rawValue);
    case "string_list"
        if lower(strtrim(rawValue)) == "all"
            value = "all";
        else
            value = string(split(rawValue, {';','|',','}))';
            value = strtrim(value);
            value = value(value ~= "");
        end
    otherwise
        error('load_motif_dyad_feature_extraction_config:BadType', ...
            'Unknown config type "%s" for parameter "%s".', valueType, parameterName);
end
end
