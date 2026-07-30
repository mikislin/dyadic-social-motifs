function params = load_motif_candidate_behavioral_validation_config(configPath)
%LOAD_MOTIF_CANDIDATE_BEHAVIORAL_VALIDATION_CONFIG Load run_10 config.

if nargin < 1 || strlength(string(configPath)) == 0
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    configPath = fullfile(repoRoot, 'config', ...
        'motif_candidate_behavioral_validation_config.csv');
end

opts = detectImportOptions(configPath, 'FileType', 'text', ...
    'Delimiter', ',', 'TextType', 'string');
opts.DataLines = [2, Inf];
opts = setvartype(opts, opts.VariableNames, 'string');
T = readtable(configPath, opts);
required = ["parameter","value","type","env_override","description"];
missing = setdiff(required, string(T.Properties.VariableNames));
assert(isempty(missing), ...
    'load_motif_candidate_behavioral_validation_config:MissingColumn', ...
    'Config is missing columns: %s', strjoin(missing, ', '));
assert(numel(unique(T.parameter)) == height(T), ...
    'load_motif_candidate_behavioral_validation_config:DuplicateParameter', ...
    'Every run_10 config parameter must be unique.');

params = struct();
envUsed = false(height(T), 1);
for i = 1:height(T)
    name = matlab.lang.makeValidName(T.parameter(i));
    valueText = string(T.value(i));
    envName = strtrim(string(T.env_override(i)));
    if ismissing(envName)
        envName = "";
    end
    if envName ~= ""
        envValue = strtrim(string(getenv(envName)));
        if envValue ~= ""
            valueText = envValue;
            envUsed(i) = true;
        end
    end
    params.(name) = i_parse(valueText, lower(T.type(i)), T.parameter(i));
end

T.env_override_used = envUsed;
[params, T] = i_apply_mode(params, T);
T.effective_value = strings(height(T), 1);
for i = 1:height(T)
    name = matlab.lang.makeValidName(T.parameter(i));
    T.effective_value(i) = i_to_string(params.(name));
end
params.config_path = string(configPath);
params.config_table = T;
params.feature_specs = i_feature_specs(params);
i_validate(params);
end

function [params, T] = i_apply_mode(params, T)
smokeRoot = "derived/motif_candidate_behavioral_validation_smoke";
fullRoot = "derived/motif_candidate_behavioral_validation";
explicitOutput = any(T.parameter == "output_dir" & T.env_override_used);
if params.run_mode == "full"
    if ~explicitOutput && any(params.output_dir == [smokeRoot fullRoot])
        params.output_dir = fullRoot;
    end
    params.active_permutation_count = params.permutation_count_full;
    params.active_bootstrap_count = params.bootstrap_count_full;
    params.active_maximum_rows_per_cell = ...
        params.maximum_rows_per_candidate_session_scale_full;
    params.active_minimum_nodes_primary = params.minimum_nodes_primary_full;
    params.active_smoke_session_count = Inf;
else
    if ~explicitOutput && any(params.output_dir == [smokeRoot fullRoot])
        params.output_dir = smokeRoot;
    end
    params.active_permutation_count = params.permutation_count_smoke;
    params.active_bootstrap_count = params.bootstrap_count_smoke;
    params.active_maximum_rows_per_cell = ...
        params.maximum_rows_per_candidate_session_scale_smoke;
    params.active_minimum_nodes_primary = params.minimum_nodes_primary_smoke;
    params.active_smoke_session_count = params.smoke_session_count;
end
end

function Specs = i_feature_specs(params)
names = string(fieldnames(params));
names = sort(names(startsWith(names, "feature_") & endsWith(names, "_spec")));
n = numel(names);
featureId = strings(n, 1);
featureName = strings(n, 1);
unit = strings(n, 1);
equation = strings(n, 1);
rawNodes = strings(n, 1);
lineageClass = strings(n, 1);
rationale = strings(n, 1);
for i = 1:n
    featureId(i) = upper(extractBetween(names(i), "feature_", "_spec"));
    parts = string(split(string(params.(names(i))), "|"));
    assert(numel(parts) == 6, ...
        'load_motif_candidate_behavioral_validation_config:BadFeatureSpec', ...
        'Feature spec %s must have six pipe-delimited fields.', names(i));
    featureName(i) = parts(1);
    unit(i) = parts(2);
    equation(i) = parts(3);
    rawNodes(i) = parts(4);
    lineageClass(i) = parts(5);
    rationale(i) = parts(6);
end
Specs = table(featureId, featureName, unit, equation, rawNodes, ...
    lineageClass, rationale, ...
    'VariableNames', {'feature_id','feature_name','unit','equation', ...
    'raw_pose_nodes','lineage_class','inclusion_rationale'});
end

function i_validate(params)
assert(any(params.run_mode == ["smoke","full"]), ...
    'load_motif_candidate_behavioral_validation_config:BadRunMode', ...
    'run_mode must be smoke or full.');
assert(any(params.phase == ["prepare","analyze","full"]), ...
    'load_motif_candidate_behavioral_validation_config:BadPhase', ...
    'phase must be prepare, analyze, or full.');
output = lower(replace(params.output_dir, "\", "/"));
if params.run_mode == "full"
    assert(~contains(output, "smoke"), ...
        'load_motif_candidate_behavioral_validation_config:BadFullOutput', ...
        'Full output root cannot contain smoke.');
else
    assert(contains(output, "smoke"), ...
        'load_motif_candidate_behavioral_validation_config:BadSmokeOutput', ...
        'Smoke output must use a separate smoke root.');
end
assert(params.expected_candidate_freeze_id == ...
    "run09_candidates_2e58d214683b_8817d942_0e956a2c", ...
    'load_motif_candidate_behavioral_validation_config:BadFreezeId', ...
    'The issued run_09 freeze ID is immutable.');
assert(strlength(params.expected_membership_sha256) == 64 && ...
    ~isempty(regexp(char(params.expected_membership_sha256), ...
    '^[0-9a-f]{64}$', 'once')), ...
    'load_motif_candidate_behavioral_validation_config:BadHash', ...
    'Expected membership SHA-256 must be lowercase hexadecimal.');
assert(numel(params.expected_eligible_candidate_ids) == 9 && ...
    numel(unique(params.expected_eligible_candidate_ids)) == 9, ...
    'load_motif_candidate_behavioral_validation_config:BadEligibleSet', ...
    'Exactly nine unique graph-eligible candidate IDs are required.');
assert(isempty(intersect(params.discovery_value_pose_nodes, ...
    params.validation_pose_nodes)), ...
    'load_motif_candidate_behavioral_validation_config:PoseNodeOverlap', ...
    'Validation pose nodes must be disjoint from discovery-value pose nodes.');
assert(height(params.feature_specs) == 16 && ...
    numel(unique(params.feature_specs.feature_name)) == 16 && ...
    all(params.feature_specs.lineage_class == ...
    "independent_behavioral_measure"), ...
    'load_motif_candidate_behavioral_validation_config:BadFeaturePanel', ...
    'The frozen v1 panel requires 16 unique independent measures.');
assert(params.validation_window_sec > 0 && ...
    params.minimum_feature_finite_fraction > 0 && ...
    params.minimum_feature_finite_fraction <= 1 && ...
    params.minimum_correlation_frames >= 3 && ...
    params.minimum_finite_feature_count >= 1 && ...
    params.minimum_finite_feature_count <= height(params.feature_specs), ...
    'load_motif_candidate_behavioral_validation_config:BadFeatureRules', ...
    'Feature window or finite-data rules are invalid.');
assert(params.grouped_fold_count >= 2 && ...
    params.active_permutation_count >= 1 && ...
    params.active_maximum_rows_per_cell >= 1, ...
    'load_motif_candidate_behavioral_validation_config:BadDesignCounts', ...
    'Fold, permutation, and sampling counts must be positive.');
assert(params.random_stream_algorithm == "mt19937ar", ...
    'load_motif_candidate_behavioral_validation_config:BadRng', ...
    'The validated random stream is mt19937ar.');
assert(params.supported_required_primary_gates == 3 && ...
    params.partial_required_primary_gates == 2, ...
    'load_motif_candidate_behavioral_validation_config:BadStatusRule', ...
    'Status v1 requires three gates for supported and two for partial.');
end

function value = i_parse(textValue, valueType, parameter)
textValue = string(textValue);
switch valueType
    case "number"
        if lower(strtrim(textValue)) == "inf"
            value = Inf;
        else
            value = str2double(textValue);
        end
        assert(~isnan(value), ...
            'load_motif_candidate_behavioral_validation_config:BadNumber', ...
            '%s must be numeric.', parameter);
    case "logical"
        x = lower(strtrim(textValue));
        assert(any(x == ["true","false","1","0","yes","no"]), ...
            'load_motif_candidate_behavioral_validation_config:BadLogical', ...
            '%s must be logical-like.', parameter);
        value = any(x == ["true","1","yes"]);
    case "string_list"
        value = strtrim(string(split(textValue, ";")))';
        value = value(value ~= "");
    case {"string","path"}
        value = strtrim(textValue);
    otherwise
        error('load_motif_candidate_behavioral_validation_config:BadType', ...
            'Unsupported type %s for %s.', valueType, parameter);
end
end

function value = i_to_string(x)
if islogical(x)
    value = string(x);
elseif isnumeric(x)
    value = strjoin(compose('%.17g', x(:)'), ';');
elseif isstring(x) && numel(x) > 1
    value = strjoin(x(:)', ';');
else
    value = string(x);
end
end
