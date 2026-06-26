function cfg = load_preprocessing_pipeline_config(configPath)
%LOAD_PREPROCESSING_PIPELINE_CONFIG Read preprocessing/QC settings.

if nargin < 1 || isempty(configPath)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    configPath = fullfile(repoRoot, 'config', 'preprocessing_pipeline_config.csv');
end

assert(isfile(configPath), 'load_preprocessing_pipeline_config:MissingConfig', ...
    'Missing preprocessing pipeline config: %s', configPath);

T = readtable(configPath, 'TextType', 'string');
required = ["section","key","value","type"];
missing = setdiff(required, string(T.Properties.VariableNames));
assert(isempty(missing), 'load_preprocessing_pipeline_config:BadSchema', ...
    'Missing required config columns: %s', strjoin(missing, ', '));

cfg = struct();
cfg.configPath = string(configPath);
for i = 1:height(T)
    section = matlab.lang.makeValidName(T.section(i));
    key = matlab.lang.makeValidName(T.key(i));
    cfg.(section).(key) = local_parse_value(T.value(i), T.type(i));
end

cfg = local_attach_bodypoint_metadata(cfg, configPath);
end

function value = local_parse_value(rawValue, rawType)
rawValue = strtrim(string(rawValue));
rawType = lower(strtrim(string(rawType)));

switch rawType
    case "string"
        value = rawValue;
    case "double"
        value = str2double(rawValue);
        assert(isfinite(value), 'load_preprocessing_pipeline_config:BadDouble', ...
            'Could not parse numeric config value: %s', rawValue);
    case "integer"
        value = str2double(rawValue);
        assert(isfinite(value) && mod(value, 1) == 0, ...
            'load_preprocessing_pipeline_config:BadInteger', ...
            'Could not parse integer config value: %s', rawValue);
    case "logical"
        value = any(strcmpi(rawValue, ["true","1","yes"]));
        assert(value || any(strcmpi(rawValue, ["false","0","no"])), ...
            'load_preprocessing_pipeline_config:BadLogical', ...
            'Could not parse logical config value: %s', rawValue);
    case "string_vector"
        value = string(strtrim(split(rawValue, ';')))';
        value = value(value ~= "");
    case "double_vector"
        parts = string(strtrim(split(rawValue, ';')))';
        value = str2double(parts);
        assert(all(isfinite(value)), 'load_preprocessing_pipeline_config:BadDoubleVector', ...
            'Could not parse numeric vector config value: %s', rawValue);
    otherwise
        error('load_preprocessing_pipeline_config:UnknownType', ...
            'Unknown config type: %s', rawType);
end
end

function cfg = local_attach_bodypoint_metadata(cfg, configPath)
if ~isfield(cfg, 'bodypoints') || ...
        ~isfield(cfg.bodypoints, 'bodypoint_metadata_file')
    return
end

metadataPath = char(cfg.bodypoints.bodypoint_metadata_file);
if ~is_absolute_path(metadataPath)
    metadataPath = fullfile(fileparts(configPath), metadataPath);
end
assert(isfile(metadataPath), 'load_preprocessing_pipeline_config:MissingBodypointMetadata', ...
    'Missing SLEAP bodypart metadata file: %s', metadataPath);

M = readtable(metadataPath, 'TextType', 'string');
required = ["node_index","bodypart_name","display_label","matlab_field"];
missing = setdiff(required, string(M.Properties.VariableNames));
assert(isempty(missing), 'load_preprocessing_pipeline_config:BadBodypointMetadata', ...
    'Missing bodypoint metadata columns: %s', strjoin(missing, ', '));

M = sortrows(M, 'node_index');
nodeIndex = M.node_index;
assert(all(isfinite(nodeIndex)) && all(mod(nodeIndex, 1) == 0) && ...
    all(nodeIndex(:)' == 1:height(M)), ...
    'load_preprocessing_pipeline_config:BadBodypointMetadata', ...
    'Bodypoint node_index values must be consecutive integers starting at 1.');
assert(numel(unique(M.bodypart_name)) == height(M), ...
    'load_preprocessing_pipeline_config:BadBodypointMetadata', ...
    'Bodypoint names must be unique.');
assert(numel(unique(M.matlab_field)) == height(M), ...
    'load_preprocessing_pipeline_config:BadBodypointMetadata', ...
    'Bodypoint MATLAB field names must be unique.');

cfg.bodypoints.bodypoint_metadata_path = string(metadataPath);
cfg.bodypoints.bodypoint_metadata = M;
cfg.bodypoints.bodypoint_names = string(M.bodypart_name(:))';
cfg.bodypoints.bodypoint_labels = string(M.display_label(:))';
cfg.bodypoints.bodypoint_matlab_fields = string(M.matlab_field(:))';
end

function tf = is_absolute_path(pathText)
pathText = char(pathText);
tf = startsWith(pathText, filesep) || ...
    ~isempty(regexp(pathText, '^[A-Za-z]:[\\/]', 'once'));
end
