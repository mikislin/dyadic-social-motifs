function styles = load_experimental_group_styles(stylePath, manifest)
%LOAD_EXPERIMENTAL_GROUP_STYLES Read stable condition labels and colors.

if nargin < 1 || isempty(stylePath)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    stylePath = fullfile(repoRoot, 'config', 'experimental_group_styles.csv');
end

assert(isfile(stylePath), 'load_experimental_group_styles:MissingConfig', ...
    'Missing experimental-group style config: %s', stylePath);

T = readtable(stylePath, 'TextType', 'string');
required = ["condition_id","condition_label","display_label","color_hex","marker","plot_order"];
missing = setdiff(required, string(T.Properties.VariableNames));
assert(isempty(missing), 'load_experimental_group_styles:BadSchema', ...
    'Missing required style columns: %s', strjoin(missing, ', '));

assert(numel(unique(T.condition_id)) == height(T), ...
    'load_experimental_group_styles:DuplicateCondition', ...
    'Each condition_id must appear once in %s.', stylePath);

for i = 1:height(T)
    local_hex_to_rgb(T.color_hex(i));
end

if nargin >= 2 && ~isempty(manifest)
    if ismember('condition_id', manifest.Properties.VariableNames)
        manifestConditions = unique(string(manifest.condition_id));
    else
        manifestConditions = unique(string(manifest.condition));
    end
    missingConditions = setdiff(manifestConditions, string(T.condition_id));
    assert(isempty(missingConditions), ...
        'load_experimental_group_styles:MissingConditionStyle', ...
        'Missing style rows for condition_id values: %s', ...
        strjoin(missingConditions, ', '));
end

T = sortrows(T, 'plot_order');
styles = struct();
styles.path = string(stylePath);
styles.table = T;
styles.conditionIds = string(T.condition_id);
styles.displayLabels = string(T.display_label);
styles.markers = string(T.marker);
styles.colors = zeros(height(T), 3);
for i = 1:height(T)
    styles.colors(i,:) = local_hex_to_rgb(T.color_hex(i));
end
end

function rgb = local_hex_to_rgb(hex)
hex = char(strtrim(string(hex)));
assert(numel(hex) == 7 && hex(1) == '#', ...
    'load_experimental_group_styles:BadHexColor', ...
    'Hex colors must use #RRGGBB format: %s', hex);
rgb = sscanf(hex(2:end), '%2x%2x%2x', [1 3]) ./ 255;
assert(numel(rgb) == 3 && all(isfinite(rgb)), ...
    'load_experimental_group_styles:BadHexColor', ...
    'Hex colors must use #RRGGBB format: %s', hex);
end
