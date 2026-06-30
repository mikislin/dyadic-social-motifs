function [nodeMap, partNames, bodypartMeta] = default_sleap_node_map(configPath)
%DEFAULT_SLEAP_NODE_MAP Build the canonical SLEAP node map from config CSV.
%
% Bodypart names, display labels, MATLAB aliases, and node indices are read
% from config/sleap_bodypart_metadata.csv through the preprocessing config
% loader. This is the single supported route for paper feature extraction.

if nargin < 1 || isempty(configPath)
    cfg = load_preprocessing_pipeline_config();
else
    cfg = load_preprocessing_pipeline_config(configPath);
end

bodypartMeta = cfg.bodypoints.bodypoint_metadata;
partNames = cellstr(bodypartMeta.bodypart_name)';

nodeMap = struct();
for iNode = 1:height(bodypartMeta)
    fieldName = char(bodypartMeta.matlab_field(iNode));
    nodeMap.(fieldName) = bodypartMeta.node_index(iNode);
end

nodeMap = local_add_alias(nodeMap, 'leftEar',  'left_ear');
nodeMap = local_add_alias(nodeMap, 'rightEar', 'right_ear');
nodeMap = local_add_alias(nodeMap, 'lfPaw',    'lf_paw');
nodeMap = local_add_alias(nodeMap, 'rfPaw',    'rf_paw');
nodeMap = local_add_alias(nodeMap, 'lhPaw',    'lh_paw');
nodeMap = local_add_alias(nodeMap, 'rhPaw',    'rh_paw');
nodeMap = local_add_alias(nodeMap, 'body',     'body_pos');
nodeMap = local_add_alias(nodeMap, 'tailBase', 'tail_base');
nodeMap = local_add_alias(nodeMap, 'tailMid',  'mid_tail');
nodeMap = local_add_alias(nodeMap, 'tailTip',  'tail_tip');
nodeMap = local_add_alias(nodeMap, 'midBody',  'mid_body');

nodeMap.head = nodeMap.nose;
nodeMap.body_center = nodeMap.midBody;
nodeMap.center = nodeMap.midBody;

requiredFields = ["nose","neck","body","tailBase","tailMid","tailTip","midBody"];
missing = setdiff(requiredFields, string(fieldnames(nodeMap)));
assert(isempty(missing), 'default_sleap_node_map:MissingRequiredField', ...
    'SLEAP metadata is missing required MATLAB fields: %s', strjoin(missing, ', '));
end

function nodeMap = local_add_alias(nodeMap, sourceField, aliasField)
if isfield(nodeMap, sourceField)
    nodeMap.(aliasField) = nodeMap.(sourceField);
end
end
