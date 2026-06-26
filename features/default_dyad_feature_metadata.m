function [featureNames, featureMeta] = default_dyad_feature_metadata()
%DEFAULT_DYAD_FEATURE_METADATA Canonical dyadic feature names and metadata.
%   [featureNames, featureMeta] = default_dyad_feature_metadata()
%   returns the ordered feature list used by compute_dyad_features and by
%   synthetic tests for downstream chunking, embedding, and scale analyses.

featureNames = {
    'centroid_dist'
    'body2body_dist'
    'head2head_dist'
    'tailbase2tailbase_dist'
    'nose1_to_body2_dist'
    'nose2_to_body1_dist'
    'nose1_to_tail2_dist'
    'nose2_to_tail1_dist'
    'partner_long_1'
    'partner_lat_1'
    'partner_long_2'
    'partner_lat_2'
    'facing_1_to_2'
    'facing_2_to_1'
    'mutual_facing'
    'heading_diff_deg'
    'cos_head_alignment'
    'radial_speed_12'
    'tangential_speed_12'
    'approach_speed_1'
    'approach_speed_2'
    'speed_alignment'
    'accel_alignment'
    'nose_bearing_1_deg'
    'nose_bearing_2_deg'
    'body_bearing_1_deg'
    'body_bearing_2_deg'
    'in_contact'
    'head_close'
    'body_close'
    'close_pair'
    'mutual_approach'
    'withdrawal'
    'asym_investigate'
    }';

families = {
    'distance'
    'distance'
    'distance'
    'distance'
    'distance'
    'distance'
    'distance'
    'distance'
    'egocentric'
    'egocentric'
    'egocentric'
    'egocentric'
    'orientation'
    'orientation'
    'orientation'
    'orientation'
    'orientation'
    'kinematics'
    'kinematics'
    'kinematics'
    'kinematics'
    'coupling'
    'coupling'
    'orientation'
    'orientation'
    'orientation'
    'orientation'
    'contact'
    'contact'
    'contact'
    'contact'
    'interaction_logic'
    'interaction_logic'
    'interaction_logic'
    }';

isCircular = contains(featureNames, 'bearing_')' | contains(featureNames, 'heading_diff')';
isBoolean = ismember(featureNames, ...
    {'in_contact','head_close','body_close','close_pair','mutual_approach','withdrawal'});
isDirected = contains(featureNames, '_1') | contains(featureNames, '_2') | ...
    contains(featureNames, '1_to_') | contains(featureNames, '2_to_');

transformHint = repmat("none", numel(featureNames), 1);
for k = 1:numel(featureNames)
    if contains(featureNames{k}, 'dist')
        transformHint(k) = "log1p";
    elseif isBoolean(k)
        transformHint(k) = "binary";
    elseif isCircular(k)
        transformHint(k) = "circular";
    else
        transformHint(k) = "zscore";
    end
end

featureMeta = table(featureNames(:), families(:), isDirected(:), isCircular(:), ...
    isBoolean(:), transformHint, ...
    'VariableNames', {'Name','Family','IsDirected','IsCircular','IsBoolean','TransformHint'});
end
