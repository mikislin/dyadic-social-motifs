function [featureNames, featureMeta] = default_dyad_feature_metadata()
%DEFAULT_DYAD_FEATURE_METADATA Canonical dyadic feature dictionary.
%   [featureNames, featureMeta] = default_dyad_feature_metadata() returns the
%   ordered feature list used by compute_dyad_features, chunking, embedding,
%   and reviewer-facing feature audit tables.

rows = {
    'centroid_dist',          'distance',          'Inter-animal distance',      'mm',           'millimeters',       false, 'symmetric_pair',        'none',     'arena_xy',                  false, false, 'log1p',    '[0, Inf)',  'Distance between trunk centroids of the two animals.'
    'body2body_dist',         'distance',          'Inter-animal distance',      'mm',           'millimeters',       false, 'symmetric_pair',        'none',     'arena_xy',                  false, false, 'log1p',    '[0, Inf)',  'Distance between SLEAP body points.'
    'head2head_dist',         'distance',          'Inter-animal distance',      'mm',           'millimeters',       false, 'symmetric_pair',        'none',     'arena_xy',                  false, false, 'log1p',    '[0, Inf)',  'Nose-to-nose distance, used as a head proximity proxy.'
    'tailbase2tailbase_dist', 'distance',          'Inter-animal distance',      'mm',           'millimeters',       false, 'symmetric_pair',        'none',     'arena_xy',                  false, false, 'log1p',    '[0, Inf)',  'Distance between tail-base points.'
    'nose1_to_body2_dist',    'distance',          'Directed target distance',   'mm',           'millimeters',       true,  'animal1_to_animal2',    'animal_1', 'arena_xy',                  false, false, 'log1p',    '[0, Inf)',  'Animal 1 nose distance to animal 2 body.'
    'nose2_to_body1_dist',    'distance',          'Directed target distance',   'mm',           'millimeters',       true,  'animal2_to_animal1',    'animal_2', 'arena_xy',                  false, false, 'log1p',    '[0, Inf)',  'Animal 2 nose distance to animal 1 body.'
    'nose1_to_tail2_dist',    'distance',          'Directed target distance',   'mm',           'millimeters',       true,  'animal1_to_animal2',    'animal_1', 'arena_xy',                  false, false, 'log1p',    '[0, Inf)',  'Animal 1 nose distance to animal 2 tail base.'
    'nose2_to_tail1_dist',    'distance',          'Directed target distance',   'mm',           'millimeters',       true,  'animal2_to_animal1',    'animal_2', 'arena_xy',                  false, false, 'log1p',    '[0, Inf)',  'Animal 2 nose distance to animal 1 tail base.'
    'partner_long_1',         'egocentric',        'Egocentric partner position','mm',           'millimeters',       true,  'animal2_relative_to_1', 'animal_1', 'animal1_tailbase_to_nose', false, false, 'zscore',   'signed',     'Animal 2 centroid position along animal 1 head-body axis; positive is in front of animal 1.'
    'partner_lat_1',          'egocentric',        'Egocentric partner position','mm',           'millimeters',       true,  'animal2_relative_to_1', 'animal_1', 'animal1_tailbase_to_nose', false, false, 'zscore',   'signed',     'Animal 2 centroid position lateral to animal 1 heading; positive is +90 degrees in track coordinates.'
    'partner_long_2',         'egocentric',        'Egocentric partner position','mm',           'millimeters',       true,  'animal1_relative_to_2', 'animal_2', 'animal2_tailbase_to_nose', false, false, 'zscore',   'signed',     'Animal 1 centroid position along animal 2 head-body axis; positive is in front of animal 2.'
    'partner_lat_2',          'egocentric',        'Egocentric partner position','mm',           'millimeters',       true,  'animal1_relative_to_2', 'animal_2', 'animal2_tailbase_to_nose', false, false, 'zscore',   'signed',     'Animal 1 centroid position lateral to animal 2 heading; positive is +90 degrees in track coordinates.'
    'facing_1_to_2',          'orientation',       'Facing score',              'unitless',     'cosine score',      true,  'animal1_to_animal2',    'animal_1', 'animal1_tailbase_to_nose', false, false, 'zscore',   '[-1, 1]',    'Cosine of animal 1 heading toward animal 2 centroid; larger means more directly facing.'
    'facing_2_to_1',          'orientation',       'Facing score',              'unitless',     'cosine score',      true,  'animal2_to_animal1',    'animal_2', 'animal2_tailbase_to_nose', false, false, 'zscore',   '[-1, 1]',    'Cosine of animal 2 heading toward animal 1 centroid; larger means more directly facing.'
    'mutual_facing',          'orientation',       'Facing score',              'unitless',     'cosine score',      false, 'symmetric_minimum',     'both',     'animal_heading_axes',       false, false, 'zscore',   '[-1, 1]',    'Minimum of the two directed facing scores; high only when both animals face each other.'
    'heading_diff_deg',       'orientation',       'Relative orientation',       'degrees',      'degrees',           false, 'signed_pair_angle',     'both',     'animal_heading_axes',       true,  false, 'circular', '[-180, 180]', 'Signed angle from animal 1 heading to animal 2 heading.'
    'cos_head_alignment',     'orientation',       'Relative orientation',       'unitless',     'cosine score',      false, 'symmetric_pair',        'both',     'animal_heading_axes',       false, false, 'zscore',   '[-1, 1]',    'Cosine alignment of the two head-body axes; 1 is parallel, -1 is antiparallel.'
    'radial_speed_12',        'kinematics',        'Relative velocity',          'mm_per_s',     'millimeters/second',true,  'animal2_relative_to_1', 'both',     'centroid_line_1_to_2',      false, false, 'zscore',   'signed',     'Relative speed of animal 2 versus animal 1 along the centroid line; positive means separation increasing.'
    'tangential_speed_12',    'kinematics',        'Relative velocity',          'mm_per_s',     'millimeters/second',true,  'animal2_relative_to_1', 'both',     'centroid_line_1_to_2',      false, false, 'zscore',   'signed',     'Relative speed perpendicular to the centroid line; sign follows +90 degrees from animal 1-to-2 direction.'
    'approach_speed_1',       'kinematics',        'Approach velocity',          'mm_per_s',     'millimeters/second',true,  'animal1_to_animal2',    'animal_1', 'centroid_line_1_to_2',      false, false, 'zscore',   'signed',     'Animal 1 approach speed toward animal 2; positive means animal 1 closes distance.'
    'approach_speed_2',       'kinematics',        'Approach velocity',          'mm_per_s',     'millimeters/second',true,  'animal2_to_animal1',    'animal_2', 'centroid_line_2_to_1',      false, false, 'zscore',   'signed',     'Animal 2 approach speed toward animal 1; positive means animal 2 closes distance.'
    'speed_alignment',        'coupling',          'Movement coupling',          'unitless',     'cosine score',      false, 'symmetric_pair',        'both',     'centroid_velocity',         false, false, 'zscore',   '[-1, 1]',    'Cosine alignment of centroid velocity vectors.'
    'accel_alignment',        'coupling',          'Movement coupling',          'unitless',     'cosine score',      false, 'symmetric_pair',        'both',     'centroid_acceleration',     false, false, 'zscore',   '[-1, 1]',    'Cosine alignment of centroid acceleration vectors.'
    'nose_bearing_1_deg',     'orientation',       'Egocentric bearing',         'degrees',      'degrees',           true,  'animal1_to_animal2',    'animal_1', 'animal1_tailbase_to_nose', true,  false, 'circular', '[-180, 180]', 'Signed egocentric bearing from animal 1 heading to animal 2 nose.'
    'nose_bearing_2_deg',     'orientation',       'Egocentric bearing',         'degrees',      'degrees',           true,  'animal2_to_animal1',    'animal_2', 'animal2_tailbase_to_nose', true,  false, 'circular', '[-180, 180]', 'Signed egocentric bearing from animal 2 heading to animal 1 nose.'
    'body_bearing_1_deg',     'orientation',       'Egocentric bearing',         'degrees',      'degrees',           true,  'animal1_to_animal2',    'animal_1', 'animal1_tailbase_to_nose', true,  false, 'circular', '[-180, 180]', 'Signed egocentric bearing from animal 1 heading to animal 2 body.'
    'body_bearing_2_deg',     'orientation',       'Egocentric bearing',         'degrees',      'degrees',           true,  'animal2_to_animal1',    'animal_2', 'animal2_tailbase_to_nose', true,  false, 'circular', '[-180, 180]', 'Signed egocentric bearing from animal 2 heading to animal 1 body.'
    'in_contact',             'contact',           'Proximity state',            'boolean',      '0/1 indicator',     false, 'symmetric_pair',        'both',     'centroid_distance',         false, true,  'binary',   '{0, 1}',     'Whether centroid distance is at or below the contact threshold.'
    'head_close',             'contact',           'Proximity state',            'boolean',      '0/1 indicator',     false, 'symmetric_pair',        'both',     'nose_distance',             false, true,  'binary',   '{0, 1}',     'Whether nose-to-nose distance is at or below the contact threshold.'
    'body_close',             'contact',           'Proximity state',            'boolean',      '0/1 indicator',     false, 'symmetric_pair',        'both',     'body_distance',             false, true,  'binary',   '{0, 1}',     'Whether body-to-body distance is at or below the contact threshold.'
    'close_pair',             'contact',           'Proximity state',            'boolean',      '0/1 indicator',     false, 'symmetric_pair',        'both',     'centroid_distance',         false, true,  'binary',   '{0, 1}',     'Whether centroid distance is at or below the broader close-pair threshold.'
    'mutual_approach',        'interaction_logic', 'Interaction state',          'boolean',      '0/1 indicator',     false, 'symmetric_pair',        'both',     'approach_speeds',           false, true,  'binary',   '{0, 1}',     'Whether both directed approach speeds are positive.'
    'withdrawal',             'interaction_logic', 'Interaction state',          'boolean',      '0/1 indicator',     false, 'symmetric_pair',        'both',     'approach_speeds',           false, true,  'binary',   '{0, 1}',     'Whether both directed approach speeds are negative.'
    'asym_investigate',       'interaction_logic', 'Asymmetry index',            'signed_index', '-1/0/1 index',      true,  'antisymmetric_pair',    'both',     'nose_to_body_distances',    false, false, 'zscore',   '{-1,0,1}',   'Signed asymmetry in nose-to-partner-body distance; +1 means animal 1 is closer to animal 2 body than animal 2 is to animal 1 body.'
    };

featureMeta = table();
featureMeta.Name = string(rows(:,1));
featureMeta.Family = string(rows(:,2));
featureMeta.FamilyLabel = string(rows(:,3));
featureMeta.FeatureFamilyRole = repmat( ...
    "computational_measurement_group_not_biological_motif_family", ...
    size(featureMeta.Name));
featureMeta.FamilyDefinition = local_family_definition(featureMeta.Family);
featureMeta.Unit = string(rows(:,4));
featureMeta.UnitLabel = string(rows(:,5));
featureMeta.IsDirected = cell2mat(rows(:,6));
featureMeta.Directionality = string(rows(:,7));
featureMeta.FocalAnimal = string(rows(:,8));
featureMeta.ReferenceFrame = string(rows(:,9));
featureMeta.IsCircular = cell2mat(rows(:,10));
featureMeta.IsBoolean = cell2mat(rows(:,11));
featureMeta.TransformHint = string(rows(:,12));
featureMeta.ExpectedRange = string(rows(:,13));
featureMeta.BiologicalInterpretation = string(rows(:,14));
featureMeta.ClusteringCandidate = true(height(featureMeta), 1);
featureMeta.FeatureLayerRole = repmat("condition_blind_frame_feature_for_chunking", height(featureMeta), 1);
featureMeta.WindowSummaryStats = repmat("mean;std;min;max;delta;slope", height(featureMeta), 1);
featureMeta.WindowSummaryStats(featureMeta.IsCircular) = "circmean;resultant";
featureMeta.WindowSummaryStats(featureMeta.IsBoolean) = "occupancy;onset_rate";
featureMeta.FeatureIndex = (1:height(featureMeta))';
featureMeta = movevars(featureMeta, 'FeatureIndex', 'Before', 'Name');

assert(numel(unique(featureMeta.Name)) == height(featureMeta), ...
    'default_dyad_feature_metadata:DuplicateName', ...
    'Feature names must be unique.');
assert(all(strlength(featureMeta.Unit) > 0), ...
    'default_dyad_feature_metadata:MissingUnit', ...
    'Every feature must have an explicit unit.');
assert(all(strlength(featureMeta.BiologicalInterpretation) > 0), ...
    'default_dyad_feature_metadata:MissingInterpretation', ...
    'Every feature must have a biological interpretation.');

featureNames = cellstr(featureMeta.Name)';
end

function definitions = local_family_definition(families)
families = string(families);
definitions = strings(size(families));
for i = 1:numel(families)
    switch families(i)
        case "distance"
            definitions(i) = "Euclidean inter-animal distances from canonical SLEAP bodypoints.";
        case "egocentric"
            definitions(i) = "Partner position expressed in a focal animal-centered heading frame.";
        case "orientation"
            definitions(i) = "Heading, facing, alignment, and bearing relationships between animals.";
        case "kinematics"
            definitions(i) = "Relative velocity components and directed approach speeds.";
        case "coupling"
            definitions(i) = "Similarity of movement direction or acceleration between animals.";
        case "contact"
            definitions(i) = "Thresholded proximity states using configured millimeter cutoffs.";
        case "interaction_logic"
            definitions(i) = "Composite signed or Boolean states derived from directed dyadic features.";
        otherwise
            definitions(i) = "Unclassified computational feature group.";
    end
end
end
