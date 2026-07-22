function [Definition, Membership, Seed] = define_condition_blind_rare_strata(graphRoot, varargin)
%DEFINE_CONDITION_BLIND_RARE_STRATA Define post-fit label-free graph strata.
%
% Strata use only run_08 graph geometry, run_07 row identifiers and
% split-half stability, and run_06 condition-blind event summaries. Arena and
% experimental labels are never selected from source CSV files and never
% enter thresholds or membership rules.

p = inputParser;
p.addParameter('outputDir', "", @(x)ischar(x) || isstring(x));
p.addParameter('embeddingRoot', "", @(x)ischar(x) || isstring(x));
p.addParameter('chunkRoot', "", @(x)ischar(x) || isstring(x));
p.addParameter('eventSummaryFile', "primary_chunk_event_summary_audit.csv", @(x)ischar(x) || isstring(x));
p.addParameter('nodeManifest', table(), @istable);
p.addParameter('degreeAudit', table(), @istable);
p.addParameter('rowManifest', table(), @istable);
p.addParameter('eventSummary', table(), @istable);
p.addParameter('stabilityAudit', table(), @istable);
p.addParameter('lowDensityQuantile', 0.90, @(x)isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('peripheryInDegreeQuantile', 0.10, @(x)isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('highRadialQuantile', 0.90, @(x)isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('undercoveredCellQuantile', 0.10, @(x)isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('shortInstabilityScaleQuantile', 0.50, @(x)isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
p.addParameter('shortInstabilityRadiusQuantile', 0.75, @(x)isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('shortInstabilityTransitionQuantile', 0.75, @(x)isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('writeOutputs', true, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});
P = p.Results;

graphRoot = string(graphRoot);
outputDir = string(P.outputDir);
if outputDir == ""
    outputDir = graphRoot;
end
embeddingRoot = string(P.embeddingRoot);
chunkRoot = string(P.chunkRoot);

Node = local_input_table(P.nodeManifest, graphRoot, 'graph_node_manifest.csv', ...
    ["graph_node_id","embedding_row_id","scale_index","chunk_sec", ...
    "hierarchical_role","session_index","raw_index","session_id", ...
    "anchor_frame","anchor_time_s","primary_anchor_global_id", ...
    "expanded_anchor_global_id"]);
Degree = local_input_table(P.degreeAudit, graphRoot, 'graph_degree_audit.csv', ...
    ["graph_node_id","embedding_row_id","scale_index","chunk_sec", ...
    "knn_in_degree","undirected_degree","mutual_neighbor_fraction", ...
    "knn_radius","labels_used_for_degree_audit", ...
    "arena_used_for_degree_audit","condition_used_for_degree_audit"]);
Rows = local_input_table(P.rowManifest, embeddingRoot, 'embedding_row_manifest.csv', ...
    ["embedding_row_id","scale_index","primary_scale_specific_anchor_id", ...
    "primary_anchor_global_id","expanded_anchor_global_id","feature_row_index", ...
    "session_index","raw_index","session_id","anchor_frame","anchor_time_s"]);
Event = local_input_table(P.eventSummary, chunkRoot, string(P.eventSummaryFile), ...
    ["scale_index","primary_anchor_global_id","expanded_anchor_global_id", ...
    "event_valid_fraction","contact_dwell_fraction","contact_transition_count", ...
    "approach_withdraw_transition_count","heading_large_turn_count", ...
    "radial_speed_mean_mm_s","labels_used_for_event_summary", ...
    "arena_used_for_event_summary","condition_used_for_event_summary"]);
Stability = local_input_table(P.stabilityAudit, embeddingRoot, 'embedding_stability_audit.csv', ...
    ["scale_index","chunk_sec","hierarchical_role","median_subspace_similarity", ...
    "labels_used_for_stability","arena_used_for_stability","condition_used_for_stability"]);

local_validate_inputs(Node, Degree, Rows, Event, Stability);
[Node, Degree] = local_align_by_id(Node, Degree, 'embedding_row_id');
[Node, Rows] = local_align_by_id(Node, Rows, 'embedding_row_id');
Event = local_align_events(Node, Event);

n = height(Node);
strata = ["low_density_tail","graph_periphery","contact_present", ...
    "contact_transition","large_heading_turn","high_radial_speed", ...
    "undercovered_scale_session","short_motif_instability"];
nStrata = numel(strata);
flag = false(n, nStrata);
score = nan(n, nStrata);
threshold = nan(n, nStrata);
threshold2 = nan(n, nStrata);

rules = [ ...
    "knn_radius_ge_within_scale_quantile"; ...
    "knn_in_degree_le_within_scale_quantile"; ...
    "contact_dwell_fraction_gt_zero"; ...
    "contact_transition_count_gt_zero"; ...
    "heading_large_turn_count_gt_zero"; ...
    "radial_speed_mean_ge_within_scale_quantile"; ...
    "scale_session_node_count_le_within_scale_quantile"; ...
    "motif_scale_low_stability_and_high_radius_and_transition_burden"];
sources = [ ...
    "run08_graph_degree_audit.knn_radius"; ...
    "run08_graph_degree_audit.knn_in_degree"; ...
    "run06_event_summary.contact_dwell_fraction"; ...
    "run06_event_summary.contact_transition_count"; ...
    "run06_event_summary.heading_large_turn_count"; ...
    "run06_event_summary.radial_speed_mean_mm_s"; ...
    "run08_graph_node_manifest.scale_session_support"; ...
    "run07_embedding_stability_plus_run08_radius_plus_run06_transition_burden"];

radius = double(Degree.knn_radius);
inDegree = double(Degree.knn_in_degree);
contact = double(Event.contact_dwell_fraction);
contactTransition = double(Event.contact_transition_count);
headingTurn = double(Event.heading_large_turn_count);
radial = double(Event.radial_speed_mean_mm_s);
approachTransition = double(Event.approach_withdraw_transition_count);
transitionBurden = contactTransition + headingTurn + approachTransition;
eventValid = double(Event.event_valid_fraction) > 0;

stabilityByNode = nan(n, 1);
[tfStab, locStab] = ismember(double(Node.scale_index), double(Stability.scale_index));
stabilityByNode(tfStab) = double(Stability.median_subspace_similarity(locStab(tfStab)));

cellCount = zeros(n, 1);
cellKey = local_cell_key(Node.scale_index, Node.session_index);
[uCell, ~, cellGroup] = unique(cellKey, 'stable'); %#ok<ASGLU>
counts = accumarray(cellGroup, 1);
cellCount = counts(cellGroup);

Definition = table();
scales = unique(double(Node.scale_index), 'stable')';
motifScale = local_motif_scale_mask(Node, Stability);
motifStability = unique(stabilityByNode(motifScale & isfinite(stabilityByNode)));
if isempty(motifStability)
    motifStabilityThreshold = -Inf;
else
    motifStabilityThreshold = quantile(motifStability, P.shortInstabilityScaleQuantile);
end

for s = scales
    idx = double(Node.scale_index) == s;
    chunkSec = double(Node.chunk_sec(find(idx, 1)));
    radiusThreshold = local_finite_quantile(radius(idx), P.lowDensityQuantile, Inf);
    peripheryThreshold = local_finite_quantile(inDegree(idx), P.peripheryInDegreeQuantile, -Inf);
    radialThreshold = local_finite_quantile(radial(idx & eventValid), P.highRadialQuantile, Inf);
    countThreshold = local_finite_quantile(double(unique(cellCount(idx))), P.undercoveredCellQuantile, -Inf);
    instabilityRadiusThreshold = local_finite_quantile(radius(idx), P.shortInstabilityRadiusQuantile, Inf);
    transitionThreshold = local_finite_quantile(transitionBurden(idx & eventValid), ...
        P.shortInstabilityTransitionQuantile, Inf);

    flag(idx, 1) = isfinite(radius(idx)) & radius(idx) >= radiusThreshold;
    score(idx, 1) = radius(idx) ./ max(radiusThreshold, eps);
    threshold(idx, 1) = radiusThreshold;

    flag(idx, 2) = isfinite(inDegree(idx)) & inDegree(idx) <= peripheryThreshold;
    score(idx, 2) = (peripheryThreshold + 1) ./ (inDegree(idx) + 1);
    threshold(idx, 2) = peripheryThreshold;

    flag(idx, 3) = eventValid(idx) & isfinite(contact(idx)) & contact(idx) > 0;
    score(idx, 3) = contact(idx);
    threshold(idx, 3) = 0;

    flag(idx, 4) = eventValid(idx) & isfinite(contactTransition(idx)) & contactTransition(idx) > 0;
    score(idx, 4) = contactTransition(idx);
    threshold(idx, 4) = 0;

    flag(idx, 5) = eventValid(idx) & isfinite(headingTurn(idx)) & headingTurn(idx) > 0;
    score(idx, 5) = headingTurn(idx);
    threshold(idx, 5) = 0;

    flag(idx, 6) = eventValid(idx) & isfinite(radial(idx)) & radial(idx) >= radialThreshold;
    score(idx, 6) = radial(idx) ./ max(abs(radialThreshold), eps);
    threshold(idx, 6) = radialThreshold;

    flag(idx, 7) = isfinite(cellCount(idx)) & cellCount(idx) <= countThreshold;
    score(idx, 7) = (countThreshold + 1) ./ (cellCount(idx) + 1);
    threshold(idx, 7) = countThreshold;

    shortFlag = motifScale(idx) & isfinite(stabilityByNode(idx)) & ...
        stabilityByNode(idx) <= motifStabilityThreshold & ...
        isfinite(radius(idx)) & radius(idx) >= instabilityRadiusThreshold & ...
        isfinite(transitionBurden(idx)) & transitionBurden(idx) >= transitionThreshold;
    flag(idx, 8) = shortFlag;
    stabilityDeficit = max(motifStabilityThreshold - stabilityByNode(idx), 0) + eps;
    score(idx, 8) = stabilityDeficit .* ...
        (radius(idx) ./ max(instabilityRadiusThreshold, eps)) .* ...
        (transitionBurden(idx) ./ max(transitionThreshold, 1));
    threshold(idx, 8) = motifStabilityThreshold;
    threshold2(idx, 8) = instabilityRadiusThreshold;

    scaleThreshold = [radiusThreshold; peripheryThreshold; 0; 0; 0; ...
        radialThreshold; countThreshold; motifStabilityThreshold];
    secondary = [NaN; NaN; NaN; NaN; NaN; NaN; NaN; instabilityRadiusThreshold];
    tertiary = [NaN; NaN; NaN; NaN; NaN; NaN; NaN; transitionThreshold];
    for j = 1:nStrata
        one = table();
        one.rare_stratum_id = strata(j);
        one.rare_stratum_rule = rules(j);
        one.rare_stratum_source = sources(j);
        one.scale_index = s;
        one.chunk_sec = chunkSec;
        one.rare_stratum_threshold = scaleThreshold(j);
        one.rare_stratum_secondary_threshold = secondary(j);
        one.rare_stratum_tertiary_threshold = tertiary(j);
        one.threshold_scope = "within_scale_condition_blind_baseline";
        one.n_source_nodes = nnz(idx);
        one.n_member_nodes = nnz(flag(idx, j));
        one.member_fraction = nnz(flag(idx, j)) ./ max(nnz(idx), 1);
        one.labels_used_for_rare_strata = "none";
        one.arena_used_for_rare_strata = false;
        one.condition_used_for_rare_strata = false;
        Definition = [Definition; one]; %#ok<AGROW>
    end
end

Base = table();
Base.graph_node_id = double(Node.graph_node_id);
Base.embedding_row_id = double(Node.embedding_row_id);
Base.scale_index = double(Node.scale_index);
Base.chunk_sec = double(Node.chunk_sec);
Base.source_primary_scale_specific_anchor_id = local_numeric_column(Rows, ...
    'primary_scale_specific_anchor_id', n, NaN);
Base.primary_anchor_global_id = local_numeric_column(Rows, 'primary_anchor_global_id', n, NaN);
Base.expanded_anchor_global_id = local_numeric_column(Rows, 'expanded_anchor_global_id', n, NaN);
Base.feature_row_index = local_numeric_column(Rows, 'feature_row_index', n, NaN);
Base.session_index = double(Node.session_index);
Base.session_raw_index = double(Node.raw_index);
Base.source_session_id = string(Node.session_id);
Base.anchor_frame = double(Node.anchor_frame);
Base.anchor_time_s = double(Node.anchor_time_s);

Membership = table();
for j = 1:nStrata
    idx = find(flag(:, j));
    if isempty(idx)
        continue
    end
    one = Base(idx, :);
    one.rare_stratum_id = repmat(strata(j), numel(idx), 1);
    one.rare_stratum_rule = repmat(rules(j), numel(idx), 1);
    one.rare_stratum_source = repmat(sources(j), numel(idx), 1);
    one.rare_stratum_threshold = threshold(idx, j);
    one.rare_stratum_secondary_threshold = threshold2(idx, j);
    one.rare_stratum_score = score(idx, j);
    one.labels_used_for_rare_strata = repmat("none", numel(idx), 1);
    one.arena_used_for_rare_strata = false(numel(idx), 1);
    one.condition_used_for_rare_strata = false(numel(idx), 1);
    Membership = [Membership; one]; %#ok<AGROW>
end
Membership = sortrows(Membership, {'scale_index','session_raw_index','anchor_frame','rare_stratum_id'});
Membership.rare_strata_membership_id = (1:height(Membership))';
Membership = movevars(Membership, 'rare_strata_membership_id', 'Before', 1);

Seed = Membership;
Seed.seed_anchor_role = repmat("baseline_condition_blind_rare_stratum_seed", height(Seed), 1);
Seed.baseline_graph_root = repmat(graphRoot, height(Seed), 1);
Seed.baseline_embedding_root = repmat(embeddingRoot, height(Seed), 1);
Seed.baseline_chunk_root = repmat(chunkRoot, height(Seed), 1);

if logical(P.writeOutputs)
    local_ensure_dir(outputDir);
    writetable(Definition, fullfile(outputDir, 'rare_strata_definition.csv'));
    writetable(Membership, fullfile(outputDir, 'rare_strata_node_membership.csv'));
    writetable(Seed, fullfile(outputDir, 'rare_strata_seed_anchor_manifest.csv'));
end
end

function T = local_input_table(provided, rootDir, fileName, wanted)
if ~isempty(provided)
    T = provided;
    keep = intersect(wanted, string(T.Properties.VariableNames), 'stable');
    T = T(:, cellstr(keep));
    return
end
assert(strlength(string(rootDir)) > 0, ...
    'define_condition_blind_rare_strata:MissingInputRoot', ...
    'An input table or its source root is required for %s.', fileName);
pathText = fullfile(rootDir, char(fileName));
assert(isfile(pathText), 'define_condition_blind_rare_strata:MissingInputFile', ...
    'Missing required rare-strata input: %s', pathText);
opts = detectImportOptions(pathText, 'FileType', 'text', 'Delimiter', ',', 'TextType', 'string');
keep = intersect(wanted, string(opts.VariableNames), 'stable');
opts.SelectedVariableNames = cellstr(keep);
T = readtable(pathText, opts);
end

function local_validate_inputs(Node, Degree, Rows, Event, Stability)
local_require(Node, ["graph_node_id","embedding_row_id","scale_index","chunk_sec", ...
    "session_index","raw_index","session_id","anchor_frame","anchor_time_s"]);
local_require(Degree, ["embedding_row_id","knn_in_degree","knn_radius"]);
local_require(Rows, ["embedding_row_id"]);
local_require(Event, ["event_valid_fraction","contact_dwell_fraction", ...
    "contact_transition_count","approach_withdraw_transition_count", ...
    "heading_large_turn_count","radial_speed_mean_mm_s"]);
local_require(Stability, ["scale_index","median_subspace_similarity"]);
local_assert_none(Degree, 'labels_used_for_degree_audit');
local_assert_false(Degree, 'arena_used_for_degree_audit');
local_assert_false(Degree, 'condition_used_for_degree_audit');
local_assert_none(Event, 'labels_used_for_event_summary');
local_assert_false(Event, 'arena_used_for_event_summary');
local_assert_false(Event, 'condition_used_for_event_summary');
local_assert_none(Stability, 'labels_used_for_stability');
local_assert_false(Stability, 'arena_used_for_stability');
local_assert_false(Stability, 'condition_used_for_stability');
end

function local_require(T, names)
missing = setdiff(names, string(T.Properties.VariableNames));
assert(isempty(missing), 'define_condition_blind_rare_strata:MissingColumn', ...
    'Required columns are missing: %s', strjoin(missing, ', '));
end

function local_assert_none(T, name)
if ismember(name, T.Properties.VariableNames)
    assert(all(string(T.(name)) == "none"), ...
        'define_condition_blind_rare_strata:LabelLeak', '%s must equal none.', name);
end
end

function local_assert_false(T, name)
if ismember(name, T.Properties.VariableNames)
    assert(~any(logical(T.(name))), ...
        'define_condition_blind_rare_strata:LabelLeak', '%s must be false.', name);
end
end

function [Left, Right] = local_align_by_id(Left, Right, keyName)
[tf, loc] = ismember(double(Left.(keyName)), double(Right.(keyName)));
assert(all(tf), 'define_condition_blind_rare_strata:AlignmentFailure', ...
    'Rows could not be aligned by %s.', keyName);
Right = Right(loc, :);
end

function Event = local_align_events(Node, Event)
if ismember('expanded_anchor_global_id', Node.Properties.VariableNames) && ...
        ismember('expanded_anchor_global_id', Event.Properties.VariableNames) && ...
        all(isfinite(double(Node.expanded_anchor_global_id)))
    key = 'expanded_anchor_global_id';
elseif ismember('primary_anchor_global_id', Node.Properties.VariableNames) && ...
        ismember('primary_anchor_global_id', Event.Properties.VariableNames)
    key = 'primary_anchor_global_id';
else
    error('define_condition_blind_rare_strata:MissingEventJoinKey', ...
        'No shared condition-blind anchor identifier is available for event summaries.');
end
[tf, loc] = ismember(double(Node.(key)), double(Event.(key)));
assert(all(tf), 'define_condition_blind_rare_strata:MissingEventRows', ...
    'Event summary rows could not be matched by %s.', key);
Event = Event(loc, :);
end

function tf = local_motif_scale_mask(Node, Stability)
tf = false(height(Node), 1);
if ismember('hierarchical_role', Node.Properties.VariableNames)
    tf = string(Node.hierarchical_role) == "motif";
elseif ismember('hierarchical_role', Stability.Properties.VariableNames)
    [ok, loc] = ismember(double(Node.scale_index), double(Stability.scale_index));
    tf(ok) = string(Stability.hierarchical_role(loc(ok))) == "motif";
else
    tf = double(Node.chunk_sec) >= 0.8 & double(Node.chunk_sec) < 2.5;
end
end

function q = local_finite_quantile(x, p, fallback)
x = double(x(:));
x = x(isfinite(x));
if isempty(x)
    q = fallback;
else
    q = quantile(x, p);
end
end

function key = local_cell_key(scaleIndex, sessionIndex)
key = string(round(double(scaleIndex(:)))) + "_" + string(round(double(sessionIndex(:))));
end

function x = local_numeric_column(T, name, n, defaultValue)
if ismember(name, T.Properties.VariableNames)
    x = double(T.(name));
else
    x = repmat(defaultValue, n, 1);
end
end

function local_ensure_dir(pathText)
if ~exist(pathText, 'dir')
    mkdir(pathText);
end
end
