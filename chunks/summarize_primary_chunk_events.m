function eventSummary = summarize_primary_chunk_events(repoRoot, sessionTable, anchorManifest, varargin)
%SUMMARIZE_PRIMARY_CHUNK_EVENTS Condition-blind within-chunk event summaries.
%
% Summaries are computed from canonical dyadic features and frame masks only.
% Biological labels are copied in anchorManifest by earlier steps but are not
% used here.

p = inputParser;
p.addParameter('contactThreshold', 0.5, @(x)isnumeric(x) && isscalar(x));
p.addParameter('stateThreshold', 0.5, @(x)isnumeric(x) && isscalar(x));
p.addParameter('turnThresholdDeg', 45, @(x)isnumeric(x) && isscalar(x) && x > 0);
p.parse(varargin{:});
P = p.Results;

if isempty(anchorManifest)
    eventSummary = table();
    return
end

required = ["feature_row_index", "start_frame", "stop_frame", "anchor_frame"];
missing = setdiff(required, string(anchorManifest.Properties.VariableNames));
assert(isempty(missing), 'summarize_primary_chunk_events:MissingColumn', ...
    'anchorManifest missing required columns: %s', strjoin(missing, ', '));

% Preallocate the output once. The previous row-wise table concatenation was
% acceptable for the 13,000-row baseline but scaled poorly for the expanded
% rare-strata bank. Numeric definitions and per-window calculations below are
% unchanged.
n = height(anchorManifest);
eventSummary = anchorManifest;
numericNames = ["chunk_duration_s","n_event_valid_frames","event_valid_fraction", ...
    "contact_dwell_fraction","contact_transition_count","first_contact_latency_s", ...
    "mutual_approach_dwell_fraction","withdrawal_dwell_fraction", ...
    "approach_withdraw_transition_count","asym_positive_dwell_fraction", ...
    "asym_negative_dwell_fraction","role_asymmetry_bias_mean", ...
    "centroid_distance_mean_mm","centroid_distance_min_mm", ...
    "centroid_distance_delta_mm","radial_speed_mean_mm_s", ...
    "radial_speed_sign_change_count","approach1_positive_fraction", ...
    "approach2_positive_fraction","role_approach_imbalance", ...
    "mutual_facing_mean","heading_abs_change_deg","heading_large_turn_count"];
for j = 1:numel(numericNames)
    eventSummary.(numericNames(j)) = nan(n, 1);
end

for sess = unique(anchorManifest.feature_row_index(:))'
    rows = find(anchorManifest.feature_row_index == sess);
    dyad = local_load_dyad(sessionTable, sess, string(repoRoot));
    frameMask = local_frame_mask(dyad);
    fps = local_fps(dyad);

    F = struct();
    names = ["centroid_dist","in_contact","mutual_approach","withdrawal", ...
        "asym_investigate","radial_speed_12","approach_speed_1", ...
        "approach_speed_2","heading_diff_deg","mutual_facing"];
    for i = 1:numel(names)
        F.(names(i)) = local_feature(dyad, names(i));
    end

    for rr = rows(:)'
        st = max(1, round(anchorManifest.start_frame(rr)));
        en = min(numel(frameMask), round(anchorManifest.stop_frame(rr)));
        if st > en
            continue
        end
        idx = st:en;
        valid = logical(frameMask(idx));
        nFrames = numel(idx);
        durationSec = nFrames ./ fps;

        contact = local_binary(F.in_contact(idx), valid, P.contactThreshold);
        mutualApproach = local_binary(F.mutual_approach(idx), valid, P.stateThreshold);
        withdrawal = local_binary(F.withdrawal(idx), valid, P.stateThreshold);
        asym = F.asym_investigate(idx);
        radial = F.radial_speed_12(idx);
        approach1 = F.approach_speed_1(idx);
        approach2 = F.approach_speed_2(idx);
        dist = F.centroid_dist(idx);
        heading = F.heading_diff_deg(idx);
        mutualFacing = F.mutual_facing(idx);

        eventSummary.chunk_duration_s(rr) = durationSec;
        eventSummary.n_event_valid_frames(rr) = nnz(valid);
        eventSummary.event_valid_fraction(rr) = nnz(valid) ./ max(nFrames, 1);
        eventSummary.contact_dwell_fraction(rr) = local_fraction_true(contact, valid);
        eventSummary.contact_transition_count(rr) = local_transition_count(contact, valid);
        eventSummary.first_contact_latency_s(rr) = local_first_latency(contact, valid, fps);
        eventSummary.mutual_approach_dwell_fraction(rr) = local_fraction_true(mutualApproach, valid);
        eventSummary.withdrawal_dwell_fraction(rr) = local_fraction_true(withdrawal, valid);
        eventSummary.approach_withdraw_transition_count(rr) = ...
            local_state_transition_count(mutualApproach, withdrawal, valid);
        eventSummary.asym_positive_dwell_fraction(rr) = local_fraction_condition(asym > 0.5, valid & isfinite(asym));
        eventSummary.asym_negative_dwell_fraction(rr) = local_fraction_condition(asym < -0.5, valid & isfinite(asym));
        eventSummary.role_asymmetry_bias_mean(rr) = mean(asym(valid & isfinite(asym)), 'omitnan');
        eventSummary.centroid_distance_mean_mm(rr) = local_mean_valid(dist, valid);
        eventSummary.centroid_distance_min_mm(rr) = local_min_valid(dist, valid);
        eventSummary.centroid_distance_delta_mm(rr) = local_late_early_delta(dist, valid);
        eventSummary.radial_speed_mean_mm_s(rr) = local_mean_valid(radial, valid);
        eventSummary.radial_speed_sign_change_count(rr) = local_sign_change_count(radial, valid);
        eventSummary.approach1_positive_fraction(rr) = local_fraction_condition(approach1 > 0, valid & isfinite(approach1));
        eventSummary.approach2_positive_fraction(rr) = local_fraction_condition(approach2 > 0, valid & isfinite(approach2));
        eventSummary.role_approach_imbalance(rr) = eventSummary.approach1_positive_fraction(rr) - ...
            eventSummary.approach2_positive_fraction(rr);
        eventSummary.mutual_facing_mean(rr) = local_mean_valid(mutualFacing, valid);
        eventSummary.heading_abs_change_deg(rr) = local_heading_abs_change(heading, valid);
        eventSummary.heading_large_turn_count(rr) = local_heading_turn_count(heading, valid, P.turnThresholdDeg);
    end
end
eventSummary.event_summary_rule = repmat("canonical_dyadic_event_features_frame_mask_only_no_labels", n, 1);
eventSummary.labels_used_for_event_summary = repmat("none", n, 1);
eventSummary.arena_used_for_event_summary = false(n, 1);
eventSummary.condition_used_for_event_summary = false(n, 1);
end

function dyad = local_load_dyad(sessionTable, rowIdx, repoRoot)
featurePath = resolve_repo_path(repoRoot, string(sessionTable.feature_output_file(rowIdx)));
S = load(featurePath, 'dyad', 'status');
assert(isfield(S, 'dyad'), 'summarize_primary_chunk_events:MissingDyad', ...
    'Missing dyad in feature file: %s', featurePath);
if isfield(S, 'status')
    assert(string(S.status) == "success", ...
        'summarize_primary_chunk_events:BadFeatureStatus', ...
        'Feature file does not report success: %s', featurePath);
end
dyad = S.dyad;
end

function fps = local_fps(dyad)
if isfield(dyad, 'fps') && ~isempty(dyad.fps)
    fps = double(dyad.fps);
elseif isfield(dyad, 'time_s')
    fps = 1 ./ median(diff(double(dyad.time_s(:))), 'omitnan');
else
    fps = 80;
end
end

function mask = local_frame_mask(dyad)
if isfield(dyad, 'frameMask') && ~isempty(dyad.frameMask)
    mask = logical(dyad.frameMask(:));
else
    T = size(dyad.X, 1);
    mask = true(T, 1);
end
end

function x = local_feature(dyad, featureName)
featureName = string(featureName);
featureChar = char(featureName);
if isfield(dyad, 'raw') && isfield(dyad.raw, featureChar)
    x = dyad.raw.(featureChar);
elseif isfield(dyad, featureChar)
    x = dyad.(featureChar);
else
    idx = find(string(dyad.featureNames(:)) == featureName, 1);
    assert(~isempty(idx), 'summarize_primary_chunk_events:MissingFeature', ...
        'Missing feature %s.', featureName);
    x = dyad.X(:, idx);
end
x = double(x(:));
end

function value = local_mean_valid(x, valid)
ok = valid(:) & isfinite(x(:));
if ~any(ok)
    value = NaN;
else
    value = mean(x(ok), 'omitnan');
end
end

function value = local_min_valid(x, valid)
ok = valid(:) & isfinite(x(:));
if ~any(ok)
    value = NaN;
else
    value = min(x(ok), [], 'omitnan');
end
end

function b = local_binary(x, valid, threshold)
b = x > threshold;
b(~valid | ~isfinite(x)) = false;
end

function frac = local_fraction_true(b, valid)
frac = local_fraction_condition(b, valid);
end

function frac = local_fraction_condition(mask, valid)
den = nnz(valid);
if den == 0
    frac = NaN;
else
    frac = nnz(mask(:) & valid(:)) ./ den;
end
end

function n = local_transition_count(b, valid)
ok = valid(:);
x = b(:);
x = x(ok);
if numel(x) < 2
    n = 0;
else
    n = sum(abs(diff(double(x))) > 0);
end
end

function n = local_state_transition_count(a, b, valid)
ok = valid(:);
state = zeros(nnz(ok), 1);
aa = a(ok);
bb = b(ok);
state(aa) = 1;
state(bb) = -1;
if numel(state) < 2
    n = 0;
else
    n = sum(abs(diff(state)) > 1);
end
end

function latency = local_first_latency(b, valid, fps)
idx = find(b(:) & valid(:), 1, 'first');
if isempty(idx)
    latency = NaN;
else
    latency = (idx - 1) ./ fps;
end
end

function n = local_sign_change_count(x, valid)
ok = valid(:) & isfinite(x(:));
xx = x(ok);
xx(abs(xx) < eps) = 0;
sgn = sign(xx);
sgn = sgn(sgn ~= 0);
if numel(sgn) < 2
    n = 0;
else
    n = sum(abs(diff(sgn)) > 0);
end
end

function delta = local_late_early_delta(x, valid)
ok = valid(:) & isfinite(x(:));
if nnz(ok) < 2
    delta = NaN;
    return
end
xx = x(:);
n = numel(xx);
third = max(1, floor(n / 3));
early = false(n, 1);
late = false(n, 1);
early(1:third) = true;
late((n - third + 1):n) = true;
delta = mean(xx(late & ok), 'omitnan') - mean(xx(early & ok), 'omitnan');
end

function value = local_heading_abs_change(heading, valid)
ok = valid(:) & isfinite(heading(:));
if nnz(ok) < 2
    value = NaN;
    return
end
h = unwrap(deg2rad(heading(ok)));
value = rad2deg(sum(abs(diff(h)), 'omitnan'));
end

function n = local_heading_turn_count(heading, valid, thresholdDeg)
ok = valid(:) & isfinite(heading(:));
if nnz(ok) < 2
    n = 0;
    return
end
h = unwrap(deg2rad(heading(ok)));
d = abs(rad2deg(diff(h)));
n = nnz(d >= thresholdDeg);
end
