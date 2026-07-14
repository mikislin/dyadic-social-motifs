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

eventSummary = table();
if isempty(anchorManifest)
    return
end

required = ["feature_row_index", "start_frame", "stop_frame", "anchor_frame"];
missing = setdiff(required, string(anchorManifest.Properties.VariableNames));
assert(isempty(missing), 'summarize_primary_chunk_events:MissingColumn', ...
    'anchorManifest missing required columns: %s', strjoin(missing, ', '));

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

        one = anchorManifest(rr, :);
        one.chunk_duration_s = durationSec;
        one.n_event_valid_frames = nnz(valid);
        one.event_valid_fraction = nnz(valid) ./ max(nFrames, 1);
        one.contact_dwell_fraction = local_fraction_true(contact, valid);
        one.contact_transition_count = local_transition_count(contact, valid);
        one.first_contact_latency_s = local_first_latency(contact, valid, fps);
        one.mutual_approach_dwell_fraction = local_fraction_true(mutualApproach, valid);
        one.withdrawal_dwell_fraction = local_fraction_true(withdrawal, valid);
        one.approach_withdraw_transition_count = ...
            local_state_transition_count(mutualApproach, withdrawal, valid);
        one.asym_positive_dwell_fraction = local_fraction_condition(asym > 0.5, valid & isfinite(asym));
        one.asym_negative_dwell_fraction = local_fraction_condition(asym < -0.5, valid & isfinite(asym));
        one.role_asymmetry_bias_mean = mean(asym(valid & isfinite(asym)), 'omitnan');
        one.centroid_distance_mean_mm = local_mean_valid(dist, valid);
        one.centroid_distance_min_mm = local_min_valid(dist, valid);
        one.centroid_distance_delta_mm = local_late_early_delta(dist, valid);
        one.radial_speed_mean_mm_s = local_mean_valid(radial, valid);
        one.radial_speed_sign_change_count = local_sign_change_count(radial, valid);
        one.approach1_positive_fraction = local_fraction_condition(approach1 > 0, valid & isfinite(approach1));
        one.approach2_positive_fraction = local_fraction_condition(approach2 > 0, valid & isfinite(approach2));
        one.role_approach_imbalance = one.approach1_positive_fraction - one.approach2_positive_fraction;
        one.mutual_facing_mean = local_mean_valid(mutualFacing, valid);
        one.heading_abs_change_deg = local_heading_abs_change(heading, valid);
        one.heading_large_turn_count = local_heading_turn_count(heading, valid, P.turnThresholdDeg);
        one.event_summary_rule = "canonical_dyadic_event_features_frame_mask_only_no_labels";
        one.labels_used_for_event_summary = "none";
        one.arena_used_for_event_summary = false;
        one.condition_used_for_event_summary = false;
        eventSummary = [eventSummary; one]; %#ok<AGROW>
    end
end
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
