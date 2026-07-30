function Measurement = extract_run10_independent_behavioral_measures( ...
        repoRoot, params, Validation, Registry, Sample)
%EXTRACT_RUN10_INDEPENDENT_BEHAVIORAL_MEASURES Disjoint-landmark panel.

assert(Validation.membership_sha256 == params.expected_membership_sha256 && ...
    all(Registry.included_in_primary_validation), ...
    'extract_run10_independent_behavioral_measures:GateNotPassed', ...
    'Freeze and lineage gates must pass before pose files are read.');
assert(all(Sample.experimental_labels_used == "none"), ...
    'extract_run10_independent_behavioral_measures:SamplingFirewallFailure', ...
    'The sample manifest violates the no-label contract.');

repoRoot = string(repoRoot);
Selected = Sample(Sample.sample_selected, :);
Selected = sortrows(Selected, {'preprocess_output_file','anchor_frame', ...
    'graph_node_id'});
featureNames = string(Registry.feature_name(:))';
expectedNames = [ ...
    "mean_forepaw_span_mm","absdiff_forepaw_span_mm", ...
    "mean_hindpaw_span_mm","absdiff_hindpaw_span_mm", ...
    "mean_paw_stance_length_mm","absdiff_paw_stance_length_mm", ...
    "mean_normalized_paw_area","absdiff_normalized_paw_area", ...
    "mean_axial_curvature_deg","absdiff_axial_curvature_deg", ...
    "median_tailtip_pair_distance_mm","q10_cross_paw_distance_mm", ...
    "q10_forepaw_partner_tailtip_distance_mm","paw_area_synchrony", ...
    "axial_curvature_synchrony","stance_change_synchrony"];
assert(isequal(featureNames, expectedNames), ...
    'extract_run10_independent_behavioral_measures:PanelMismatch', ...
    'Extractor equations do not match the frozen feature registry order.');

n = height(Selected);
X = nan(n, numel(featureNames));
poseWindowValidFraction = nan(n, 1);
windowStartFrame = nan(n, 1);
windowStopFrame = nan(n, 1);
windowFrameCount = nan(n, 1);
sourceSha256 = strings(n, 1);
fpsValue = nan(n, 1);
pixelSizeMM = nan(n, 1);

[nodeMap, ~, ~] = default_sleap_node_map();
requiredNodeNames = params.validation_pose_nodes;
for i = 1:numel(requiredNodeNames)
    assert(isfield(nodeMap, requiredNodeNames(i)), ...
        'extract_run10_independent_behavioral_measures:MissingNodeMap', ...
        'Canonical node map lacks %s.', requiredNodeNames(i));
end
requiredNodeIndices = zeros(1, numel(requiredNodeNames));
for i = 1:numel(requiredNodeNames)
    requiredNodeIndices(i) = nodeMap.(requiredNodeNames(i));
end

[fileNames, ~, fileGroup] = unique(Selected.preprocess_output_file, ...
    'stable');
for iFile = 1:numel(fileNames)
    rows = find(fileGroup == iFile);
    pathText = string(resolve_repo_path(repoRoot, fileNames(iFile)));
    assert(isfile(pathText), ...
        'extract_run10_independent_behavioral_measures:MissingPoseFile', ...
        'Missing preprocessed pose input: %s', pathText);
    fileHash = compute_file_sha256(pathText);
    S = load(pathText, 'sessionPreproc');
    assert(isfield(S, 'sessionPreproc') && ...
        isfield(S.sessionPreproc, 'clean') && ...
        isfield(S.sessionPreproc.clean, 'tracks') && ...
        isfield(S.sessionPreproc, 'qc') && ...
        isfield(S.sessionPreproc.qc, 'badframes'), ...
        'extract_run10_independent_behavioral_measures:BadPoseSchema', ...
        'Preprocessed pose file has an unexpected schema: %s', pathText);
    P = S.sessionPreproc;
    tracks = double(P.clean.tracks);
    assert(ndims(tracks) == 4 && size(tracks,3) == 2 && ...
        size(tracks,4) == 2 && ...
        size(tracks,2) >= max(requiredNodeIndices), ...
        'extract_run10_independent_behavioral_measures:BadTrackShape', ...
        'Pose tracks must be frame-by-node-by-xy-by-two-animals.');
    bad = logical(P.qc.badframes);
    if isvector(bad)
        bad = bad(:);
    end
    badPair = any(bad(:, 1:min(2, size(bad,2))), 2);
    fps = double(P.params.data.fps);
    px2mm = double(P.params.data.pixel_size_mm);
    assert(isfinite(fps) && fps > 0 && isfinite(px2mm) && px2mm > 0, ...
        'extract_run10_independent_behavioral_measures:BadUnits', ...
        'Pose fps and pixel size must be finite and positive.');
    F = i_frame_features(tracks, badPair, nodeMap, px2mm);
    halfFrames = round(0.5 * params.validation_window_sec * fps);

    for k = 1:numel(rows)
        row = rows(k);
        anchor = round(Selected.anchor_frame(row));
        st = max(1, anchor - halfFrames);
        en = min(size(tracks,1), anchor + halfFrames);
        idx = st:en;
        windowStartFrame(row) = st;
        windowStopFrame(row) = en;
        windowFrameCount(row) = numel(idx);
        poseWindowValidFraction(row) = mean(F.pose_valid(idx));
        sourceSha256(row) = fileHash;
        fpsValue(row) = fps;
        pixelSizeMM(row) = px2mm;

        X(row,1) = i_summary_mean(F.mean_forepaw_span(idx), ...
            params.minimum_feature_finite_fraction);
        X(row,2) = i_summary_mean(F.absdiff_forepaw_span(idx), ...
            params.minimum_feature_finite_fraction);
        X(row,3) = i_summary_mean(F.mean_hindpaw_span(idx), ...
            params.minimum_feature_finite_fraction);
        X(row,4) = i_summary_mean(F.absdiff_hindpaw_span(idx), ...
            params.minimum_feature_finite_fraction);
        X(row,5) = i_summary_mean(F.mean_stance(idx), ...
            params.minimum_feature_finite_fraction);
        X(row,6) = i_summary_mean(F.absdiff_stance(idx), ...
            params.minimum_feature_finite_fraction);
        X(row,7) = i_summary_mean(F.mean_paw_area(idx), ...
            params.minimum_feature_finite_fraction);
        X(row,8) = i_summary_mean(F.absdiff_paw_area(idx), ...
            params.minimum_feature_finite_fraction);
        X(row,9) = i_summary_mean(F.mean_curvature(idx), ...
            params.minimum_feature_finite_fraction);
        X(row,10) = i_summary_mean(F.absdiff_curvature(idx), ...
            params.minimum_feature_finite_fraction);
        X(row,11) = i_summary_quantile(F.tailtip_distance(idx), 50, ...
            params.minimum_feature_finite_fraction);
        X(row,12) = i_summary_quantile(F.cross_paw_distance(idx), 10, ...
            params.minimum_feature_finite_fraction);
        X(row,13) = i_summary_quantile( ...
            F.forepaw_partner_tailtip_distance(idx), 10, ...
            params.minimum_feature_finite_fraction);
        X(row,14) = i_corr(F.paw_area_1(idx), F.paw_area_2(idx), ...
            params.minimum_correlation_frames, ...
            params.minimum_feature_finite_fraction);
        X(row,15) = i_corr(F.curvature_1(idx), F.curvature_2(idx), ...
            params.minimum_correlation_frames, ...
            params.minimum_feature_finite_fraction);
        X(row,16) = i_corr(diff(F.stance_1(idx)), ...
            diff(F.stance_2(idx)), ...
            max(3, params.minimum_correlation_frames - 1), ...
            params.minimum_feature_finite_fraction);
    end
end

Measurement = Selected(:, {'graph_node_id','motif_candidate_id', ...
    'candidate_local_index','eligible_for_behavioral_interpretation', ...
    'embedding_row_id','scale_index','chunk_sec','session_index', ...
    'session_id','raw_index','anchor_frame','anchor_time_s', ...
    'preprocess_output_file'});
Measurement.window_start_frame = windowStartFrame;
Measurement.window_stop_frame = windowStopFrame;
Measurement.window_frame_count = windowFrameCount;
Measurement.pose_window_valid_fraction = poseWindowValidFraction;
Measurement.source_pose_sha256 = sourceSha256;
Measurement.fps = fpsValue;
Measurement.pixel_size_mm = pixelSizeMM;
FeatureTable = array2table(X, 'VariableNames', cellstr(featureNames));
Measurement = [Measurement FeatureTable];
Measurement.finite_feature_count = sum(isfinite(X), 2);
Measurement.eligible_for_automated_analysis = ...
    Measurement.finite_feature_count >= ...
    params.minimum_finite_feature_count;
Measurement.measurement_status = repmat("sufficient_independent_features", ...
    n, 1);
Measurement.measurement_status( ...
    ~Measurement.eligible_for_automated_analysis) = ...
    "insufficient_independent_features";
Measurement.validation_feature_panel_version = repmat( ...
    params.validation_feature_panel_version, n, 1);
Measurement.validation_window_sec = repmat( ...
    params.validation_window_sec, n, 1);
Measurement.candidate_freeze_id = repmat( ...
    params.expected_candidate_freeze_id, n, 1);
Measurement.frozen_membership_sha256 = repmat( ...
    params.expected_membership_sha256, n, 1);
Measurement.run05_to_run09_values_used = repmat("none", n, 1);
Measurement.experimental_labels_used = repmat("none", n, 1);
Measurement.arena_used = false(n, 1);
Measurement.annotation_used = false(n, 1);
Measurement.graph_stability_used = false(n, 1);

assert(all(Measurement.run05_to_run09_values_used == "none") && ...
    all(Measurement.experimental_labels_used == "none") && ...
    ~any(Measurement.arena_used) && ~any(Measurement.annotation_used) && ...
    ~any(Measurement.graph_stability_used), ...
    'extract_run10_independent_behavioral_measures:FirewallFailure', ...
    'Non-independent or forbidden metadata entered measurement extraction.');
Measurement = sortrows(Measurement, 'graph_node_id');
end

function F = i_frame_features(tracks, badPair, nodeMap, px2mm)
node = @(name, animal) i_xy(tracks, nodeMap.(name), animal, badPair);
lf1 = node('lfPaw',1); rf1 = node('rfPaw',1);
lh1 = node('lhPaw',1); rh1 = node('rhPaw',1);
lf2 = node('lfPaw',2); rf2 = node('rfPaw',2);
lh2 = node('lhPaw',2); rh2 = node('rhPaw',2);
neck1 = node('neck',1); neck2 = node('neck',2);
tm1 = node('tailMid',1); tm2 = node('tailMid',2);
tt1 = node('tailTip',1); tt2 = node('tailTip',2);

fore1 = i_dist(lf1, rf1) * px2mm;
fore2 = i_dist(lf2, rf2) * px2mm;
hind1 = i_dist(lh1, rh1) * px2mm;
hind2 = i_dist(lh2, rh2) * px2mm;
frontMid1 = (lf1 + rf1) ./ 2;
frontMid2 = (lf2 + rf2) ./ 2;
hindMid1 = (lh1 + rh1) ./ 2;
hindMid2 = (lh2 + rh2) ./ 2;
stance1 = i_dist(frontMid1, hindMid1) * px2mm;
stance2 = i_dist(frontMid2, hindMid2) * px2mm;
area1 = i_polygon_area(lf1, rf1, rh1, lh1) * px2mm^2;
area2 = i_polygon_area(lf2, rf2, rh2, lh2) * px2mm^2;
areaNorm1 = area1 ./ max(stance1.^2, eps);
areaNorm2 = area2 ./ max(stance2.^2, eps);
curv1 = i_angle_between(tm1 - neck1, tt1 - tm1);
curv2 = i_angle_between(tm2 - neck2, tt2 - tm2);

paws1 = {lf1,rf1,lh1,rh1};
paws2 = {lf2,rf2,lh2,rh2};
crossPaw = inf(size(fore1));
for i = 1:4
    for j = 1:4
        crossPaw = min(crossPaw, i_dist(paws1{i}, paws2{j}) * px2mm);
    end
end
crossPaw(~isfinite(crossPaw)) = NaN;
foreTail = min([i_dist(lf1,tt2), i_dist(rf1,tt2), ...
    i_dist(lf2,tt1), i_dist(rf2,tt1)], [], 2) * px2mm;

F = struct();
F.mean_forepaw_span = (fore1 + fore2) ./ 2;
F.absdiff_forepaw_span = abs(fore1 - fore2);
F.mean_hindpaw_span = (hind1 + hind2) ./ 2;
F.absdiff_hindpaw_span = abs(hind1 - hind2);
F.mean_stance = (stance1 + stance2) ./ 2;
F.absdiff_stance = abs(stance1 - stance2);
F.mean_paw_area = (areaNorm1 + areaNorm2) ./ 2;
F.absdiff_paw_area = abs(areaNorm1 - areaNorm2);
F.mean_curvature = (curv1 + curv2) ./ 2;
F.absdiff_curvature = abs(curv1 - curv2);
F.tailtip_distance = i_dist(tt1, tt2) * px2mm;
F.cross_paw_distance = crossPaw;
F.forepaw_partner_tailtip_distance = foreTail;
F.paw_area_1 = areaNorm1;
F.paw_area_2 = areaNorm2;
F.curvature_1 = curv1;
F.curvature_2 = curv2;
F.stance_1 = stance1;
F.stance_2 = stance2;
allXY = [lf1 rf1 lh1 rh1 lf2 rf2 lh2 rh2 neck1 neck2 tm1 tm2 tt1 tt2];
F.pose_valid = ~badPair & all(isfinite(allXY), 2);
end

function xy = i_xy(tracks, nodeIndex, animal, badPair)
xy = reshape(tracks(:, nodeIndex, :, animal), size(tracks,1), 2);
xy(badPair, :) = NaN;
end

function d = i_dist(a, b)
d = sqrt(sum((a - b).^2, 2));
end

function area = i_polygon_area(a, b, c, d)
x = [a(:,1) b(:,1) c(:,1) d(:,1)];
y = [a(:,2) b(:,2) c(:,2) d(:,2)];
area = 0.5 .* abs(sum(x .* y(:,[2 3 4 1]) - ...
    y .* x(:,[2 3 4 1]), 2));
bad = any(~isfinite(x) | ~isfinite(y), 2);
area(bad) = NaN;
end

function angle = i_angle_between(a, b)
den = sqrt(sum(a.^2,2)) .* sqrt(sum(b.^2,2));
c = sum(a .* b, 2) ./ den;
c = max(-1, min(1, c));
angle = acosd(c);
angle(~isfinite(den) | den <= 0) = NaN;
end

function value = i_summary_mean(x, minFraction)
ok = isfinite(x);
if mean(ok) < minFraction
    value = NaN;
else
    value = mean(x(ok));
end
end

function value = i_summary_quantile(x, percentile, minFraction)
ok = isfinite(x);
if mean(ok) < minFraction
    value = NaN;
else
    value = prctile(x(ok), percentile);
end
end

function value = i_corr(x, y, minFrames, minFraction)
ok = isfinite(x) & isfinite(y);
if nnz(ok) < minFrames || mean(ok) < minFraction || ...
        std(x(ok)) <= eps || std(y(ok)) <= eps
    value = NaN;
else
    C = corrcoef(x(ok), y(ok));
    value = C(1,2);
end
end
