% demo_run06_multiscale_chunking_principles
% Condition-blind toy walkthrough for run_06 multiscale chunking.
%
% This script does not use project data. It creates synthetic dyadic
% features with known toy interaction phases so the figures can explain the
% run_06 principles. The toy phase labels are for visualization only and are
% never used for transforms, anchor selection, scale scoring, or PCA.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
outRoot = fullfile(repoRoot, 'derived', 'demo_run06_multiscale_chunking_principles');
figRoot = fullfile(outRoot, 'figures');
if ~exist(figRoot, 'dir')
    mkdir(figRoot);
end

rng(106);
fps = 80;
durationSec = 240;
nFrames = durationSec * fps;
t = (0:nFrames-1)' ./ fps;
primaryMicroSec = [0.2215 0.3011 0.4092 0.5561 0.7557]';
primaryMotifSec = [0.9272 1.2601 1.7125 2.3274]';
primaryContextSec = [2.5780 4.2986 5.8420 7.9395]';
surveyOnlyContextSec = [10.7900 17.9917 24.4513 30.0000]';
scaleSec = [primaryMicroSec; primaryMotifSec; primaryContextSec; surveyOnlyContextSec];
strideSec = 0.25;
strideFrames = round(strideSec * fps);
minValidFrac = 0.85;
minFeatureFiniteFrac = 0.95;
maxPcaChunks = 180;
summaryMicroTemporalBins = 6;
summaryMicroDctCoeffs = 4;
summaryMotifTemporalBins = 12;
summaryMotifDctCoeffs = 8;
summaryContextTemporalBins = 6;
summaryContextDctCoeffs = 4;
scorePcsRetained = 12;
embeddingDimMicroMaxPcs = 24;
embeddingDimMotifMaxPcs = 48;
embeddingDimContextMaxPcs = 64;
dimensionGuardMinEffectiveDim = 2;
dimensionGuardMaxPc1Pct = 80;
fixedFivePcVarianceTargetPct = 80;
exampleScaleSec = [0.2215 1.2601 5.8420 24.4513]';

featureNames = [
    "centroid_dist"
    "body2body_dist"
    "head2head_dist"
    "tailbase2tailbase_dist"
    "nose1_to_body2_dist"
    "nose2_to_body1_dist"
    "nose1_to_tail2_dist"
    "nose2_to_tail1_dist"
    "partner_long_1"
    "partner_lat_1"
    "partner_long_2"
    "partner_lat_2"
    "facing_1_to_2"
    "facing_2_to_1"
    "mutual_facing"
    "heading_diff_deg_cos"
    "heading_diff_deg_sin"
    "cos_head_alignment"
    "radial_speed_12"
    "tangential_speed_12"
    "approach_speed_1"
    "approach_speed_2"
    "speed_alignment"
    "accel_alignment"
    "nose_bearing_1_deg_cos"
    "nose_bearing_1_deg_sin"
    "nose_bearing_2_deg_cos"
    "nose_bearing_2_deg_sin"
    "body_bearing_1_deg_cos"
    "body_bearing_1_deg_sin"
    "body_bearing_2_deg_cos"
    "body_bearing_2_deg_sin"
    "in_contact"
    "head_close"
    "body_close"
    "close_pair"
    "mutual_approach"
    "withdrawal"
    "asym_investigate"
    ];
nFeature = numel(featureNames);

stateNames = ["apart","approach","contact","investigate","parallel","withdraw"];
statePattern = [1 2 3 4 5 6 1 5 2 3 6 1 4 2 3 5 6 1];
durationPatternSec = [12 4 6 7 10 4 14 12 3 7 5 11 8 4 6 14 5 13];
stateId = zeros(nFrames, 1);
segments = table();
frame0 = 1;
seg = 1;
while frame0 <= nFrames
    state = statePattern(1 + mod(seg-1, numel(statePattern)));
    durFrames = round(durationPatternSec(1 + mod(seg-1, numel(durationPatternSec))) * fps);
    frame1 = min(nFrames, frame0 + durFrames - 1);
    stateId(frame0:frame1) = state;
    one = table();
    one.segment_index = seg;
    one.start_frame = frame0;
    one.stop_frame = frame1;
    one.start_time_s = (frame0 - 1) / fps;
    one.stop_time_s = (frame1 - 1) / fps;
    one.toy_state = stateNames(state);
    one.toy_state_used_for_run06_logic = false;
    segments = [segments; one]; %#ok<AGROW>
    frame0 = frame1 + 1;
    seg = seg + 1;
end

dist = nan(nFrames, 1);
mutualFacing = nan(nFrames, 1);
radialSpeed = nan(nFrames, 1);
tangentialSpeed = nan(nFrames, 1);
headingDeg = nan(nFrames, 1);
asym = nan(nFrames, 1);
inContact = false(nFrames, 1);

for i = 1:height(segments)
    idx = segments.start_frame(i):segments.stop_frame(i);
    u = linspace(0, 1, numel(idx))';
    switch string(segments.toy_state(i))
        case "apart"
            d = 260 + 12*sin(2*pi*u);
            face = -0.35 + 0.1*cos(2*pi*u);
            rs = 2*sin(2*pi*u);
            ts = 18*cos(2*pi*u);
            hd = 150 + 20*sin(2*pi*u);
            a = 0;
            contact = false(numel(idx), 1);
        case "approach"
            d = 250 - 165*u + 8*sin(6*pi*u);
            face = -0.1 + 0.9*u;
            rs = -85 + 12*sin(4*pi*u);
            ts = 25*sin(2*pi*u);
            hd = 120 - 80*u;
            a = 0.2;
            contact = false(numel(idx), 1);
        case "contact"
            d = 58 + 9*sin(4*pi*u);
            face = 0.75 + 0.15*sin(2*pi*u);
            rs = 5*sin(5*pi*u);
            ts = 15*sin(7*pi*u);
            hd = 30 + 25*sin(4*pi*u);
            a = 0.1*sin(2*pi*u);
            contact = true(numel(idx), 1);
        case "investigate"
            d = 85 + 20*sin(3*pi*u);
            face = 0.65 + 0.25*sin(5*pi*u);
            rs = 20*sin(6*pi*u);
            ts = 40*cos(3*pi*u);
            hd = 55 + 85*sin(3*pi*u);
            a = 0.85*sin(pi*u);
            contact = d < 70;
        case "parallel"
            d = 130 + 15*sin(2*pi*u);
            face = 0.15 + 0.1*sin(2*pi*u);
            rs = 0.5*sin(2*pi*u);
            ts = 55 + 10*sin(4*pi*u);
            hd = 10 + 15*sin(2*pi*u);
            a = -0.15;
            contact = false(numel(idx), 1);
        otherwise
            d = 70 + 180*u + 8*sin(6*pi*u);
            face = 0.7 - 0.8*u;
            rs = 80 + 15*sin(4*pi*u);
            ts = 20*sin(2*pi*u);
            hd = 35 + 120*u;
            a = -0.25;
            contact = false(numel(idx), 1);
    end
    dist(idx) = d;
    mutualFacing(idx) = face;
    radialSpeed(idx) = rs;
    tangentialSpeed(idx) = ts;
    headingDeg(idx) = hd;
    asym(idx) = a;
    inContact(idx) = contact;
end

slowContext = 1 + 0.12*sin(2*pi*t/80);
dist = dist .* slowContext + 4*randn(nFrames, 1);
mutualFacing = max(-1, min(1, mutualFacing + 0.08*randn(nFrames, 1)));
radialSpeed = radialSpeed + 8*randn(nFrames, 1);
tangentialSpeed = tangentialSpeed + 6*randn(nFrames, 1);
headingDeg = mod(headingDeg + 12*randn(nFrames, 1) + 180, 360) - 180;
asym = max(-1, min(1, asym + 0.12*randn(nFrames, 1)));

X = zeros(nFrames, nFeature);
X(:,1) = dist;
X(:,2) = dist - 12 + 5*randn(nFrames,1);
X(:,3) = dist - 18 + 8*randn(nFrames,1);
X(:,4) = dist + 15 + 5*randn(nFrames,1);
X(:,5) = dist - 5 + 9*randn(nFrames,1);
X(:,6) = dist - 3 + 9*randn(nFrames,1);
X(:,7) = dist + 22 + 10*randn(nFrames,1);
X(:,8) = dist + 20 + 10*randn(nFrames,1);
X(:,9) = 0.75*dist.*cosd(headingDeg) + 25*randn(nFrames,1);
X(:,10) = 0.75*dist.*sind(headingDeg) + 25*randn(nFrames,1);
X(:,11) = 0.70*dist.*cosd(headingDeg + 15) + 25*randn(nFrames,1);
X(:,12) = 0.70*dist.*sind(headingDeg - 10) + 25*randn(nFrames,1);
X(:,13) = max(-1, min(1, mutualFacing + 0.12*randn(nFrames,1)));
X(:,14) = max(-1, min(1, mutualFacing + 0.12*randn(nFrames,1)));
X(:,15) = mutualFacing;
X(:,16) = cosd(headingDeg);
X(:,17) = sind(headingDeg);
X(:,18) = cosd(headingDeg);
X(:,19) = radialSpeed;
X(:,20) = tangentialSpeed;
X(:,21) = -radialSpeed;
X(:,22) = -radialSpeed + 3*randn(nFrames,1);
X(:,23) = max(-1, min(1, 0.6*cosd(headingDeg) + 0.2*randn(nFrames,1)));
X(:,24) = max(-1, min(1, 0.4*cosd(headingDeg) + 0.25*randn(nFrames,1)));
X(:,25) = cosd(headingDeg + 25);
X(:,26) = sind(headingDeg + 25);
X(:,27) = cosd(headingDeg - 25);
X(:,28) = sind(headingDeg - 25);
X(:,29) = cosd(headingDeg + 45);
X(:,30) = sind(headingDeg + 45);
X(:,31) = cosd(headingDeg - 45);
X(:,32) = sind(headingDeg - 45);
X(:,33) = double(inContact);
X(:,34) = double(dist < 80);
X(:,35) = double(dist < 95);
X(:,36) = double(dist < 110);
X(:,37) = double(radialSpeed < -35 & mutualFacing > 0.2);
X(:,38) = double(radialSpeed > 35);
X(:,39) = asym;

validMask = true(nFrames, 1);
gapSec = [35 37; 110 114; 176 183];
for g = 1:size(gapSec,1)
    validMask((gapSec(g,1)*fps):(gapSec(g,2)*fps)) = false;
end
X(~validMask, :) = NaN;

isBoolean = ismember(featureNames, ["in_contact","head_close","body_close","close_pair","mutual_approach","withdrawal"]);
isDistance = contains(featureNames, "dist");
Xtrans = X;
transformRows = table();
for j = 1:nFeature
    x = X(:,j);
    fitMask = validMask & isfinite(x);
    transformName = "robust_zscore";
    if isBoolean(j)
        xt = x;
        center = 0;
        scale = 1;
        transformName = "binary_identity";
    else
        if isDistance(j)
            x = log1p(max(x, 0));
            transformName = "log1p_then_robust_zscore";
        end
        fitVals = x(fitMask);
        center = median(fitVals, 'omitnan');
        scale = iqr(fitVals);
        if ~(isfinite(scale) && scale > 0)
            scale = std(fitVals, 0, 'omitnan');
        end
        if ~(isfinite(scale) && scale > 0)
            scale = 1;
        end
        xt = (x - center) ./ scale;
    end
    Xtrans(:,j) = xt;

    row = table();
    row.channel_index = j;
    row.feature_name = featureNames(j);
    row.transform_name = transformName;
    row.robust_center = center;
    row.robust_scale = scale;
    row.n_valid_fit_rows = nnz(fitMask);
    row.frame_mask_used_for_fit = true;
    row.toy_state_used_for_fit = false;
    row.provenance_label_used_for_fit = false;
    transformRows = [transformRows; row]; %#ok<AGROW>
end

anchorRows = table();
pcaRows = table();
anchorsByScale = cell(numel(scaleSec), 1);
chunkExamples = table();
featuresToPlot = ["centroid_dist","mutual_facing","radial_speed_12","in_contact"];
[~, featurePlotIdx] = ismember(featuresToPlot, featureNames);

pcaWarnState = warning('off', 'stats:pca:ColRankDefX');
pcaWarnCleanup = onCleanup(@() warning(pcaWarnState)); %#ok<NASGU>

for s = 1:numel(scaleSec)
    chunkFrames = max(1, round(scaleSec(s) * fps));
    leftFrames = chunkFrames - 1;
    candidateFrames = (1 + leftFrames):strideFrames:nFrames;
    keep = false(numel(candidateFrames), 1);
    validFrac = nan(numel(candidateFrames), 1);
    featureFrac = nan(numel(candidateFrames), 1);

    for i = 1:numel(candidateFrames)
        f = candidateFrames(i);
        idx = (f-leftFrames):f;
        validFrac(i) = mean(validMask(idx));
        featureFrac(i) = mean(isfinite(Xtrans(idx,:)), 'all');
        keep(i) = validMask(f) && validFrac(i) >= minValidFrac && featureFrac(i) >= minFeatureFiniteFrac;
    end

    anchors = candidateFrames(keep)';
    anchorsByScale{s} = anchors(:);
    band = "context";
    if scaleSec(s) < 0.8
        band = "micro";
    elseif scaleSec(s) < 2.5
        band = "motif";
    end
    if band == "micro"
        summaryTemporalBins = summaryMicroTemporalBins;
        summaryDctCoeffs = summaryMicroDctCoeffs;
        recommendedPcCap = embeddingDimMicroMaxPcs;
    elseif band == "motif"
        summaryTemporalBins = summaryMotifTemporalBins;
        summaryDctCoeffs = summaryMotifDctCoeffs;
        recommendedPcCap = embeddingDimMotifMaxPcs;
    else
        summaryTemporalBins = summaryContextTemporalBins;
        summaryDctCoeffs = summaryContextDctCoeffs;
        recommendedPcCap = embeddingDimContextMaxPcs;
    end
    summaryProfile = sprintf('%s_adaptive_bins%d_dct%d', band, summaryTemporalBins, summaryDctCoeffs);
    productionPrimaryStatus = ismember(scaleSec(s), [primaryMicroSec; primaryMotifSec; primaryContextSec]);
    productionScoreSelectedStatus = productionPrimaryStatus || ismember(scaleSec(s), surveyOnlyContextSec(1:3));
    primaryExclusionReason = "primary";
    if productionScoreSelectedStatus && ~productionPrimaryStatus
        primaryExclusionReason = "score_selected_but_rejected_by_dimension_or_stability_guard";
    elseif ~productionScoreSelectedStatus
        primaryExclusionReason = "survey_only_not_score_selected";
    end

    one = table();
    one.scale_index = s;
    one.chunk_sec = scaleSec(s);
    one.chunk_frames = chunkFrames;
    one.temporal_band = band;
    one.summary_profile = string(summaryProfile);
    one.summary_temporal_bins = summaryTemporalBins;
    one.summary_dct_coeffs = summaryDctCoeffs;
    one.production_score_selected_scale_analog = productionScoreSelectedStatus;
    one.production_primary_scale_analog = productionPrimaryStatus;
    one.primary_exclusion_reason_analog = primaryExclusionReason;
    one.n_stride_candidates = numel(candidateFrames);
    one.n_scale_specific_candidate_anchors = numel(anchors);
    one.n_candidates_removed_by_frame_mask = nnz(~keep);
    one.mean_valid_fraction = mean(validFrac(keep), 'omitnan');
    one.anchor_rule = "frameMask plus finite transformed features plus fixed stride";
    one.toy_state_used_for_anchor_selection = false;
    one.provenance_label_used_for_anchor_selection = false;
    anchorRows = [anchorRows; one]; %#ok<AGROW>

    pcaAnchors = anchors;
    if numel(pcaAnchors) > maxPcaChunks
        pick = unique(round(linspace(1, numel(pcaAnchors), maxPcaChunks)));
        pcaAnchors = pcaAnchors(pick);
    end

    n = numel(pcaAnchors);
    Xchunk = nan(n, chunkFrames, nFeature);
    for i = 1:n
        idx = (pcaAnchors(i)-leftFrames):pcaAnchors(i);
        Xi = Xtrans(idx, :);
        med = median(Xi, 1, 'omitnan');
        for j = 1:nFeature
            miss = ~isfinite(Xi(:,j));
            Xi(miss,j) = med(j);
        end
        Xchunk(i,:,:) = Xi;
    end

    binEdges = unique(round(linspace(1, chunkFrames + 1, summaryTemporalBins + 1)));
    relT = linspace(-1, 1, chunkFrames);
    relT = relT - mean(relT);
    slopeDen = sum(relT .^ 2);
    dctBasis = zeros(summaryDctCoeffs, chunkFrames);
    nTime = 0:(chunkFrames - 1);
    for k = 1:summaryDctCoeffs
        dctBasis(k,:) = sqrt(2 / chunkFrames) .* cos(pi .* (nTime + 0.5) .* k ./ chunkFrames);
    end

    Xsummary = [];
    for j = 1:nFeature
        Xc = squeeze(Xchunk(:,:,j));
        rowMean = mean(Xc, 2, 'omitnan');
        rowStd = std(Xc, 0, 2, 'omitnan');
        rowMedian = median(Xc, 2, 'omitnan');
        rowMean(~isfinite(rowMean)) = 0;
        rowStd(~isfinite(rowStd)) = 0;
        rowMedian(~isfinite(rowMedian)) = rowMean(~isfinite(rowMedian));
        Xfill = Xc;
        for ii = 1:n
            miss = ~isfinite(Xfill(ii,:));
            if any(miss)
                Xfill(ii,miss) = rowMedian(ii);
            end
        end
        q10 = prctile(Xfill, 10, 2);
        q90 = prctile(Xfill, 90, 2);
        slope = (Xfill * relT(:)) ./ max(slopeDen, eps);
        earlyIdx = 1:max(1, floor(chunkFrames / 3));
        lateIdx = max(1, chunkFrames - floor(chunkFrames / 3) + 1):chunkFrames;
        deltaLateEarly = mean(Xfill(:, lateIdx), 2) - mean(Xfill(:, earlyIdx), 2);
        Xsummary = [Xsummary, rowMean, rowStd, rowMedian, q10, q90, slope, deltaLateEarly]; %#ok<AGROW>

        for b = 1:(numel(binEdges) - 1)
            idxBin = binEdges(b):(binEdges(b + 1) - 1);
            idxBin = idxBin(idxBin >= 1 & idxBin <= chunkFrames);
            Xsummary = [Xsummary, mean(Xc(:, idxBin), 2, 'omitnan')]; %#ok<AGROW>
        end

        Xcentered = Xfill - mean(Xfill, 2);
        Xsummary = [Xsummary, Xcentered * dctBasis']; %#ok<AGROW>

        if isBoolean(j)
            xb = Xfill > 0.5;
            transitionRate = sum(abs(diff(double(xb), 1, 2)), 2) ./ max(chunkFrames - 1, 1);
            Xsummary = [Xsummary, transitionRate]; %#ok<AGROW>
        end
    end

    Xz = zscore(Xsummary, 0, 1);
    Xz(~isfinite(Xz)) = 0;
    [~, ~, ~, ~, explained] = pca(Xz, 'Centered', false);
    cumExplained = cumsum(explained);
    n90 = find(cumExplained >= 90, 1, 'first');
    n95 = find(cumExplained >= 95, 1, 'first');
    if isempty(n90), n90 = NaN; end
    if isempty(n95), n95 = NaN; end
    eigFrac = explained ./ max(sum(explained), eps);
    pc1 = explained(1);

    prow = table();
    prow.scale_index = s;
    prow.chunk_sec = scaleSec(s);
    prow.temporal_band = band;
    prow.n_observed_chunks = numel(anchors);
    prow.n_pca_chunks = n;
    prow.raw_flattened_dimensions = chunkFrames * nFeature;
    prow.summary_dimensions = size(Xsummary, 2);
    prow.compression_ratio_raw_to_summary = prow.raw_flattened_dimensions ./ max(prow.summary_dimensions, 1);
    prow.summary_profile = string(summaryProfile);
    prow.summary_temporal_bins = summaryTemporalBins;
    prow.summary_dct_coeffs = summaryDctCoeffs;
    prow.cum5_explained = sum(explained(1:min(5,numel(explained))));
    prow.cum12_explained = sum(explained(1:min(scorePcsRetained,numel(explained))));
    prow.pc1_explained = pc1;
    prow.n_pcs_90pct = n90;
    prow.n_pcs_95pct = n95;
    prow.effective_dim = 1 ./ max(sum(eigFrac.^2), eps);
    prow.score_pcs_retained = scorePcsRetained;
    prow.recommended_max_pcs_by_temporal_role = recommendedPcCap;
    prow.recommended_pcs_for_next_embedding = min(max(n90, scorePcsRetained), recommendedPcCap);
    prow.recommendation_reaches_90pct = n90 <= recommendedPcCap;
    prow.dimension_guard_min_effective_dim = dimensionGuardMinEffectiveDim;
    prow.dimension_guard_max_pc1_explained_pct = dimensionGuardMaxPc1Pct;
    prow.dominant_pc_warning = pc1 > dimensionGuardMaxPc1Pct | prow.effective_dim < dimensionGuardMinEffectiveDim;
    prow.passes_dimension_guard = ~prow.dominant_pc_warning;
    prow.representation_mode = "multiresolution";
    prow.production_score_selected_scale_analog = productionScoreSelectedStatus;
    prow.production_primary_scale_analog = productionPrimaryStatus;
    prow.primary_exclusion_reason_analog = primaryExclusionReason;
    prow.labels_used_for_pca = "none";
    prow.fixed_five_pc_variance_target_pct = fixedFivePcVarianceTargetPct;
    prow.fixed_five_pc_is_sufficient = prow.cum5_explained >= fixedFivePcVarianceTargetPct;
    pcaRows = [pcaRows; prow]; %#ok<AGROW>
end

commonAnchors = anchorsByScale{end};
tolFrames = ceil(0.51 * strideFrames);
keepCommon = true(numel(commonAnchors), 1);
for s = 1:(numel(anchorsByScale)-1)
    A = anchorsByScale{s};
    for i = 1:numel(commonAnchors)
        keepCommon(i) = keepCommon(i) && any(abs(A - commonAnchors(i)) <= tolFrames);
    end
end
commonAnchors = commonAnchors(keepCommon);
anchorRows.n_common_all_scale_candidate_anchors = repmat(numel(commonAnchors), height(anchorRows), 1);
anchorRows.common_anchor_retention_fraction = anchorRows.n_common_all_scale_candidate_anchors ./ ...
    max(anchorRows.n_scale_specific_candidate_anchors, 1);

if isempty(commonAnchors)
    candidateExampleAnchors = anchorsByScale{end};
else
    candidateExampleAnchors = commonAnchors;
end
contactScore = zeros(numel(candidateExampleAnchors), 1);
longFrames = round(30 * fps);
for i = 1:numel(candidateExampleAnchors)
    f = candidateExampleAnchors(i);
    idx = max(1, f - longFrames + 1):f;
    contactScore(i) = mean(inContact(idx)) + 0.5 * double(inContact(f));
end
[~, bestExample] = max(contactScore);
exampleAnchor = candidateExampleAnchors(bestExample);
for s = 1:numel(exampleScaleSec)
    [~, scaleIdx] = min(abs(scaleSec - exampleScaleSec(s)));
    chunkFrames = max(1, round(scaleSec(scaleIdx) * fps));
    leftFrames = chunkFrames - 1;
    idx = (exampleAnchor-leftFrames):exampleAnchor;
    idx = idx(idx >= 1 & idx <= nFrames);
    relTime = ((idx(:) - exampleAnchor) ./ fps);
    for f = 1:numel(featurePlotIdx)
        row = table();
        row.scale_index = repmat(scaleIdx, numel(idx), 1);
        row.chunk_sec = repmat(scaleSec(scaleIdx), numel(idx), 1);
        row.feature_name = repmat(featuresToPlot(f), numel(idx), 1);
        row.relative_time_s = relTime;
        row.value = X(idx, featurePlotIdx(f));
        row.transformed_value = Xtrans(idx, featurePlotIdx(f));
        row.valid_mask = validMask(idx);
        toyStateForRows = string(stateNames(stateId(idx(:))));
        row.toy_state = toyStateForRows(:);
        chunkExamples = [chunkExamples; row]; %#ok<AGROW>
    end
end

selectedDemoScales = table();
selectedDemoScales.scale_sec = exampleScaleSec;
selectedDemoScales.role = ["primary_micro"; "primary_motif"; "primary_context"; "score_selected_long_context_rejected"];
selectedDemoScales.what_it_can_capture = [
    "instantaneous posture and contact onset"
    "approach/contact/withdraw fragments"
    "short social episodes and local persistence"
    "slow context across many transitions, useful to audit but not automatically primary"
    ];
selectedDemoScales.condition_blind_selection_note = repmat( ...
    "selected here for explanation only; production run_06 uses condition-blind score, stability, and dimensionality audits", ...
    height(selectedDemoScales), 1);

demoDecisionRows = pcaRows;
demoDecisionRows.n_scale_specific_candidate_anchors = anchorRows.n_scale_specific_candidate_anchors;
demoDecisionRows.n_common_all_scale_candidate_anchors = anchorRows.n_common_all_scale_candidate_anchors;
demoDecisionRows.common_anchor_retention_fraction = anchorRows.common_anchor_retention_fraction;
demoDecisionRows.n_candidates_removed_by_frame_mask = anchorRows.n_candidates_removed_by_frame_mask;
demoDecisionRows.demonstration_selected_scale = ismember(demoDecisionRows.chunk_sec, exampleScaleSec);
demoDecisionRows.selection_rule_for_demo_only = repmat( ...
    "one illustrative primary micro, primary motif, primary context, and rejected long-context scale; production run_06 uses score, stability, and dimensionality audits", ...
    height(demoDecisionRows), 1);
demoDecisionRows.labels_used_for_demo_decision = repmat("none", height(demoDecisionRows), 1);
demoDecisionRows.production_like_selection_story = strings(height(demoDecisionRows), 1);
demoDecisionRows.demo_decision_rationale = strings(height(demoDecisionRows), 1);
for i = 1:height(demoDecisionRows)
    if demoDecisionRows.production_primary_scale_analog(i)
        demoDecisionRows.production_like_selection_story(i) = ...
            "primary scale analog: score-selected and promoted after stability/dimensionality guards";
    elseif demoDecisionRows.production_score_selected_scale_analog(i)
        demoDecisionRows.production_like_selection_story(i) = ...
            "score-selected analog: retained for audit but not promoted to primary";
    else
        demoDecisionRows.production_like_selection_story(i) = ...
            "survey-only analog: evaluated but not score-selected";
    end

    if demoDecisionRows.demonstration_selected_scale(i) && demoDecisionRows.production_primary_scale_analog(i)
        if demoDecisionRows.passes_dimension_guard(i)
            demoDecisionRows.demo_decision_rationale(i) = ...
                "illustrative primary scale with non-dominant compact PCA representation";
        else
            demoDecisionRows.demo_decision_rationale(i) = ...
                "illustrative primary-like scale would need a guard audit before promotion";
        end
    elseif demoDecisionRows.demonstration_selected_scale(i) && demoDecisionRows.production_score_selected_scale_analog(i)
        demoDecisionRows.demo_decision_rationale(i) = ...
            "illustrative long-context scale shows why score-selected candidates still need primary-promotion guards";
    elseif ~demoDecisionRows.fixed_five_pc_is_sufficient(i)
        demoDecisionRows.demo_decision_rationale(i) = ...
            "survey scale shows fixed five PCs would underrepresent compact-summary variance";
    else
        demoDecisionRows.demo_decision_rationale(i) = ...
            "survey scale retained in audit table but not highlighted in demo figures";
    end
end

demoEventRows = table();
for s = 1:numel(scaleSec)
    chunkFrames = max(1, round(scaleSec(s) * fps));
    leftFrames = chunkFrames - 1;
    anchors = anchorsByScale{s};
    if isempty(anchors)
        continue
    end
    if numel(anchors) > maxPcaChunks
        pick = unique(round(linspace(1, numel(anchors), maxPcaChunks)));
        anchors = anchors(pick);
    end
    for i = 1:numel(anchors)
        f = anchors(i);
        idx = (f-leftFrames):f;
        idx = idx(idx >= 1 & idx <= nFrames);
        valid = validMask(idx);
        contact = X(idx, 33) > 0.5;
        approach = X(idx, 37) > 0.5;
        withdraw = X(idx, 38) > 0.5;
        asymChunk = X(idx, 39);
        distChunk = X(idx, 1);
        radialChunk = X(idx, 19);
        facingChunk = X(idx, 15);
        headingChunk = headingDeg(idx);

        validContact = contact(valid);
        validApproach = approach(valid);
        validWithdraw = withdraw(valid);
        validAsym = asymChunk(valid & isfinite(asymChunk));
        validDist = distChunk(valid & isfinite(distChunk));
        validRadial = radialChunk(valid & isfinite(radialChunk));
        validFacing = facingChunk(valid & isfinite(facingChunk));
        validHeading = headingChunk(valid & isfinite(headingChunk));

        contactTransitions = NaN;
        if numel(validContact) > 1
            contactTransitions = sum(abs(diff(double(validContact))) > 0);
        end
        stateCode = zeros(nnz(valid), 1);
        stateCode(validApproach) = 1;
        stateCode(validWithdraw) = -1;
        approachWithdrawTransitions = NaN;
        if numel(stateCode) > 1
            approachWithdrawTransitions = sum(abs(diff(stateCode)) > 0);
        end
        firstContactLatency = NaN;
        firstContact = find(contact & valid, 1, 'first');
        if ~isempty(firstContact)
            firstContactLatency = (firstContact - 1) / fps;
        end
        turnCount = NaN;
        headingAbsChange = NaN;
        if numel(validHeading) > 1
            dHeading = mod(diff(validHeading) + 180, 360) - 180;
            turnCount = sum(abs(dHeading) > 45);
            headingAbsChange = mean(abs(dHeading), 'omitnan');
        end
        radialSignChanges = NaN;
        if numel(validRadial) > 1
            radialSignChanges = sum(abs(diff(sign(validRadial))) > 0);
        end
        distDelta = NaN;
        if ~isempty(validDist)
            third = max(1, floor(numel(validDist) / 3));
            earlyDist = validDist(1:third);
            lateDist = validDist((numel(validDist) - third + 1):numel(validDist));
            distDelta = mean(lateDist, 'omitnan') - mean(earlyDist, 'omitnan');
        end

        row = table();
        row.scale_index = s;
        row.chunk_sec = scaleSec(s);
        row.temporal_band = anchorRows.temporal_band(s);
        row.anchor_frame = f;
        row.anchor_time_s = (f - 1) / fps;
        row.event_valid_fraction = nnz(valid) ./ max(numel(idx), 1);
        row.contact_dwell_fraction = mean(validContact, 'omitnan');
        row.contact_transition_count = contactTransitions;
        row.first_contact_latency_s = firstContactLatency;
        row.mutual_approach_dwell_fraction = mean(validApproach, 'omitnan');
        row.withdrawal_dwell_fraction = mean(validWithdraw, 'omitnan');
        row.approach_withdraw_transition_count = approachWithdrawTransitions;
        row.role_asymmetry_bias_mean = mean(validAsym, 'omitnan');
        row.centroid_distance_mean = mean(validDist, 'omitnan');
        row.centroid_distance_delta = distDelta;
        row.radial_speed_mean = mean(validRadial, 'omitnan');
        row.radial_speed_sign_change_count = radialSignChanges;
        row.mutual_facing_mean = mean(validFacing, 'omitnan');
        row.heading_abs_change_deg = headingAbsChange;
        row.heading_large_turn_count = turnCount;
        row.event_summary_rule = "feature_thresholds_plus_frame_mask_no_labels";
        row.toy_state_used_for_event_summary = false;
        row.provenance_label_used_for_event_summary = false;
        demoEventRows = [demoEventRows; row]; %#ok<AGROW>
    end
end

demoEventByScale = table();
for s = 1:numel(scaleSec)
    idx = demoEventRows.scale_index == s;
    if ~any(idx)
        continue
    end
    row = table();
    row.scale_index = s;
    row.chunk_sec = scaleSec(s);
    row.temporal_band = anchorRows.temporal_band(s);
    row.n_event_chunks = nnz(idx);
    row.median_valid_fraction = median(demoEventRows.event_valid_fraction(idx), 'omitnan');
    row.median_contact_dwell_fraction = median(demoEventRows.contact_dwell_fraction(idx), 'omitnan');
    row.median_contact_transition_count = median(demoEventRows.contact_transition_count(idx), 'omitnan');
    row.median_mutual_approach_dwell_fraction = median(demoEventRows.mutual_approach_dwell_fraction(idx), 'omitnan');
    row.median_withdrawal_dwell_fraction = median(demoEventRows.withdrawal_dwell_fraction(idx), 'omitnan');
    row.median_approach_withdraw_transition_count = median(demoEventRows.approach_withdraw_transition_count(idx), 'omitnan');
    row.median_heading_large_turn_count = median(demoEventRows.heading_large_turn_count(idx), 'omitnan');
    row.median_role_asymmetry_bias_mean = median(demoEventRows.role_asymmetry_bias_mean(idx), 'omitnan');
    row.labels_used_for_event_summary = "none";
    demoEventByScale = [demoEventByScale; row]; %#ok<AGROW>
end

principleMap = table();
principleMap.demo_stage = [
    "1 synthetic frame features"
    "2 valid-frame transform fit"
    "3 scale-specific anchor eligibility"
    "4 common-anchor cross-scale audit"
    "5 band-adaptive multiresolution summaries"
    "6 PCA dimensionality audit"
    "7 score-selected to primary-scale promotion"
    "8 event-summary audit before embedding"
    ];
principleMap.production_run06_analog = [
    "run_05 dyadic feature MAT files and feature_dictionary.csv"
    "chunk_feature_transform_audit.csv"
    "chunk_anchor_manifest.csv and scale_anchor_coverage_audit.csv"
    "scale usefulness scores from aligned common anchors"
    "band-adaptive summarize_multiresolution_chunks and embedding_dimension_audit.csv"
    "scale_usefulness_scores.csv, pca_loading_stability.csv, selected_operational_scales.csv"
    "primary_operational_scales.csv and primary_scale_specific_anchor_manifest.csv"
    "primary_chunk_event_summary_audit.csv"
    ];
principleMap.labels_allowed = [
    "toy labels only for plotting"
    "none"
    "none"
    "none"
    "none"
    "none"
    "none for selection; provenance only in output tables"
    "none for summaries; provenance only in output tables"
    ];
principleMap.reviewer_question_answered = [
    "What behavior-like signals enter the chunk layer?"
    "Were transforms fit without condition labels and only on valid frames?"
    "Do frame masks and finite features control anchor eligibility?"
    "How are cross-scale scores compared without changing anchor identity?"
    "Why does motif time get richer summaries than micro/context time?"
    "Why is a fixed five-PC representation insufficient?"
    "How do 16 score-selected scales become 13 primary scales?"
    "What social-event detail remains visible before embedding?"
    ];
principleMap.demo_source_csv = [
    "demo_latent_state_segments_visual_only.csv"
    "demo_condition_blind_transform_audit.csv"
    "demo_anchor_coverage.csv"
    "demo_anchor_coverage.csv"
    "demo_embedding_dimension_audit.csv"
    "demo_scale_decision_audit.csv"
    "demo_explanatory_scales.csv"
    "demo_event_summary_by_scale.csv"
    ];

summaryProfileRows = table();
summaryProfileRows.temporal_band = ["micro"; "motif"; "context"];
summaryProfileRows.summary_temporal_bins = [summaryMicroTemporalBins; summaryMotifTemporalBins; summaryContextTemporalBins];
summaryProfileRows.summary_dct_coeffs = [summaryMicroDctCoeffs; summaryMotifDctCoeffs; summaryContextDctCoeffs];
summaryProfileRows.summary_profile = [
    "micro_adaptive_bins6_dct4"
    "motif_adaptive_bins12_dct8"
    "context_adaptive_bins6_dct4"
    ];
summaryProfileRows.n_input_channels = repmat(nFeature, 3, 1);
summaryProfileRows.n_boolean_transition_channels = repmat(nnz(isBoolean), 3, 1);
summaryProfileRows.summary_dims_for_demo_channels = nFeature .* ...
    (7 + summaryProfileRows.summary_temporal_bins + summaryProfileRows.summary_dct_coeffs) + ...
    summaryProfileRows.n_boolean_transition_channels;
summaryProfileRows.labels_used_for_profile_definition = repmat("none", 3, 1);
summaryProfileRows.production_run06_analog = repmat(true, 3, 1);

writetable(segments, fullfile(outRoot, 'demo_latent_state_segments_visual_only.csv'));
writetable(transformRows, fullfile(outRoot, 'demo_condition_blind_transform_audit.csv'));
writetable(anchorRows, fullfile(outRoot, 'demo_anchor_coverage.csv'));
writetable(summaryProfileRows, fullfile(outRoot, 'demo_summary_profiles.csv'));
writetable(pcaRows, fullfile(outRoot, 'demo_pca_by_scale.csv'));
writetable(pcaRows, fullfile(outRoot, 'demo_embedding_dimension_audit.csv'));
writetable(chunkExamples, fullfile(outRoot, 'demo_example_chunk_traces.csv'));
writetable(selectedDemoScales, fullfile(outRoot, 'demo_explanatory_scales.csv'));
writetable(demoDecisionRows, fullfile(outRoot, 'demo_scale_decision_audit.csv'));
writetable(demoEventRows, fullfile(outRoot, 'demo_event_summary_by_chunk.csv'));
writetable(demoEventByScale, fullfile(outRoot, 'demo_event_summary_by_scale.csv'));
writetable(principleMap, fullfile(outRoot, 'demo_principle_map.csv'));

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [60 60 1700 1100]);
tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
yyaxis left;
plot(t, X(:,1), 'LineWidth', 1.0);
ylabel('Distance-like feature', 'Interpreter', 'none');
yyaxis right;
stairs(t, stateId, 'LineWidth', 0.8);
hold on;
plot(t, double(validMask) .* max(stateId), 'k-', 'LineWidth', 1.0);
ylim([0.5 numel(stateNames)+0.5]);
yticks(1:numel(stateNames));
yticklabels(cellstr(stateNames));
xlabel('Time (s)', 'Interpreter', 'none');
title('Toy dyadic phases and frame-mask gaps', 'Interpreter', 'none');

nexttile;
yyaxis left;
plot(t, X(:,1), 'Color', [0.5 0.5 0.5], 'LineWidth', 0.9);
ylabel('Raw distance', 'Interpreter', 'none');
yyaxis right;
plot(t, Xtrans(:,1), 'LineWidth', 1.0);
ylabel('Transformed distance', 'Interpreter', 'none');
xlabel('Time (s)', 'Interpreter', 'none');
title('Condition-blind robust transform', 'Interpreter', 'none');
legend({'Raw distance','Transformed distance'}, 'Box', 'off', 'Location', 'best');

nexttile;
yyaxis left;
plot(anchorRows.chunk_sec, anchorRows.n_scale_specific_candidate_anchors, 'o-', 'LineWidth', 1.4);
hold on;
plot(anchorRows.chunk_sec, anchorRows.n_common_all_scale_candidate_anchors, 's--', 'LineWidth', 1.2);
ylabel('Anchor count', 'Interpreter', 'none');
yyaxis right;
plot(anchorRows.chunk_sec, anchorRows.common_anchor_retention_fraction, 'd:', 'LineWidth', 1.2);
ylabel('Common-anchor retention', 'Interpreter', 'none');
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
title('Scale-specific versus common anchors', 'Interpreter', 'none');
legend({'Scale-specific','Common to all scales','Retention'}, 'Box', 'off', 'Location', 'best');

nexttile;
plot(anchorRows.chunk_sec, anchorRows.n_candidates_removed_by_frame_mask, 'o-', 'LineWidth', 1.4);
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Rejected candidates', 'Interpreter', 'none');
title('Frame-mask propagation grows with context length', 'Interpreter', 'none');

nexttile;
yyaxis left;
plot(pcaRows.chunk_sec, pcaRows.raw_flattened_dimensions, 'o-', 'LineWidth', 1.4);
hold on;
plot(pcaRows.chunk_sec, pcaRows.summary_dimensions, 's-', 'LineWidth', 1.4);
ylabel('Representation dimensions', 'Interpreter', 'none');
set(gca, 'YScale', 'log');
yyaxis right;
plot(pcaRows.chunk_sec, pcaRows.n_pcs_90pct, 'd--', 'LineWidth', 1.2);
hold on;
plot(pcaRows.chunk_sec, pcaRows.recommended_pcs_for_next_embedding, '^-', 'LineWidth', 1.2);
plot(pcaRows.chunk_sec, repmat(5, height(pcaRows), 1), 'k:', 'LineWidth', 1.0);
ylabel('PCs on summary representation', 'Interpreter', 'none');
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
title('Multiresolution summaries compress long contexts', 'Interpreter', 'none');
legend({'Raw flattened dims','Summary dims','PCs for 90%','Recommended cap','Fixed 5 PCs'}, ...
    'Box', 'off', 'Location', 'best');

nexttile;
axis off;
summaryText = {
    'run_06 principles demonstrated here:'
    '1. Fit transforms on valid frames only, without labels.'
    '2. Build anchors from frame mask, feature availability, time, and scale.'
    '3. Use common anchors for fair scale scoring.'
    '4. Use primary scale-specific anchors for downstream coverage.'
    '5. Give motif windows richer 12-bin/8-DCT summaries.'
    '6. Promote score-selected scales only after stability/dimension guards.'
    '7. Carry provenance only for diagnostics and later statistics.'
    };
text(0, 0.95, summaryText, 'Units', 'normalized', 'Interpreter', 'none', ...
    'VerticalAlignment', 'top', 'FontName', 'Arial', 'FontSize', 10);

exportgraphics(fig, fullfile(figRoot, 'demo_run06_multiscale_walkthrough.png'), 'Resolution', 240);
exportgraphics(fig, fullfile(figRoot, 'demo_run06_multiscale_walkthrough.pdf'));
close(fig);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [70 70 1550 980]);
tiledlayout(fig, numel(exampleScaleSec), numel(featuresToPlot), ...
    'TileSpacing', 'compact', 'Padding', 'compact');
for s = 1:numel(exampleScaleSec)
    [~, scaleIdx] = min(abs(scaleSec - exampleScaleSec(s)));
    for f = 1:numel(featuresToPlot)
        nexttile;
        idx = chunkExamples.scale_index == scaleIdx & chunkExamples.feature_name == featuresToPlot(f);
        T = chunkExamples(idx, :);
        plot(T.relative_time_s, T.value, 'LineWidth', 1.0);
        hold on;
        bad = ~logical(T.valid_mask);
        if any(bad)
            scatter(T.relative_time_s(bad), T.value(bad), 10, 'filled');
        end
        title(sprintf('%s, %.1fs', featuresToPlot(f), scaleSec(scaleIdx)), 'Interpreter', 'none');
        if f == 1
            ylabel(sprintf('Scale %.1fs', scaleSec(scaleIdx)), 'Interpreter', 'none');
        end
        if s == numel(exampleScaleSec)
            xlabel('Time from anchor (s)', 'Interpreter', 'none');
        end
    end
end
exportgraphics(fig, fullfile(figRoot, 'demo_run06_example_chunks.png'), 'Resolution', 240);
exportgraphics(fig, fullfile(figRoot, 'demo_run06_example_chunks.pdf'));
close(fig);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1700 1050]);
tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(pcaRows.chunk_sec, pcaRows.raw_flattened_dimensions, 'o-', 'LineWidth', 1.3);
hold on;
plot(pcaRows.chunk_sec, pcaRows.summary_dimensions, 's-', 'LineWidth', 1.3);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Dimensions', 'Interpreter', 'none');
title('Raw long-window flattening is not the analysis layer', 'Interpreter', 'none');
legend({'Raw frame x channel','Multiresolution summary'}, 'Box', 'off', 'Location', 'best');

nexttile;
plot(pcaRows.chunk_sec, pcaRows.cum5_explained, 'o-', 'LineWidth', 1.3);
hold on;
plot(pcaRows.chunk_sec, pcaRows.cum12_explained, 's-', 'LineWidth', 1.3);
yline(fixedFivePcVarianceTargetPct, 'k--', 'LineWidth', 1.0);
set(gca, 'XScale', 'log');
ylim([0 100]);
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Variance explained (%)', 'Interpreter', 'none');
title('Fixed five PCs is an audit failure, not a default', 'Interpreter', 'none');
legend({'First 5 PCs','Scoring PCs','80% reference'}, 'Box', 'off', 'Location', 'best');

nexttile;
plot(pcaRows.chunk_sec, pcaRows.n_pcs_90pct, 'o-', 'LineWidth', 1.3);
hold on;
plot(pcaRows.chunk_sec, pcaRows.recommended_pcs_for_next_embedding, 's-', 'LineWidth', 1.3);
plot(pcaRows.chunk_sec, repmat(scorePcsRetained, height(pcaRows), 1), 'k:', 'LineWidth', 1.0);
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('PCA components', 'Interpreter', 'none');
title('Scale-specific dimensionality should guide embedding', 'Interpreter', 'none');
legend({'PCs for 90%','Recommended cap','Scoring PCs'}, 'Box', 'off', 'Location', 'best');

nexttile;
plot(demoEventByScale.chunk_sec, demoEventByScale.median_contact_dwell_fraction, 'o-', 'LineWidth', 1.3);
hold on;
plot(demoEventByScale.chunk_sec, demoEventByScale.median_mutual_approach_dwell_fraction, 's-', 'LineWidth', 1.3);
plot(demoEventByScale.chunk_sec, demoEventByScale.median_withdrawal_dwell_fraction, 'd-', 'LineWidth', 1.3);
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Median dwell fraction', 'Interpreter', 'none');
title('Event dwell summaries preserve social content', 'Interpreter', 'none');
legend({'Contact','Mutual approach','Withdrawal'}, 'Box', 'off', 'Location', 'best');

nexttile;
plot(demoEventByScale.chunk_sec, demoEventByScale.median_contact_transition_count, 'o-', 'LineWidth', 1.3);
hold on;
plot(demoEventByScale.chunk_sec, demoEventByScale.median_approach_withdraw_transition_count, 's-', 'LineWidth', 1.3);
plot(demoEventByScale.chunk_sec, demoEventByScale.median_heading_large_turn_count, 'd-', 'LineWidth', 1.3);
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Median event count', 'Interpreter', 'none');
title('Longer windows capture transitions, not just state averages', 'Interpreter', 'none');
legend({'Contact changes','Approach/withdraw changes','Large turns'}, 'Box', 'off', 'Location', 'best');

nexttile;
axis off;
detailText = {
    'How to read the production step:'
    'A. Scale scores survey the space; they are not biological tests.'
    'B. Primary scales define the downstream bank after predefined guards.'
    'C. Scale-specific anchors feed embedding.'
    'D. Common anchors stay as cross-scale audit rows.'
    'E. Event summaries ask whether social detail survives.'
    'F. Later clustering should use scale-specific PC caps.'
    'G. Arena/QC diagnostics remain sensitivity audits.'
    };
text(0, 0.95, detailText, 'Units', 'normalized', 'Interpreter', 'none', ...
    'VerticalAlignment', 'top', 'FontName', 'Arial', 'FontSize', 9);

exportgraphics(fig, fullfile(figRoot, 'demo_run06_scale_decision_event_detail.png'), 'Resolution', 240);
exportgraphics(fig, fullfile(figRoot, 'demo_run06_scale_decision_event_detail.pdf'));
close(fig);

figureManifest = table();
figureManifest.figure_id = [
    "demo_run06_multiscale_walkthrough"
    "demo_run06_example_chunks"
    "demo_run06_scale_decision_event_detail"
    ];
figureManifest.figure_file = [
    string(fullfile(figRoot, 'demo_run06_multiscale_walkthrough.png'))
    string(fullfile(figRoot, 'demo_run06_example_chunks.png'))
    string(fullfile(figRoot, 'demo_run06_scale_decision_event_detail.png'))
    ];
figureManifest.source_csv = [
    string(fullfile(outRoot, 'demo_anchor_coverage.csv')) + "; " + ...
        string(fullfile(outRoot, 'demo_embedding_dimension_audit.csv')) + "; " + ...
        string(fullfile(outRoot, 'demo_summary_profiles.csv')) + "; " + ...
        string(fullfile(outRoot, 'demo_condition_blind_transform_audit.csv'))
    string(fullfile(outRoot, 'demo_example_chunk_traces.csv'))
    string(fullfile(outRoot, 'demo_scale_decision_audit.csv')) + "; " + ...
        string(fullfile(outRoot, 'demo_event_summary_by_scale.csv')) + "; " + ...
        string(fullfile(outRoot, 'demo_principle_map.csv'))
    ];
figureManifest.principle = [
    "condition-blind transforms, two anchor banks, band-adaptive summaries, and PC dimensionality audits"
    "same anchor viewed at primary micro, primary motif, primary context, and rejected long-context scales"
    "why scale-specific dimensionality and event-detail audits should guide the next embedding step"
    ];
figureManifest.result_interpretation = [
    "The toy example shows that frame-mask gaps increasingly constrain long windows, while the current band-adaptive summaries keep representations auditable: micro/context stay compact and motif windows receive richer 12-bin/8-DCT temporal detail."
    "The same anchor contains different information at different scales: short windows emphasize instantaneous posture/contact, motif windows preserve local approach-contact-withdraw order, primary context windows capture local persistence, and long context windows are useful audits but not automatically primary."
    "The detail audit separates four ideas that should not be conflated: raw flattened dimensionality, band-adaptive summary dimensionality, PCA dimensionality, and social event information. This supports using run_06 as a scale survey followed by scale-specific embedding rather than a single fixed-PC representation."
    ];
figureManifest.toy_labels_used_for_fitting = [false; false; false];
writetable(figureManifest, fullfile(outRoot, 'demo_figure_manifest.csv'));

stepFigureIds = strings(9, 1);
stepFigureFiles = strings(9, 1);
stepSourceCsv = strings(9, 1);
stepPrinciples = strings(9, 1);
stepInterpretations = strings(9, 1);

stepFigureIds(1) = "demo_run06_step01_feature_traces_valid_mask";
stepFigureFiles(1) = string(fullfile(figRoot, stepFigureIds(1) + ".png"));
stepSourceCsv(1) = "demo_latent_state_segments_visual_only.csv";
stepPrinciples(1) = "Start with time-series features and a frame-validity mask.";
stepInterpretations(1) = "The same dyadic session contains feature dynamics and invalid frame gaps. Run_06 treats the mask as part of the data quality model, so chunk eligibility depends on whether enough valid frames surround each candidate anchor.";

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1180 660]);
yyaxis left;
plot(t, X(:,1), 'Color', [0.1 0.35 0.75], 'LineWidth', 1.0);
ylabel('Centroid distance', 'Interpreter', 'none');
yyaxis right;
stairs(t, stateId, 'Color', [0.55 0.55 0.55], 'LineWidth', 0.8);
hold on;
plot(t, double(validMask) .* max(stateId), 'k-', 'LineWidth', 1.0);
ylim([0.5 numel(stateNames)+0.5]);
yticks(1:numel(stateNames));
yticklabels(cellstr(stateNames));
xlabel('Time (s)', 'Interpreter', 'none');
title('Step 1: feature traces plus valid-frame mask', 'Interpreter', 'none');
legend({'Distance feature','Toy phase label for display','Valid-frame mask'}, ...
    'Box', 'off', 'Location', 'southoutside', 'Orientation', 'horizontal');
exportgraphics(fig, fullfile(figRoot, stepFigureIds(1) + ".png"), 'Resolution', 240);
exportgraphics(fig, fullfile(figRoot, stepFigureIds(1) + ".pdf"));
close(fig);

stepFigureIds(2) = "demo_run06_step02_condition_blind_transform";
stepFigureFiles(2) = string(fullfile(figRoot, stepFigureIds(2) + ".png"));
stepSourceCsv(2) = "demo_condition_blind_transform_audit.csv";
stepPrinciples(2) = "Fit robust feature transforms without labels.";
stepInterpretations(2) = "Robust centering and scaling put heterogeneous dyadic features on comparable units. The toy labels are not used for the fit, matching the condition-blind production rule.";

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1050 520]);
yyaxis left;
plot(t, X(:,1), 'Color', [0.65 0.65 0.65], 'LineWidth', 0.9);
ylabel('Raw distance', 'Interpreter', 'none');
yyaxis right;
plot(t, Xtrans(:,1), 'Color', [0.1 0.45 0.35], 'LineWidth', 1.0);
ylabel('Transformed distance', 'Interpreter', 'none');
xlabel('Time (s)', 'Interpreter', 'none');
title('Step 2: condition-blind robust transform', 'Interpreter', 'none');
legend({'Raw feature','Transformed feature'}, 'Box', 'off', 'Location', 'best');
exportgraphics(fig, fullfile(figRoot, stepFigureIds(2) + ".png"), 'Resolution', 240);
exportgraphics(fig, fullfile(figRoot, stepFigureIds(2) + ".pdf"));
close(fig);

stepFigureIds(3) = "demo_run06_step03_anchor_eligibility";
stepFigureFiles(3) = string(fullfile(figRoot, stepFigureIds(3) + ".png"));
stepSourceCsv(3) = "demo_anchor_coverage.csv";
stepPrinciples(3) = "Candidate anchors survive only if their surrounding chunk is valid.";
stepInterpretations(3) = "For a fixed stride, candidate anchor times are filtered by anchor-frame validity, valid-frame fraction, and finite transformed features. Longer chunks are more likely to overlap gaps, so this gate matters before any PCA or scale scoring.";

scaleForAnchorDemo = find(abs(scaleSec - 5.8420) < 1e-6, 1);
chunkFramesDemo = max(1, round(scaleSec(scaleForAnchorDemo) * fps));
leftFramesDemo = chunkFramesDemo - 1;
candidateFramesDemo = (1 + leftFramesDemo):strideFrames:nFrames;
keepDemo = false(numel(candidateFramesDemo), 1);
validFracDemo = nan(numel(candidateFramesDemo), 1);
for i = 1:numel(candidateFramesDemo)
    idx = (candidateFramesDemo(i)-leftFramesDemo):candidateFramesDemo(i);
    validFracDemo(i) = mean(validMask(idx));
    keepDemo(i) = validMask(candidateFramesDemo(i)) && validFracDemo(i) >= minValidFrac && ...
        mean(isfinite(Xtrans(idx,:)), 'all') >= minFeatureFiniteFrac;
end

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1050 520]);
plot(t, Xtrans(:,1), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
hold on;
scatter(t(candidateFramesDemo(~keepDemo)), Xtrans(candidateFramesDemo(~keepDemo),1), ...
    12, [0.75 0.25 0.2], 'filled');
scatter(t(candidateFramesDemo(keepDemo)), Xtrans(candidateFramesDemo(keepDemo),1), ...
    18, [0.1 0.45 0.75], 'filled');
xlabel('Time (s)', 'Interpreter', 'none');
ylabel('Transformed distance', 'Interpreter', 'none');
title(sprintf('Step 3: anchor eligibility at %.1f s', scaleSec(scaleForAnchorDemo)), ...
    'Interpreter', 'none');
legend({'Feature trace','Rejected candidates','Eligible anchors'}, ...
    'Box', 'off', 'Location', 'best');
exportgraphics(fig, fullfile(figRoot, stepFigureIds(3) + ".png"), 'Resolution', 240);
exportgraphics(fig, fullfile(figRoot, stepFigureIds(3) + ".pdf"));
close(fig);

stepFigureIds(4) = "demo_run06_step04_common_vs_scale_specific_anchors";
stepFigureFiles(4) = string(fullfile(figRoot, stepFigureIds(4) + ".png"));
stepSourceCsv(4) = "demo_anchor_coverage.csv";
stepPrinciples(4) = "Use common anchors for fair scale comparison and scale-specific anchors for coverage.";
stepInterpretations(4) = "Common all-scale anchors are the subset valid at every candidate scale, which makes the scale survey fair but conservative. Scale-specific anchors preserve more valid material for downstream embedding after primary scales have been chosen.";

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1050 520]);
yyaxis left;
plot(anchorRows.chunk_sec, anchorRows.n_scale_specific_candidate_anchors, 'o-', 'LineWidth', 1.4);
hold on;
plot(anchorRows.chunk_sec, anchorRows.n_common_all_scale_candidate_anchors, 's--', 'LineWidth', 1.3);
ylabel('Anchor count', 'Interpreter', 'none');
yyaxis right;
plot(anchorRows.chunk_sec, anchorRows.common_anchor_retention_fraction, 'd:', 'LineWidth', 1.3);
ylabel('Common-anchor retention', 'Interpreter', 'none');
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
title('Step 4: common anchors versus scale-specific anchors', 'Interpreter', 'none');
legend({'Scale-specific candidates','Common all-scale anchors','Retention'}, ...
    'Box', 'off', 'Location', 'best');
exportgraphics(fig, fullfile(figRoot, stepFigureIds(4) + ".png"), 'Resolution', 240);
exportgraphics(fig, fullfile(figRoot, stepFigureIds(4) + ".pdf"));
close(fig);

stepFigureIds(5) = "demo_run06_step05_multiresolution_summary";
stepFigureFiles(5) = string(fullfile(figRoot, stepFigureIds(5) + ".png"));
stepSourceCsv(5) = "demo_example_chunk_traces.csv";
stepPrinciples(5) = "Use band-adaptive summaries: richer temporal detail for motif windows.";
stepInterpretations(5) = "The motif-band example uses 12 temporal bins and 8 low-frequency DCT coefficients, matching the current run_06 setting. This preserves local order in approach-contact-withdraw fragments without storing every frame as a separate dimension.";

scaleForSummaryDemo = find(abs(scaleSec - 1.2601) < 1e-6, 1);
Tsum = chunkExamples(chunkExamples.scale_index == scaleForSummaryDemo & ...
    chunkExamples.feature_name == "centroid_dist", :);
xRel = Tsum.relative_time_s;
yTrace = Tsum.value;
if any(~isfinite(yTrace))
    yTrace(~isfinite(yTrace)) = median(yTrace, 'omitnan');
end
demoSummaryTemporalBins = summaryMotifTemporalBins;
demoSummaryDctCoeffs = summaryMotifDctCoeffs;
binEdgesDemo = unique(round(linspace(1, numel(yTrace) + 1, demoSummaryTemporalBins + 1)));
binMeanX = nan(demoSummaryTemporalBins, 1);
binMeanY = nan(demoSummaryTemporalBins, 1);
for b = 1:demoSummaryTemporalBins
    idx = binEdgesDemo(b):(binEdgesDemo(b+1)-1);
    idx = idx(idx >= 1 & idx <= numel(yTrace));
    binMeanX(b) = mean(xRel(idx), 'omitnan');
    binMeanY(b) = mean(yTrace(idx), 'omitnan');
end
dctBasisDemo = zeros(demoSummaryDctCoeffs, numel(yTrace));
mDemo = 0:(numel(yTrace)-1);
for k = 1:demoSummaryDctCoeffs
    dctBasisDemo(k,:) = cos(pi * (mDemo + 0.5) * (k - 1) / numel(yTrace));
end
dctBasisDemo(1,:) = dctBasisDemo(1,:) ./ sqrt(numel(yTrace));
if demoSummaryDctCoeffs > 1
    dctBasisDemo(2:end,:) = dctBasisDemo(2:end,:) .* sqrt(2 / numel(yTrace));
end
yCentered = yTrace(:)' - mean(yTrace, 'omitnan');
dctCoeffDemo = yCentered * dctBasisDemo';
ySmooth = dctCoeffDemo * dctBasisDemo + mean(yTrace, 'omitnan');

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1050 520]);
plot(xRel, yTrace, 'Color', [0.35 0.35 0.35], 'LineWidth', 1.0);
hold on;
plot(xRel, ySmooth(:), 'Color', [0.05 0.45 0.70], 'LineWidth', 1.8);
scatter(binMeanX, binMeanY, 55, [0.85 0.45 0.10], 'filled');
xlabel('Time from anchor (s)', 'Interpreter', 'none');
ylabel('Transformed distance', 'Interpreter', 'none');
title('Step 5: motif-band 12-bin/8-DCT summary', 'Interpreter', 'none');
legend({'Frame trace','Low-frequency DCT shape','Temporal-bin means'}, ...
    'Box', 'off', 'Location', 'best');
exportgraphics(fig, fullfile(figRoot, stepFigureIds(5) + ".png"), 'Resolution', 240);
exportgraphics(fig, fullfile(figRoot, stepFigureIds(5) + ".pdf"));
close(fig);

stepFigureIds(6) = "demo_run06_step06_dimensionality_compression";
stepFigureFiles(6) = string(fullfile(figRoot, stepFigureIds(6) + ".png"));
stepSourceCsv(6) = "demo_embedding_dimension_audit.csv";
stepPrinciples(6) = "Long-window raw flattening grows quickly, while summaries stay compact.";
stepInterpretations(6) = "Raw frame-by-channel dimensionality increases with scale. The summary representation stays bounded, with a deliberate motif-band bump because motif windows carry 12 bins and 8 DCT coefficients instead of the micro/context 6-bin/4-DCT profile.";

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1050 520]);
plot(pcaRows.chunk_sec, pcaRows.raw_flattened_dimensions, 'o-', 'LineWidth', 1.4);
hold on;
plot(pcaRows.chunk_sec, pcaRows.summary_dimensions, 's-', 'LineWidth', 1.4);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
ylabel('Dimensions', 'Interpreter', 'none');
title('Step 6: compact summaries prevent dimensionality blow-up', 'Interpreter', 'none');
legend({'Raw frame x channel vector','Multiresolution summary'}, ...
    'Box', 'off', 'Location', 'northwest');
exportgraphics(fig, fullfile(figRoot, stepFigureIds(6) + ".png"), 'Resolution', 240);
exportgraphics(fig, fullfile(figRoot, stepFigureIds(6) + ".pdf"));
close(fig);

stepFigureIds(7) = "demo_run06_step07_scale_specific_pca";
stepFigureFiles(7) = string(fullfile(figRoot, stepFigureIds(7) + ".png"));
stepSourceCsv(7) = "demo_embedding_dimension_audit.csv";
stepPrinciples(7) = "Audit PCA dimensionality by scale instead of using a fixed five-PC rule.";
stepInterpretations(7) = "The number of PCs needed to represent each scale changes with temporal context. A fixed five-PC representation can be too small, so run_06 records scale-specific dimensionality guidance for run_07.";

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1050 520]);
yyaxis left;
plot(pcaRows.chunk_sec, pcaRows.cum5_explained, 'o-', 'LineWidth', 1.4);
hold on;
plot(pcaRows.chunk_sec, pcaRows.cum12_explained, 's-', 'LineWidth', 1.4);
yline(fixedFivePcVarianceTargetPct, 'k--', 'LineWidth', 1.0);
ylabel('Variance explained (%)', 'Interpreter', 'none');
yyaxis right;
plot(pcaRows.chunk_sec, pcaRows.n_pcs_90pct, 'd-', 'LineWidth', 1.2);
ylabel('PCs for 90% variance', 'Interpreter', 'none');
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
title('Step 7: scale-specific PCA audit', 'Interpreter', 'none');
legend({'First 5 PCs','Scoring PCs','80% reference','PCs for 90%'}, ...
    'Box', 'off', 'Location', 'best');
exportgraphics(fig, fullfile(figRoot, stepFigureIds(7) + ".png"), 'Resolution', 240);
exportgraphics(fig, fullfile(figRoot, stepFigureIds(7) + ".pdf"));
close(fig);

stepFigureIds(8) = "demo_run06_step08_scale_decision_guards";
stepFigureFiles(8) = string(fullfile(figRoot, stepFigureIds(8) + ".png"));
stepSourceCsv(8) = "demo_scale_decision_audit.csv";
stepPrinciples(8) = "Promote scales only after predefined dimensionality guards.";
stepInterpretations(8) = "High apparent scale scores are not enough. The demo mirrors the production idea that score-selected scales can be rejected before primary promotion if stability or dimensionality guards are not satisfied.";

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1050 520]);
yyaxis left;
plot(demoDecisionRows.chunk_sec, demoDecisionRows.effective_dim, 'o-', 'LineWidth', 1.4);
hold on;
primaryMask = logical(demoDecisionRows.production_primary_scale_analog);
scoreRejectedMask = logical(demoDecisionRows.production_score_selected_scale_analog) & ~primaryMask;
scatter(demoDecisionRows.chunk_sec(primaryMask), demoDecisionRows.effective_dim(primaryMask), ...
    55, [0.0 0.55 0.35], 'filled');
scatter(demoDecisionRows.chunk_sec(scoreRejectedMask), demoDecisionRows.effective_dim(scoreRejectedMask), ...
    55, [0.35 0.35 0.35], 'filled');
yline(dimensionGuardMinEffectiveDim, 'k--', 'LineWidth', 1.0);
ylabel('Effective dimension', 'Interpreter', 'none');
ylim([0 max(8, 1.15 * max(demoDecisionRows.effective_dim, [], 'omitnan'))]);
yyaxis right;
plot(demoDecisionRows.chunk_sec, demoDecisionRows.pc1_explained, 's-', 'LineWidth', 1.4);
yline(dimensionGuardMaxPc1Pct, 'r--', 'LineWidth', 1.0);
ylabel('PC1 variance explained (%)', 'Interpreter', 'none');
ylim([0 85]);
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
title('Step 8: 16 score-selected analogs become 13 primary analogs', 'Interpreter', 'none');
lgd = legend({'Eff. dimension','Primary analog','Rejected selected analog', ...
    'Min eff. dimension','PC1 variance','Max PC1 variance'}, ...
    'Box', 'off', 'Location', 'southoutside', 'Orientation', 'horizontal');
if isprop(lgd, 'NumColumns')
    lgd.NumColumns = 3;
end
exportgraphics(fig, fullfile(figRoot, stepFigureIds(8) + ".png"), 'Resolution', 240);
exportgraphics(fig, fullfile(figRoot, stepFigureIds(8) + ".pdf"));
close(fig);

stepFigureIds(9) = "demo_run06_step09_event_detail_audit";
stepFigureFiles(9) = string(fullfile(figRoot, stepFigureIds(9) + ".png"));
stepSourceCsv(9) = "demo_event_summary_by_scale.csv";
stepPrinciples(9) = "Audit whether social-event detail remains visible before embedding.";
stepInterpretations(9) = "Dwell and transition summaries help explain what social information different scales can carry. They are detail audits at run_06, not condition effects or final motif definitions.";

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1050 520]);
yyaxis left;
plot(demoEventByScale.chunk_sec, demoEventByScale.median_contact_dwell_fraction, 'o-', ...
    'Color', [0.05 0.35 0.75], 'LineWidth', 1.3);
hold on;
plot(demoEventByScale.chunk_sec, demoEventByScale.median_mutual_approach_dwell_fraction, 's-', ...
    'Color', [0.1 0.55 0.35], 'LineWidth', 1.3);
plot(demoEventByScale.chunk_sec, demoEventByScale.median_withdrawal_dwell_fraction, 'd-', ...
    'Color', [0.55 0.25 0.65], 'LineWidth', 1.3);
ylabel('Median dwell fraction', 'Interpreter', 'none');
ylim([0 0.22]);
yyaxis right;
plot(demoEventByScale.chunk_sec, demoEventByScale.median_heading_large_turn_count, '^-', ...
    'Color', [0.85 0.30 0.05], 'LineWidth', 1.2);
ylabel('Median large-turn count', 'Interpreter', 'none');
ylim([0 26]);
set(gca, 'XScale', 'log');
xlabel('Chunk scale (s)', 'Interpreter', 'none');
title('Step 9: social-event detail audit', 'Interpreter', 'none');
legend({'Contact dwell','Approach dwell','Withdrawal dwell','Large turns'}, ...
    'Box', 'off', 'Location', 'southoutside', 'Orientation', 'horizontal');
exportgraphics(fig, fullfile(figRoot, stepFigureIds(9) + ".png"), 'Resolution', 240);
exportgraphics(fig, fullfile(figRoot, stepFigureIds(9) + ".pdf"));
close(fig);

stepFigureManifest = table();
stepFigureManifest.figure_id = stepFigureIds;
stepFigureManifest.figure_file = stepFigureFiles;
stepFigureManifest.source_csv = stepSourceCsv;
stepFigureManifest.principle = stepPrinciples;
stepFigureManifest.result_interpretation = stepInterpretations;
stepFigureManifest.toy_labels_used_for_fitting = false(height(stepFigureManifest), 1);
writetable(stepFigureManifest, fullfile(outRoot, 'demo_stepwise_figure_manifest.csv'));

stepCaptions = [
    "# run_06 Stepwise Demo Captions"
    ""
    "These figures are synthetic teaching figures. They do not use project data. Toy phase labels appear only to orient the reader and are not used for transforms, anchor selection, PCA, scale scoring, or event audits."
    ""
    "## Step 1: Feature traces plus valid-frame mask"
    "A dyadic feature trace is shown together with toy phase labels and the valid-frame mask. The key idea is that the data are not just numbers through time; each time point also has a quality state. A chunk can only be trusted if enough frames around its anchor are valid."
    ""
    "## Step 2: Condition-blind robust transform"
    "Raw features can have very different units. Robust centering and scaling put them on comparable units using valid frames only. Labels are deliberately not used, which protects later unsupervised discovery from learning the experimental design."
    ""
    "## Step 3: Anchor eligibility"
    "Candidate anchors are placed on a fixed time grid and then filtered. A candidate survives only when the anchor frame is valid, the surrounding window has enough valid frames, and the transformed features are finite."
    ""
    "## Step 4: Common versus scale-specific anchors"
    "Common all-scale anchors are valid for every candidate scale, so they are fair for scale comparison. Scale-specific anchors keep more valid material once primary scales are known. The production run_06 uses both ideas for different purposes."
    ""
    "## Step 5: Band-adaptive multiresolution summary"
    "The figure shows a motif-scale chunk as a trace, a smooth low-frequency shape, and 12 temporal-bin means. Motif windows get more bins and DCT coefficients than micro/context windows because the order of approach, contact, and withdrawal matters at this time scale."
    ""
    "## Step 6: Dimensionality compression"
    "Raw frame-by-channel dimensionality grows as windows get longer. Summary dimensions stay compact, with a visible motif-band bump from the 12-bin/8-DCT profile. That is why long-context analysis remains feasible on a local workstation while motif-scale order is represented more richly."
    ""
    "## Step 7: Scale-specific PCA"
    "Different temporal scales need different numbers of PCs. A fixed five-PC rule can miss substantial structure, so run_06 records scale-specific dimensionality guidance for the embedding step."
    ""
    "## Step 8: Decision guards"
    "A scale should not be promoted just because it scores highly. The current production run_06 selected 16 operational scales but promoted 13 primary scales after stability and dimensionality guards. The demo marks the same idea with primary and rejected score-selected analogs."
    ""
    "## Step 9: Event-detail audit"
    "Event dwell and transition summaries show what kinds of social detail each scale can carry. At run_06 these are audits and explanations, not condition effects and not final motif definitions."
    ""
    "## Educational Link"
    "MouseMotionMapper illustrates behavior mapping as staged processing from pose or feature trajectories to dimensionality reduction, time-frequency summaries, embedding, clustering, and fingerprinting. The run_06 demo follows the same teaching style but adapts it to dyadic social features and condition-blind multiscale chunking."
    ""
    "Reference: https://github.com/PrincetonUniversity/MouseMotionMapper"
    ""
    "Scientific anchors:"
    "- Berman et al. introduced unsupervised behavioral space mapping from posture/movement time series: https://arxiv.org/abs/1310.4249"
    "- Klibaite et al. extended unsupervised behavioral repertoire mapping to interacting individuals: https://arxiv.org/abs/1609.09345"
    ];
writelines(stepCaptions, fullfile(outRoot, 'demo_stepwise_figure_captions.md'));

interpretation = [
    "# run_06 Demo Figure Interpretations"
    ""
    "## demo_run06_multiscale_walkthrough"
    "This synthetic walkthrough illustrates the run_06 design without using project data. Toy phase labels are plotted only to orient the reader. The transform panel shows that robust scaling is fit from valid frames only. The anchor panel shows the distinction between the conservative common-anchor bank used for fair scale scoring and the scale-specific bank used downstream. The dimensionality panel shows the current patch logic: micro/context windows use compact 6-bin/4-DCT summaries, motif windows use richer 12-bin/8-DCT summaries, and fixed five-PC summaries are not a defensible default."
    ""
    "## demo_run06_example_chunks"
    "The same anchor is visualized at a primary micro, primary motif, primary context, and rejected long-context scale. The micro window captures instantaneous geometry near the anchor. The motif-scale window captures a local approach/contact fragment. The primary context window captures local persistence. The long-context window captures multiple transitions, which can be useful for audit and interpretation but should not be promoted automatically if it becomes too low-dimensional or unstable."
    ""
    "## demo_run06_scale_decision_event_detail"
    "This figure connects dimensionality, band-adaptive summaries, and behavioral detail. Raw frame-by-channel vectors grow linearly with window length, whereas summaries remain compact with a deliberate motif-band increase. The five-PC reference often falls short of the variance target, so downstream embedding should use scale-specific dimensionality targets. The event panels show that longer chunks preserve transitions and direction changes that are not captured by short-window state averages alone."
    ""
    "## What This Demo Does Not Claim"
    "The toy labels are not used for fitting, anchoring, PCA, or scale choice. The primary and rejected scale analogs are explanatory, not a substitute for production run_06 scale scores, bootstrap stability, or dimension guards. The demo is intended to make the condition-blind logic inspectable before downstream embedding and clustering."
    ];
writelines(interpretation, fullfile(outRoot, 'demo_result_interpretations.md'));

walkthrough = [
    "# run_06 Multiscale Chunking Demo"
    ""
    "This demo is synthetic and does not use project data. Toy phase labels are used only to make the figure understandable; they are never used for transforms, anchors, scale scoring, or PCA."
    ""
    "## Files"
    "- demo_condition_blind_transform_audit.csv: robust transform fit on valid frames only."
    "- demo_anchor_coverage.csv: scale-specific anchors, all-scale common anchors, and retention."
    "- demo_summary_profiles.csv: micro, motif, and context summary profiles used by the current run_06 logic."
    "- demo_pca_by_scale.csv: compact-summary PCA audit by scale."
    "- demo_embedding_dimension_audit.csv: raw flattened dimensions, band-adaptive summary dimensions, compression, and suggested PC caps."
    "- demo_scale_decision_audit.csv: condition-blind audit columns explaining fixed-five-PC insufficiency, dimension guards, primary analogs, rejected score-selected analogs, and illustrative scale choices."
    "- demo_event_summary_by_chunk.csv: event summaries for each synthetic anchor and scale."
    "- demo_event_summary_by_scale.csv: scale-level event dwell and transition summaries."
    "- demo_principle_map.csv: reviewer-facing map from demo panels to production run_06 outputs."
    "- demo_example_chunk_traces.csv: example primary micro, primary motif, primary context, and rejected long-context chunks."
    "- demo_figure_manifest.csv: figure-to-source mapping and brief interpretation for each generated figure."
    "- demo_stepwise_figure_manifest.csv: one-row-per-teaching-step figure manifest."
    "- demo_stepwise_figure_captions.md: high-school-readable captions for the stepwise figures."
    "- figures/demo_run06_multiscale_walkthrough.png: overview of the principles."
    "- figures/demo_run06_example_chunks.png: what the same anchor looks like at multiple scales."
    "- figures/demo_run06_scale_decision_event_detail.png: dimensionality and event-detail audit."
    "- figures/demo_run06_step01_feature_traces_valid_mask.png: feature trace plus frame mask."
    "- figures/demo_run06_step02_condition_blind_transform.png: label-free robust transform."
    "- figures/demo_run06_step03_anchor_eligibility.png: valid and rejected anchor candidates."
    "- figures/demo_run06_step04_common_vs_scale_specific_anchors.png: fair scale comparison versus downstream coverage."
    "- figures/demo_run06_step05_multiresolution_summary.png: motif-band trace, 12 temporal bins, and 8-DCT low-frequency shape."
    "- figures/demo_run06_step06_dimensionality_compression.png: raw flattening versus band-adaptive compact summaries."
    "- figures/demo_run06_step07_scale_specific_pca.png: scale-specific PC audit."
    "- figures/demo_run06_step08_scale_decision_guards.png: dimension guards before scale promotion."
    "- figures/demo_run06_step09_event_detail_audit.png: event-detail audit before embedding."
    "- demo_result_interpretations.md: result-style text for the generated demo figures."
    ""
    "## Interpretation"
    "The demo mirrors the production rule: labels are not used to define the representation. As chunks get longer, frame-mask gaps remove more candidate anchors. The current run_06 representation summarizes each chunk with condition-blind band-adaptive multiresolution statistics before PCA: micro/context scales remain compact, while motif scales receive extra temporal detail. The fixed-five-PC reference line is included to make clear that a small fixed PC count is not a defensible default for longer or richer scales."
    ""
    "## Design Implication"
    "Use run_06 as a condition-blind scale survey. For final embedding, build a selected-scale bank with scale-specific anchors, stability-audited operational scales, and per-scale dimensionality targets, while retaining common-anchor outputs for cross-scale diagnostics."
    ];
writelines(walkthrough, fullfile(outRoot, 'README_run06_demo_walkthrough.md'));

fprintf('Wrote demo walkthrough outputs to %s\n', outRoot);
