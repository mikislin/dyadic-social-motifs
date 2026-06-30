function demofeatures()
%DEMOFEATURES Illustrate the canonical 34 dyadic motif-discovery features.
%
% This demo uses synthetic dyad tracks, not raw experiment data. It exists to
% make the feature stage auditable: the same canonical node metadata,
% canonical feature dictionary, units, bad-frame propagation, and publication
% feature names used by paper/run_05_extract_motif_dyad_features.m are used
% here on a small readable example.

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
cd(repoRoot);
addpath(genpath(repoRoot));

cfg = load_preprocessing_pipeline_config(fullfile(repoRoot, 'config', 'preprocessing_pipeline_config.csv'));
featureParams = load_motif_dyad_feature_extraction_config( ...
    fullfile(repoRoot, 'config', 'motif_dyad_feature_extraction_config.csv'));
[nodeMap, ~, nodeMetadata] = default_sleap_node_map(fullfile(repoRoot, 'config', 'preprocessing_pipeline_config.csv'));
[featureNames, featureMeta] = default_dyad_feature_metadata();

outRoot = fullfile(repoRoot, 'derived', 'features_demo');
figDir = fullfile(outRoot, 'figures');
if ~exist(outRoot, 'dir'); mkdir(outRoot); end
if ~exist(figDir, 'dir'); mkdir(figDir); end

fps = featureParams.fps;
pixelSizeMM = featureParams.pixel_size_mm;
[tracks, badframes, poseAudit] = local_make_synthetic_dyad(nodeMap, height(nodeMetadata), fps);

dyad = compute_dyad_features(tracks, fps, nodeMap, ...
    'badframes', badframes, ...
    'pixelSizeMM', pixelSizeMM, ...
    'contactThresholdMM', featureParams.contact_threshold_mm, ...
    'closeThresholdMM', featureParams.close_threshold_mm, ...
    'smoothSpanFrames', featureParams.smooth_span_frames);

dictionaryPath = fullfile(outRoot, 'demo_feature_dictionary.csv');
familySummaryPath = fullfile(outRoot, 'demo_feature_family_summary.csv');
timeseriesPath = fullfile(outRoot, 'demo_feature_timeseries.csv');
posePath = fullfile(outRoot, 'demo_pose_snapshots.csv');
figureManifestPath = fullfile(outRoot, 'demo_feature_figure_manifest.csv');

writetable(featureMeta, dictionaryPath);
familySummary = local_feature_family_summary(featureMeta);
writetable(familySummary, familySummaryPath);
timeSeries = local_demo_timeseries_table(dyad, featureNames);
writetable(timeSeries, timeseriesPath);
poseSnapshots = local_pose_snapshot_table(poseAudit, nodeMetadata);
writetable(poseSnapshots, posePath);

rows = struct([]);
rows = local_add_manifest_rows(rows, local_plot_demo_geometry(posePath, figDir, cfg), ...
    posePath, 'Synthetic dyad body geometry and frame snapshots used for feature illustration.');
rows = local_add_manifest_rows(rows, local_plot_demo_timeseries(timeseriesPath, figDir, cfg), ...
    timeseriesPath, 'Selected canonical feature time series with bad-frame propagation visible as gaps.');
rows = local_add_manifest_rows(rows, local_plot_demo_family_summary(familySummaryPath, figDir, cfg), ...
    familySummaryPath, 'Feature family counts and rationale for the stable 34-feature dictionary.');
figureManifest = struct2table(rows);
writetable(figureManifest, figureManifestPath);

fprintf('Wrote demo feature dictionary: %s\n', dictionaryPath);
fprintf('Wrote demo feature family summary: %s\n', familySummaryPath);
fprintf('Wrote demo feature time series: %s\n', timeseriesPath);
fprintf('Wrote demo pose snapshots: %s\n', posePath);
fprintf('Wrote demo figure manifest: %s\n', figureManifestPath);
fprintf('Feature count: %d stable dyadic features.\n', numel(featureNames));
disp(familySummary(:, {'FamilyLabel','n_features'}));
end

function [tracks, badframes, poseAudit] = local_make_synthetic_dyad(nodeMap, nNodes, fps)
T = 480;
phase = linspace(0, 1, T)';

c1 = [95 + 230*phase, 185 + 22*sin(2*pi*phase)];
sepX = 280 - 260*exp(-((phase - 0.72) ./ 0.20).^2);
sepY = 45*cos(2*pi*phase + 0.30);
c2 = c1 + [sepX, sepY];

toward12 = c2 - c1;
toward21 = -toward12;
theta1 = atan2d(toward12(:,2), toward12(:,1));
theta2 = atan2d(toward21(:,2), toward21(:,1));

v1Heading = local_velocity_heading(c1);
v2Heading = local_velocity_heading(c2);
late = phase > 0.60;
theta1(late) = v1Heading(late);
theta2(late) = v2Heading(late);

tracks = nan(T, nNodes, 2, 2);
tracks = local_place_mouse(tracks, nodeMap, 1, c1, theta1);
tracks = local_place_mouse(tracks, nodeMap, 2, c2, theta2);

badframes = false(T, 1);
badframes(round(2.65*fps):round(2.85*fps)) = true;
badframes(round(4.25*fps):round(4.35*fps)) = true;

snapshotFrames = round([0.15 0.45 0.75] .* (T-1)) + 1;
poseAudit = struct();
poseAudit.tracks = tracks;
poseAudit.snapshotFrames = snapshotFrames(:);
poseAudit.time_s = (snapshotFrames(:) - 1) ./ fps;
poseAudit.badframes = badframes;
end

function heading = local_velocity_heading(xy)
v = [xy(2,:) - xy(1,:); diff(xy, 1, 1)];
heading = atan2d(v(:,2), v(:,1));
heading = fillmissing(heading, 'nearest');
end

function tracks = local_place_mouse(tracks, nodeMap, animal, centroid, thetaDeg)
u = [cosd(thetaDeg), sind(thetaDeg)];
lat = [-u(:,2), u(:,1)];

points = struct();
points.tailTip = centroid - 42*u;
points.tailMid = centroid - 30*u;
points.tailBase = centroid - 16*u;
points.midBody = centroid;
points.body = centroid + 6*u;
points.neck = centroid + 22*u;
points.nose = centroid + 36*u;
points.leftEar = centroid + 26*u + 7*lat;
points.rightEar = centroid + 26*u - 7*lat;
points.lfPaw = centroid + 11*u + 12*lat;
points.rfPaw = centroid + 11*u - 12*lat;
points.lhPaw = centroid - 8*u + 11*lat;
points.rhPaw = centroid - 8*u - 11*lat;

fields = fieldnames(points);
for i = 1:numel(fields)
    f = fields{i};
    if isfield(nodeMap, f)
        tracks(:, nodeMap.(f), :, animal) = reshape(points.(f), size(points.(f),1), 1, 2, 1);
    end
end
end

function T = local_feature_family_summary(featureMeta)
familyLabels = unique(featureMeta.FamilyLabel, 'stable');
rationale = strings(numel(familyLabels), 1);
for i = 1:numel(familyLabels)
    switch familyLabels(i)
        case "Inter-animal distance"
            rationale(i) = "Scale-aware proximity anchors shared by most social motifs.";
        case "Directed target distance"
            rationale(i) = "Directional nose-to-partner geometry captures investigation asymmetry.";
        case "Egocentric partner position"
            rationale(i) = "Partner location in each animal's body frame separates front/side/back motifs.";
        case "Facing score"
            rationale(i) = "Line-of-sight scores quantify attention-like orientation without condition labels.";
        case "Relative orientation"
            rationale(i) = "Heading alignment and signed angular relation distinguish parallel, opposed, and crossing states.";
        case "Relative velocity"
            rationale(i) = "Radial and tangential movement describe approach, separation, and circling.";
        case "Approach velocity"
            rationale(i) = "Directed closure speeds separate who is moving toward whom.";
        case "Movement coupling"
            rationale(i) = "Velocity and acceleration alignment summarize shared movement dynamics.";
        case "Egocentric bearing"
            rationale(i) = "Bearing angles preserve where partner head/body targets fall in each animal's frame.";
        case "Proximity state"
            rationale(i) = "Thresholded states provide interpretable contact and close-pair occupancy.";
        case "Interaction state"
            rationale(i) = "Simple rule states expose mutual approach and withdrawal motifs.";
        otherwise
            rationale(i) = "Compact signed index for asymmetry in directed investigation.";
    end
end

n = zeros(numel(familyLabels), 1);
featureList = strings(numel(familyLabels), 1);
for i = 1:numel(familyLabels)
    idx = featureMeta.FamilyLabel == familyLabels(i);
    n(i) = nnz(idx);
    featureList(i) = strjoin(featureMeta.Name(idx), ";");
end
T = table(familyLabels, n, featureList, rationale, ...
    'VariableNames', {'FamilyLabel','n_features','feature_names','why_in_stable_dictionary'});
end

function T = local_demo_timeseries_table(dyad, featureNames)
T = table();
T.frame = (1:numel(dyad.time_s))';
T.time_s = dyad.time_s;
T.frame_valid = logical(dyad.frameMask);
T.preprocessing_badframe = logical(dyad.badframeMask);
for i = 1:numel(featureNames)
    T.(featureNames{i}) = dyad.X(:, i);
end
end

function T = local_pose_snapshot_table(poseAudit, nodeMetadata)
rows = struct([]);
for s = 1:numel(poseAudit.snapshotFrames)
    frame = poseAudit.snapshotFrames(s);
    for animal = 1:2
        for n = 1:height(nodeMetadata)
            row = struct();
            row.snapshot_index = s;
            row.frame = frame;
            row.time_s = poseAudit.time_s(s);
            row.animal = animal;
            row.node_index = nodeMetadata.node_index(n);
            row.bodypart_name = string(nodeMetadata.bodypart_name(n));
            row.matlab_field = string(nodeMetadata.matlab_field(n));
            row.x = poseAudit.tracks(frame, nodeMetadata.node_index(n), 1, animal);
            row.y = poseAudit.tracks(frame, nodeMetadata.node_index(n), 2, animal);
            row.preprocessing_badframe = poseAudit.badframes(frame);
            rows = [rows; row]; %#ok<AGROW>
        end
    end
end
T = struct2table(rows);
end

function files = local_plot_demo_geometry(posePath, figDir, cfg)
P = readtable(posePath, 'TextType', 'string');
snapshots = unique(P.snapshot_index, 'stable');
fig = figure('Visible','off', 'Color','w', 'Position',[80 80 1400 560]);
tl = tiledlayout(fig, 1, numel(snapshots), 'TileSpacing','compact', 'Padding','loose');
for s = 1:numel(snapshots)
    nexttile(tl);
    hold on;
    S = P(P.snapshot_index == snapshots(s), :);
    for animal = 1:2
        A = S(S.animal == animal, :);
        c = [0.00 0.45 0.70] .* (animal == 1) + [0.80 0.47 0.10] .* (animal == 2);
        local_plot_body_axis(A, c);
        scatter(A.x, A.y, 18, c, 'filled', 'MarkerFaceAlpha', 0.75);
    end
    axis equal;
    xlim([40 640]);
    ylim([120 290]);
    title(sprintf('t = %.2f s', S.time_s(1)), 'Interpreter','none');
    local_style_axes(gca, cfg);
end
xlabel(tl, 'x (pixels; feature scale from config)', 'Interpreter','none');
ylabel(tl, 'y (pixels)', 'Interpreter','none');
files = local_export_figure(fig, figDir, 'demo_dyad_feature_geometry', cfg);
close(fig);
end

function local_plot_body_axis(A, color)
ordered = ["tailTip","tailMid","tailBase","midBody","body","neck","nose"];
xy = nan(numel(ordered), 2);
for i = 1:numel(ordered)
    r = A(A.matlab_field == ordered(i), :);
    if ~isempty(r)
        xy(i,:) = [r.x(1), r.y(1)];
    end
end
plot(xy(:,1), xy(:,2), '-', 'Color', color, 'LineWidth', 2);
end

function files = local_plot_demo_timeseries(timeseriesPath, figDir, cfg)
T = readtable(timeseriesPath, 'TextType', 'string');
features = ["centroid_dist","mutual_facing","radial_speed_12", ...
    "heading_diff_deg","in_contact","asym_investigate"];
fig = figure('Visible','off', 'Color','w', 'Position',[80 80 1200 950]);
tl = tiledlayout(fig, numel(features), 1, 'TileSpacing','compact', 'Padding','compact');
for i = 1:numel(features)
    nexttile(tl);
    hold on;
    invalid = ~T.frame_valid;
    if any(invalid)
        ylPreview = [min(T.(features(i)), [], 'omitnan'), max(T.(features(i)), [], 'omitnan')];
        if all(isfinite(ylPreview)) && ylPreview(1) < ylPreview(2)
            for r = reshape(find(diff([false; invalid; false]) ~= 0), 2, [])
                patch(T.time_s([r(1) r(2)-1 r(2)-1 r(1)]), ...
                    [ylPreview(1) ylPreview(1) ylPreview(2) ylPreview(2)], ...
                    [0.86 0.86 0.86], 'EdgeColor','none', 'FaceAlpha',0.45);
            end
        end
    end
    plot(T.time_s, T.(features(i)), 'Color', [0.10 0.28 0.48], 'LineWidth', 1.4);
    ylabel(strrep(features(i), "_", " "), 'Interpreter','none');
    if i == numel(features)
        xlabel('time (s)', 'Interpreter','none');
    end
    local_style_axes(gca, cfg);
end
files = local_export_figure(fig, figDir, 'demo_dyad_feature_timeseries', cfg);
close(fig);
end

function files = local_plot_demo_family_summary(familySummaryPath, figDir, cfg)
T = readtable(familySummaryPath, 'TextType', 'string');
fig = figure('Visible','off', 'Color','w', 'Position',[80 80 1000 760]);
barh(T.n_features, 'FaceColor', [0.20 0.45 0.55], 'EdgeColor','none');
set(gca, 'YTick', 1:height(T), 'YTickLabel', T.FamilyLabel, ...
    'TickLabelInterpreter','none');
xlabel('Number of stable features', 'Interpreter','none');
title('Canonical dyadic feature families', 'Interpreter','none');
local_style_axes(gca, cfg);
files = local_export_figure(fig, figDir, 'demo_feature_family_summary', cfg);
close(fig);
end

function rows = local_add_manifest_rows(rows, files, sourceCsv, description)
for i = 1:numel(files)
    row = struct();
    row.figure_file = string(files(i));
    row.source_csv = string(sourceCsv);
    row.description = string(description);
    rows = [rows; row]; %#ok<AGROW>
end
end

function files = local_export_figure(fig, figDir, baseName, cfg)
files = strings(0,1);
if cfg.figures.export_png
    pngFile = fullfile(figDir, [baseName '.png']);
    exportgraphics(fig, pngFile, 'Resolution', cfg.figures.dpi);
    files(end+1,1) = string(pngFile);
end
if cfg.figures.export_pdf
    pdfFile = fullfile(figDir, [baseName '.pdf']);
    exportgraphics(fig, pdfFile, 'ContentType','vector');
    files(end+1,1) = string(pdfFile);
end
end

function local_style_axes(ax, cfg)
set(ax, 'FontName', char(cfg.figures.font_name), 'FontSize', cfg.figures.font_size, ...
    'LineWidth', 0.8, 'Box','off', 'TickDir','out');
grid(ax, 'on');
ax.GridColor = local_hex_to_rgb(cfg.figures.grid_color);
ax.GridAlpha = 0.35;
title(ax, ax.Title.String, 'FontSize', cfg.figures.title_font_size, ...
    'FontWeight','bold', 'Interpreter','none');
end

function rgb = local_hex_to_rgb(hex)
hex = char(strtrim(string(hex)));
rgb = sscanf(hex(2:end), '%2x%2x%2x', [1 3]) ./ 255;
end
