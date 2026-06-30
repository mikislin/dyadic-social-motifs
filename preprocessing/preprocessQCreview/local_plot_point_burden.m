function files = local_plot_point_burden(pointSummary, cfg, figDir)
if isempty(pointSummary) || height(pointSummary) == 0
    files = strings(0,1);
    return
end

metricVars = {'median_fraction_low_confidence_samples','median_fraction_jump_samples', ...
    'median_fraction_geometry_samples','median_fraction_interpolated_samples', ...
    'median_fraction_final_nan_samples','median_fraction_repaired_prediction_issue_samples', ...
    'median_fraction_unresolved_prediction_issue_samples'};
metricLabels = {'Low confidence','Jump','Geometry','Interpolated','Final NaN', ...
    'Repaired issue','Unresolved issue'};
X = zeros(height(pointSummary), numel(metricVars));
for j = 1:numel(metricVars)
    if ismember(metricVars{j}, pointSummary.Properties.VariableNames)
        X(:,j) = pointSummary.(metricVars{j}) .* 100;
    else
        X(:,j) = NaN;
    end
end

fig = figure('Visible','off', 'Color','w', 'Position',[80 80 1100 850]);
imagesc(X);
colormap(parula);
cb = colorbar;
cb.Label.String = 'Median % point samples';
set(gca, 'XTick', 1:numel(metricLabels), 'XTickLabel', metricLabels, ...
    'YTick', 1:height(pointSummary), 'YTickLabel', pointSummary.node_label);
xtickangle(35);
ylabel('SLEAP point');
title('Point-level preprocessing burden');
local_style_axes(gca, cfg);

files = local_export_figure(fig, figDir, 'preprocess_qc_point_level_bodypoints', cfg);
close(fig);
end
