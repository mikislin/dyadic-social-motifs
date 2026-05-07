function Fig = plot_motif_segment_diagnostics(Segments, varargin)
%PLOT_MOTIF_SEGMENT_DIAGNOSTICS Plot anchor labels and segment duration summaries.
p = inputParser;
p.addRequired('Segments', @isstruct);
p.addParameter('ExampleSession', [], @(x)isempty(x) || isscalar(x));
p.parse(Segments, varargin{:});

A = Segments.anchorLabelTable;
T = Segments.table;
if isempty(p.Results.ExampleSession)
    sess = A.session_index(1);
else
    sess = p.Results.ExampleSession;
end
idx = A.session_index == sess;
if ismember('anchor_time_s', A.Properties.VariableNames) && any(isfinite(A.anchor_time_s(idx)))
    x = A.anchor_time_s(idx);
    xlab = 'Time (s)';
else
    x = A.anchor_frame(idx);
    xlab = 'Frame';
end

Fig = figure('Name','Motif segment diagnostics','Color','w');
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(x, A.raw_label(idx), '.', 'MarkerSize', 6); hold on;
plot(x, A.smoothed_label(idx), '-', 'LineWidth', 1);
ylabel('Motif'); title(sprintf('Session %d raw vs smoothed anchor labels', sess));
legend({'raw','smoothed'}, 'Location','best'); grid on;

nexttile;
if ~isempty(T)
    histogram(T.duration_s(isfinite(T.duration_s)), 40);
end
xlabel('Segment duration (s)'); ylabel('Count'); title('Segment duration distribution'); grid on;

nexttile;
if ~isempty(T)
    boxchart(categorical(T.motif_id), T.duration_s);
end
xlabel('Motif'); ylabel('Duration (s)'); title('Segment durations by motif'); grid on;
xlabel(xlab);
end
