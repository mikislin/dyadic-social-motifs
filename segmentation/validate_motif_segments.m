function SegReport = validate_motif_segments(Segments, Cluster, varargin)
%VALIDATE_MOTIF_SEGMENTS Summarize and sanity-check motif segments.
p = inputParser;
p.addRequired('Segments', @isstruct);
p.addRequired('Cluster', @isstruct);
p.addParameter('Verbose', true, @(x)islogical(x) || isnumeric(x));
p.parse(Segments, Cluster, varargin{:});

T = Segments.table;
A = Segments.anchorLabelTable;
K = Cluster.NumClusters;

issues = strings(0,1);
if isempty(T)
    issues(end+1) = "No segments were produced.";
end
if ~ismember('anchor_time_s', A.Properties.VariableNames)
    issues(end+1) = "anchor_time_s missing from anchorLabelTable.";
end
if any(A.smoothed_label < 0 | A.smoothed_label > K)
    issues(end+1) = "Smoothed labels outside expected range.";
end
if ~isempty(T) && any(T.stop_frame < T.start_frame)
    issues(end+1) = "Some segments have stop_frame < start_frame.";
end
if ~isempty(T) && any(T.motif_id < 1 | T.motif_id > K)
    issues(end+1) = "Some segments have invalid motif_id.";
end

SegReport = struct();
SegReport.summary = Segments.summary;
SegReport.motifSummary = Segments.motifSummary;
SegReport.sessionSummary = Segments.sessionSummary;
SegReport.issues = issues;
SegReport.ready = isempty(issues) || all(issues == "anchor_time_s missing from anchorLabelTable.");

if p.Results.Verbose
    disp('=== Segment validation summary ===');
    disp(SegReport.summary);
    disp('=== Motif summary ===');
    disp(SegReport.motifSummary);
    if isempty(issues)
        fprintf('No validation issues. Ready for raw-vs-smoothed comparison.\n');
    else
        fprintf('Validation issues:\n');
        disp(issues);
    end
end
end
