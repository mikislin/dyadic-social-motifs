function Report = validate_scale_usefulness_scores(ScaleScore, varargin)
%VALIDATE_SCALE_USEFULNESS_SCORES Validate scale-usefulness scoring outputs.
%
% Example
%   Report = validate_scale_usefulness_scores(ScaleScore, 'makePlots', true);

p = inputParser;
p.addParameter('makePlots', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('verbose', true, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});
P = p.Results;

T = ScaleScore.scaleTable;
S = ScaleScore.selectedTable;
issues = strings(0,1);

reqVars = {'chunk_sec','initial_band','composite_micro','composite_motif','composite_context','composite_global'};
for i = 1:numel(reqVars)
    if ~ismember(reqVars{i}, T.Properties.VariableNames)
        issues(end+1,1) = "Missing required variable: " + reqVars{i}; %#ok<AGROW>
    end
end

if any(diff(T.chunk_sec) <= 0)
    issues(end+1,1) = "chunk_sec is not strictly increasing."; %#ok<AGROW>
end

if isempty(S)
    issues(end+1,1) = "No operational scales were selected."; %#ok<AGROW>
end

Report = struct();
Report.isValid = isempty(issues);
Report.issues = issues;
Report.scaleSummary = T(:, {'scale_index','chunk_sec','initial_band','pc1_explained','cum5_explained', ...
    'lag1_embedding_corr','label_run_frames','predict_short_r2','predict_long_r2', ...
    'composite_micro','composite_motif','composite_context','composite_global'});
Report.selectedSummary = S(:, intersect({'scale_index','chunk_sec','hierarchical_role','rank_within_role','composite_global'}, S.Properties.VariableNames, 'stable'));
Report.bandSummary = groupsummary(T, 'initial_band', {'mean','max'}, {'composite_micro','composite_motif','composite_context','chunk_sec'});

if P.verbose
    disp('=== Scale usefulness validation: scale summary ===');
    disp(Report.scaleSummary)
    disp('=== Scale usefulness validation: selected scales ===');
    disp(Report.selectedSummary)
    if isempty(issues)
        disp('No scale-usefulness validation issues detected.')
    else
        disp('Scale-usefulness validation issues:')
        disp(issues)
    end
end

if P.makePlots
    %plot_scale_usefulness_diagnostics(ScaleScore);
    plot_scale_usefulness_diagnostics(ScaleScore, ChunkSet)
end
end
