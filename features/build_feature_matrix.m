function matrix = build_feature_matrix(sessionWindows)
%BUILD_FEATURE_MATRIX Combine session-wise window tables into one clustering matrix.
%
% Input
%   sessionWindows : struct array with field .table from summarize_windows
%
% Output
%   matrix.X       numeric matrix of summary features only
%   matrix.table   combined table including metadata columns
%   matrix.names   summary feature names

arguments
    sessionWindows (1,:) struct
end

T = table();
for i = 1:numel(sessionWindows)
    if isfield(sessionWindows(i), 'table') && ~isempty(sessionWindows(i).table)
        T = [T; sessionWindows(i).table]; %#ok<AGROW>
    end
end

matrix = struct();
matrix.table = T;
if isempty(T)
    matrix.X = [];
    matrix.names = {};
    return;
end

metaCols = ismember(T.Properties.VariableNames, {'session_idx','center_frame','center_time_s','window_sec','valid_frame_fraction'});
summaryNames = T.Properties.VariableNames(~metaCols);
matrix.X = table2array(T(:, summaryNames));
matrix.names = summaryNames;
end
