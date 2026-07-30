function FoldManifest = build_run10_blocked_validation_folds( ...
        Measurement, params)
%BUILD_RUN10_BLOCKED_VALIDATION_FOLDS Assign whole sessions to fixed folds.

T = Measurement(Measurement.eligible_for_automated_analysis, ...
    {'session_index','session_id'});
sessions = unique(T, 'rows', 'stable');
sessions = sortrows(sessions, {'session_index','session_id'});
n = height(sessions);
assert(n >= 2, ...
    'build_run10_blocked_validation_folds:TooFewSessions', ...
    'At least two sessions are required for grouped validation.');
nFolds = min(params.grouped_fold_count, n);
nRows = zeros(n, 1);
for i = 1:n
    nRows(i) = nnz(T.session_index == sessions.session_index(i));
end

stream = RandStream(char(params.random_stream_algorithm), ...
    'Seed', params.fold_seed);
tieKey = rand(stream, n, 1);
sortTable = table(-nRows, tieKey, sessions.session_index, ...
    'VariableNames', {'negative_row_count','tie_key','session_index'});
[~, order] = sortrows(sortTable, ...
    {'negative_row_count','tie_key','session_index'});
fold = zeros(n, 1);
foldLoad = zeros(nFolds, 1);
for i = 1:n
    [~, target] = min(foldLoad);
    row = order(i);
    fold(row) = target;
    foldLoad(target) = foldLoad(target) + nRows(row);
end

FoldManifest = sessions;
FoldManifest.heldout_fold = fold;
FoldManifest.n_validation_rows = nRows;
FoldManifest.fold_assignment_rule = repmat( ...
    "whole_session_greedy_row_balance_with_fixed_seed_tie_break", n, 1);
FoldManifest.fold_seed = repmat(params.fold_seed, n, 1);
FoldManifest.session_in_exactly_one_test_fold = true(n, 1);
FoldManifest.candidate_labels_used_for_fold_assignment = false(n, 1);
FoldManifest.scale_used_for_fold_assignment = false(n, 1);
FoldManifest.experimental_labels_used = repmat("none", n, 1);
FoldManifest = sortrows(FoldManifest, 'session_index');

assert(all(FoldManifest.heldout_fold >= 1 & ...
    FoldManifest.heldout_fold <= nFolds) && ...
    numel(unique(FoldManifest.session_index)) == height(FoldManifest) && ...
    ~any(FoldManifest.candidate_labels_used_for_fold_assignment) && ...
    all(FoldManifest.experimental_labels_used == "none"), ...
    'build_run10_blocked_validation_folds:LeakageContractFailure', ...
    'Grouped session fold contract failed.');
end
