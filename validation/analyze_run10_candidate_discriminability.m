function [Distinctness, Heldout, Confusion, NullAudit] = ...
        analyze_run10_candidate_discriminability( ...
        Measurement, Registry, FoldManifest, Validation, params)
%ANALYZE_RUN10_CANDIDATE_DISCRIMINABILITY Grouped nearest-centroid test.

Data = prepare_run10_independent_feature_matrix(Measurement, Registry);
T = Data.table;
X = Data.raw;
labels = string(T.motif_candidate_id);
candidateIds = Validation.eligible_candidate_ids;
primary = ismember(labels, candidateIds);
T = T(primary,:);
X = X(primary,:);
labels = labels(primary);
assert(~isempty(labels), ...
    'analyze_run10_candidate_discriminability:NoPrimaryRows', ...
    'No eligible-candidate validation rows are available.');

[present, loc] = ismember(T.session_index, FoldManifest.session_index);
assert(all(present), ...
    'analyze_run10_candidate_discriminability:MissingFold', ...
    'Every validation row requires a grouped-session fold.');
fold = FoldManifest.heldout_fold(loc);
nFolds = max(fold);
[predicted, foldStatus] = i_crossvalidated_predict( ...
    X, labels, fold, candidateIds);

nCandidates = numel(candidateIds);
testCount = zeros(nCandidates,1);
correctCount = zeros(nCandidates,1);
recall = nan(nCandidates,1);
ciLower = nan(nCandidates,1);
ciUpper = nan(nCandidates,1);
sessionCount = zeros(nCandidates,1);
for c = 1:nCandidates
    mask = labels == candidateIds(c) & predicted ~= "";
    testCount(c) = nnz(mask);
    correctCount(c) = nnz(predicted(mask) == candidateIds(c));
    sessionCount(c) = numel(unique(T.session_index(mask)));
    if testCount(c) > 0
        recall(c) = correctCount(c) ./ testCount(c);
        [ciLower(c),ciUpper(c)] = i_wilson( ...
            correctCount(c),testCount(c),1.96);
    end
end
macro = mean(recall,'omitnan');

trueId = repelem(candidateIds,nCandidates);
predId = repmat(candidateIds,nCandidates,1);
count = zeros(nCandidates^2,1);
rowFraction = nan(nCandidates^2,1);
for i = 1:numel(count)
    mask = labels == trueId(i) & predicted ~= "";
    count(i) = nnz(mask & predicted == predId(i));
    den = nnz(mask);
    if den > 0
        rowFraction(i) = count(i)./den;
    end
end
Confusion = table(trueId,predId,count,rowFraction, ...
    repmat("grouped_heldout_session_predictions",numel(count),1), ...
    repmat("none",numel(count),1), ...
    'VariableNames', {'true_motif_candidate_id', ...
    'predicted_motif_candidate_id','count','row_fraction', ...
    'prediction_scope','experimental_labels_used'});

nPerm = params.active_permutation_count;
nullRecall = nan(nCandidates,nPerm);
nullMacro = nan(nPerm,1);
stratum = findgroups(T.session_index,T.scale_index);
for b = 1:nPerm
    stream = RandStream(char(params.random_stream_algorithm), ...
        'Seed', params.permutation_seed + 100000 + b - 1);
    permuted = labels;
    for g = 1:max(stratum)
        rows = find(stratum == g);
        if numel(rows) > 1
            permuted(rows) = labels(rows(randperm(stream,numel(rows))));
        end
    end
    permPred = i_crossvalidated_predict(X,permuted,fold,candidateIds);
    for c = 1:nCandidates
        mask = permuted == candidateIds(c) & permPred ~= "";
        if any(mask)
            nullRecall(c,b) = mean(permPred(mask) == candidateIds(c));
        end
    end
    nullMacro(b) = mean(nullRecall(:,b),'omitnan');
end
nullRecallP95 = i_row_percentile(nullRecall, ...
    100*params.heldout_null_quantile);
nullRecallMean = mean(nullRecall,2,'omitnan');
nullMacroP95 = prctile(nullMacro(isfinite(nullMacro)), ...
    100*params.heldout_null_quantile);
nullMacroMean = mean(nullMacro,'omitnan');
pValue = nan(nCandidates,1);
for c = 1:nCandidates
    x = nullRecall(c,isfinite(nullRecall(c,:)));
    if isfinite(recall(c)) && ~isempty(x)
        pValue(c) = (1+nnz(x>=recall(c)))./(1+numel(x));
    end
end
macroP = (1+nnz(nullMacro>=macro))./(1+nnz(isfinite(nullMacro)));

sufficient = testCount >= params.active_minimum_nodes_primary & ...
    sessionCount >= params.minimum_sessions_primary & isfinite(recall);
overallPass = macro >= params.heldout_min_macro_balanced_accuracy && ...
    macro > nullMacroP95;
pass = sufficient & overallPass & ...
    recall >= params.heldout_min_candidate_recall & ...
    recall > nullRecallP95;
status = repmat("primary_gate_not_passed",nCandidates,1);
status(~sufficient) = "insufficient_evidence";
status(pass) = "primary_gate_passed";
Heldout = table(candidateIds,testCount,correctCount,recall,ciLower,ciUpper, ...
    sessionCount,nullRecallMean,nullRecallP95,pValue, ...
    repmat(macro,nCandidates,1),repmat(nullMacroMean,nCandidates,1), ...
    repmat(nullMacroP95,nCandidates,1),repmat(macroP,nCandidates,1), ...
    sufficient,pass,status, ...
    repmat(string(params.classifier_id),nCandidates,1), ...
    repmat(string(params.standardization_rule),nCandidates,1), ...
    repmat(string(params.classifier_missing_data_rule),nCandidates,1), ...
    repmat("whole_session_heldout",nCandidates,1), ...
    repmat("candidate_labels_permuted_within_session_by_scale", ...
    nCandidates,1), repmat(params.expected_membership_sha256, ...
    nCandidates,1), repmat("none",nCandidates,1), ...
    'VariableNames', {'motif_candidate_id','heldout_test_count', ...
    'heldout_correct_count','heldout_recall','recall_ci_lower', ...
    'recall_ci_upper','heldout_session_count','null_mean_recall', ...
    'null_recall_quantile','permutation_p_value', ...
    'macro_balanced_accuracy','null_mean_macro_balanced_accuracy', ...
    'null_macro_quantile','macro_permutation_p_value', ...
    'sufficient_evidence','heldout_primary_gate_pass', ...
    'heldout_status','classifier_id','standardization_rule', ...
    'missing_data_rule','fold_scope','null_blocking_rule', ...
    'frozen_membership_sha256','experimental_labels_used'});

Distinctness = table(string(params.classifier_id),nFolds, ...
    sum(testCount),macro,nullMacroMean,nullMacroP95,macroP,overallPass, ...
    strjoin(foldStatus,";"), ...
    string(params.standardization_rule), ...
    string(params.classifier_missing_data_rule), ...
    "macro_mean_per_candidate_recall", ...
    "candidate_id_prediction_target_not_biological_label", ...
    "none", ...
    'VariableNames', {'classifier_id','grouped_fold_count', ...
    'heldout_prediction_count','macro_balanced_accuracy', ...
    'null_mean_macro_balanced_accuracy','null_macro_quantile', ...
    'macro_permutation_p_value','overall_discriminability_gate_pass', ...
    'fold_status','standardization_rule','missing_data_rule', ...
    'primary_metric','target_role','experimental_labels_used'});

candidateColumn = [repelem(candidateIds,nPerm); ...
    repmat("__macro__",nPerm,1)];
permutationIndex = [repmat((1:nPerm)',nCandidates,1); (1:nPerm)'];
value = [reshape(nullRecall',[],1);nullMacro];
NullAudit = table(repmat("heldout_discriminability",numel(value),1), ...
    candidateColumn,permutationIndex, ...
    params.permutation_seed + 100000 + permutationIndex - 1, ...
    value,repmat("within_session_by_scale_label_permutation", ...
    numel(value),1),repmat("candidate_membership_target_only", ...
    numel(value),1),repmat("none",numel(value),1), ...
    'VariableNames', {'analysis_component','motif_candidate_id', ...
    'permutation_index','permutation_seed','null_metric_value', ...
    'blocking_rule','permuted_field','experimental_labels_used'});
end

function [predicted, status] = i_crossvalidated_predict( ...
        X,labels,fold,candidateIds)
n = size(X,1);
predicted = strings(n,1);
nFolds = max(fold);
status = strings(nFolds,1);
for f = 1:nFolds
    test = fold == f;
    train = ~test;
    if ~any(test) || ~any(train)
        status(f) = "empty_train_or_test";
        continue
    end
    trainMedian = median(X(train,:),1,'omitnan');
    trainMedian(~isfinite(trainMedian)) = 0;
    Xtrain = X(train,:);
    Xtest = X(test,:);
    for j = 1:size(X,2)
        Xtrain(~isfinite(Xtrain(:,j)),j) = trainMedian(j);
        Xtest(~isfinite(Xtest(:,j)),j) = trainMedian(j);
    end
    mu = mean(Xtrain,1);
    sigma = std(Xtrain,0,1);
    sigma(~isfinite(sigma)|sigma<=1e-12) = 1;
    Xtrain = (Xtrain-mu)./sigma;
    Xtest = (Xtest-mu)./sigma;
    centroids = nan(numel(candidateIds),size(X,2));
    available = false(numel(candidateIds),1);
    trainLabels = labels(train);
    for c = 1:numel(candidateIds)
        mask = trainLabels == candidateIds(c);
        if any(mask)
            centroids(c,:) = mean(Xtrain(mask,:),1);
            available(c) = true;
        end
    end
    if nnz(available) < 2
        status(f) = "fewer_than_two_training_classes";
        continue
    end
    availableIds = candidateIds(available);
    C = centroids(available,:);
    D = zeros(size(Xtest,1),size(C,1));
    for c = 1:size(C,1)
        D(:,c) = sum((Xtest-C(c,:)).^2,2);
    end
    [~,best] = min(D,[],2);
    predicted(test) = availableIds(best);
    status(f) = "success";
end
end

function [lower,upper] = i_wilson(k,n,z)
p = k./n;
den = 1+z^2/n;
center = (p+z^2/(2*n))./den;
half = z*sqrt(p*(1-p)/n+z^2/(4*n^2))./den;
lower = max(0,center-half);
upper = min(1,center+half);
end

function value = i_row_percentile(X,p)
value = nan(size(X,1),1);
for i = 1:size(X,1)
    x = X(i,isfinite(X(i,:)));
    if ~isempty(x)
        value(i) = prctile(x,p);
    end
end
end
