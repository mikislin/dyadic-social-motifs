function [Audit, NullAudit] = ...
        analyze_run10_candidate_behavioral_coherence( ...
        Measurement, Registry, Validation, params)
%ANALYZE_RUN10_CANDIDATE_BEHAVIORAL_COHERENCE Blocked-null dispersion.

Data = prepare_run10_independent_feature_matrix(Measurement, Registry);
T = Data.table;
Z = Data.global_standardized;
labels = string(T.motif_candidate_id);
candidateIds = Validation.candidate_ids;
nCandidates = numel(candidateIds);
nPerm = params.active_permutation_count;

nodeCount = zeros(nCandidates,1);
sessionCount = zeros(nCandidates,1);
observed = nan(nCandidates,1);
nullValues = nan(nCandidates,nPerm);
for c = 1:nCandidates
    mask = labels == candidateIds(c);
    nodeCount(c) = nnz(mask);
    sessionCount(c) = numel(unique(T.session_index(mask)));
    observed(c) = i_dispersion(Z(mask,:));
end

stratum = findgroups(T.session_index, T.scale_index);
for b = 1:nPerm
    stream = RandStream(char(params.random_stream_algorithm), ...
        'Seed', params.permutation_seed + b - 1);
    permuted = labels;
    for g = 1:max(stratum)
        rows = find(stratum == g);
        if numel(rows) > 1
            permuted(rows) = labels(rows(randperm(stream, numel(rows))));
        end
    end
    for c = 1:nCandidates
        nullValues(c,b) = i_dispersion(Z(permuted == candidateIds(c),:));
    end
end

nullMean = mean(nullValues,2,'omitnan');
nullSd = std(nullValues,0,2,'omitnan');
nullLower = i_row_percentile(nullValues,2.5);
nullUpper = i_row_percentile(nullValues,97.5);
relativeImprovement = (nullMean - observed) ./ max(abs(nullMean), eps);
standardizedEffect = (nullMean - observed) ./ max(nullSd, eps);
pValue = nan(nCandidates,1);
for c = 1:nCandidates
    finite = isfinite(nullValues(c,:));
    if isfinite(observed(c)) && any(finite)
        pValue(c) = (1 + nnz(nullValues(c,finite) <= observed(c))) ./ ...
            (1 + nnz(finite));
    end
end
qValue = i_bh(pValue);

ciLower = nan(nCandidates,1);
ciUpper = nan(nCandidates,1);
for c = 1:nCandidates
    rows = find(labels == candidateIds(c));
    sessions = unique(T.session_index(rows));
    if numel(sessions) < 2
        continue
    end
    stream = RandStream(char(params.random_stream_algorithm), ...
        'Seed', params.bootstrap_seed + c);
    values = nan(params.active_bootstrap_count,1);
    for b = 1:params.active_bootstrap_count
        sampledSessions = sessions(randi(stream, numel(sessions), ...
            numel(sessions), 1));
        bootRows = zeros(0,1);
        for j = 1:numel(sampledSessions)
            bootRows = [bootRows; rows(T.session_index(rows) == ...
                sampledSessions(j))]; %#ok<AGROW>
        end
        values(b) = i_dispersion(Z(bootRows,:));
    end
    ciLower(c) = prctile(values(isfinite(values)),2.5);
    ciUpper(c) = prctile(values(isfinite(values)),97.5);
end

[tfTopology, topologyLoc] = ismember(candidateIds, ...
    string(Validation.topology.motif_candidate_id));
assert(all(tfTopology), ...
    'analyze_run10_candidate_behavioral_coherence:TopologyMismatch', ...
    'Every candidate requires topology provenance.');
eligible = logical(Validation.topology. ...
    eligible_for_behavioral_interpretation(topologyLoc));
sufficient = nodeCount >= params.active_minimum_nodes_primary & ...
    sessionCount >= params.minimum_sessions_primary & isfinite(observed);
pass = eligible & sufficient & ...
    relativeImprovement >= params.coherence_min_relative_improvement & ...
    qValue <= params.coherence_max_fdr_q;
status = repmat("descriptive_run09_residual",nCandidates,1);
status(eligible & ~sufficient) = "insufficient_evidence";
status(eligible & sufficient & ~pass) = "primary_gate_not_passed";
status(pass) = "primary_gate_passed";

Audit = table(candidateIds, nodeCount, sessionCount, observed, ciLower, ...
    ciUpper, nullMean, nullSd, nullLower, nullUpper, ...
    relativeImprovement, standardizedEffect, pValue, qValue, eligible, ...
    sufficient, pass, status, ...
    repmat(string(params.coherence_distance_metric),nCandidates,1), ...
    repmat(string(params.coherence_standardization_rule),nCandidates,1), ...
    repmat("candidate_labels_permuted_within_session_by_scale", ...
    nCandidates,1), repmat(nPerm,nCandidates,1), ...
    repmat(params.expected_membership_sha256,nCandidates,1), ...
    repmat("none",nCandidates,1), ...
    'VariableNames', {'motif_candidate_id','validation_node_count', ...
    'contributing_session_count','observed_within_dispersion', ...
    'observed_dispersion_ci_lower','observed_dispersion_ci_upper', ...
    'null_mean_dispersion','null_sd_dispersion', ...
    'null_dispersion_p025','null_dispersion_p975', ...
    'relative_dispersion_improvement','standardized_coherence_effect', ...
    'permutation_p_value','bh_fdr_q_value', ...
    'run09_graph_interpretation_eligible','sufficient_evidence', ...
    'coherence_primary_gate_pass','coherence_status', ...
    'distance_metric','standardization_rule','null_blocking_rule', ...
    'permutation_count','frozen_membership_sha256', ...
    'experimental_labels_used'});

candidateColumn = repelem(candidateIds,nPerm);
permutationIndex = repmat((1:nPerm)',nCandidates,1);
value = reshape(nullValues',[],1);
NullAudit = table(repmat("within_candidate_coherence", ...
    numel(value),1), candidateColumn, permutationIndex, ...
    params.permutation_seed + permutationIndex - 1, value, ...
    repmat("within_session_by_scale_label_permutation",numel(value),1), ...
    repmat("candidate_membership_target_only",numel(value),1), ...
    repmat("none",numel(value),1), ...
    'VariableNames', {'analysis_component','motif_candidate_id', ...
    'permutation_index','permutation_seed','null_metric_value', ...
    'blocking_rule','permuted_field','experimental_labels_used'});
end

function value = i_dispersion(X)
if isempty(X)
    value = NaN;
    return
end
centroid = mean(X,1,'omitnan');
value = mean(sum((X-centroid).^2,2),'omitnan');
end

function q = i_bh(p)
q = nan(size(p));
finite = find(isfinite(p));
if isempty(finite)
    return
end
[sorted, order] = sort(p(finite));
m = numel(sorted);
adjusted = sorted .* m ./ (1:m)';
for i = m-1:-1:1
    adjusted(i) = min(adjusted(i),adjusted(i+1));
end
adjusted = min(adjusted,1);
q(finite(order)) = adjusted;
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
