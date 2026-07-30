function Audit = analyze_run10_cross_session_reproducibility( ...
        Measurement, Registry, Validation, params)
%ANALYZE_RUN10_CROSS_SESSION_REPRODUCIBILITY Candidate profile reliability.

Data = prepare_run10_independent_feature_matrix(Measurement, Registry);
T = Data.table;
Z = Data.global_standardized;
labels = string(T.motif_candidate_id);
candidateIds = Validation.candidate_ids;
n = numel(candidateIds);

nodeCount = zeros(n,1);
sessionCount = zeros(n,1);
effectiveSessions = nan(n,1);
largestFraction = nan(n,1);
largestSessionIndex = nan(n,1);
splitCorr = nan(n,1);
largestRemovedCorr = nan(n,1);
medianSessionToRestCorr = nan(n,1);
for c = 1:n
    rows = find(labels == candidateIds(c));
    nodeCount(c) = numel(rows);
    sessions = unique(T.session_index(rows));
    sessionCount(c) = numel(sessions);
    if isempty(rows)
        continue
    end
    counts = zeros(numel(sessions),1);
    for j = 1:numel(sessions)
        counts(j) = nnz(T.session_index(rows) == sessions(j));
    end
    fractions = counts./sum(counts);
    effectiveSessions(c) = 1./sum(fractions.^2);
    [largestFraction(c),loc] = max(fractions);
    largestSessionIndex(c) = sessions(loc);
    fullProfile = mean(Z(rows,:),1);
    reducedRows = rows(T.session_index(rows) ~= largestSessionIndex(c));
    if ~isempty(reducedRows)
        largestRemovedCorr(c) = i_profile_corr( ...
            fullProfile,mean(Z(reducedRows,:),1));
    end
    if numel(sessions) >= 2
        stream = RandStream(char(params.random_stream_algorithm), ...
            'Seed',params.fold_seed+500+c);
        order = sessions(randperm(stream,numel(sessions)));
        aSessions = order(1:2:end);
        bSessions = order(2:2:end);
        aRows = rows(ismember(T.session_index(rows),aSessions));
        bRows = rows(ismember(T.session_index(rows),bSessions));
        if ~isempty(aRows) && ~isempty(bRows)
            splitCorr(c) = i_profile_corr( ...
                mean(Z(aRows,:),1),mean(Z(bRows,:),1));
        end
        sessionCorr = nan(numel(sessions),1);
        for j = 1:numel(sessions)
            one = rows(T.session_index(rows)==sessions(j));
            rest = rows(T.session_index(rows)~=sessions(j));
            if ~isempty(one) && ~isempty(rest)
                sessionCorr(j) = i_profile_corr( ...
                    mean(Z(one,:),1),mean(Z(rest,:),1));
            end
        end
        medianSessionToRestCorr(c) = ...
            median(sessionCorr,'omitnan');
    end
end

[tf,loc] = ismember(candidateIds, ...
    string(Validation.topology.motif_candidate_id));
assert(all(tf), ...
    'analyze_run10_cross_session_reproducibility:TopologyMismatch', ...
    'Every candidate requires topology provenance.');
eligible = logical(Validation.topology. ...
    eligible_for_behavioral_interpretation(loc));
sufficient = nodeCount >= params.active_minimum_nodes_primary & ...
    sessionCount >= params.minimum_sessions_primary & ...
    isfinite(splitCorr) & isfinite(largestRemovedCorr);
pass = eligible & sufficient & ...
    splitCorr >= params.cross_session_min_split_profile_correlation & ...
    largestRemovedCorr >= ...
    params.cross_session_min_largest_session_removed_correlation & ...
    largestFraction <= params.cross_session_max_largest_session_fraction;
status = repmat("descriptive_run09_residual",n,1);
status(eligible & ~sufficient) = "insufficient_evidence";
status(eligible & sufficient & ~pass) = "primary_gate_not_passed";
status(pass) = "primary_gate_passed";

Audit = table(candidateIds,nodeCount,sessionCount,effectiveSessions, ...
    largestSessionIndex,largestFraction,splitCorr,largestRemovedCorr, ...
    medianSessionToRestCorr,eligible,sufficient,pass,status, ...
    repmat("fixed_seed_candidate_session_split",n,1), ...
    repmat("candidate_centroids_in_global_robust_independent_feature_space", ...
    n,1), repmat(params.expected_membership_sha256,n,1), ...
    repmat("none",n,1), ...
    'VariableNames', {'motif_candidate_id','validation_node_count', ...
    'contributing_session_count','effective_session_count', ...
    'largest_contributing_session_index','largest_session_fraction', ...
    'split_session_profile_correlation', ...
    'largest_session_removed_profile_correlation', ...
    'median_session_to_rest_profile_correlation', ...
    'run09_graph_interpretation_eligible','sufficient_evidence', ...
    'cross_session_primary_gate_pass','cross_session_status', ...
    'session_split_rule','profile_definition', ...
    'frozen_membership_sha256','experimental_labels_used'});
end

function value = i_profile_corr(a,b)
ok = isfinite(a)&isfinite(b);
if nnz(ok)<3 || std(a(ok))<=eps || std(b(ok))<=eps
    value = NaN;
else
    C = corrcoef(a(ok),b(ok));
    value = C(1,2);
end
end
