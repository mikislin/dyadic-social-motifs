function Audit = analyze_run10_cross_scale_consistency( ...
        Measurement, Registry, Validation, params)
%ANALYZE_RUN10_CROSS_SCALE_CONSISTENCY Audit profile scale dependence.

Data = prepare_run10_independent_feature_matrix(Measurement, Registry);
T = Data.table;
Z = Data.global_standardized;
labels = string(T.motif_candidate_id);
candidateIds = Validation.candidate_ids;
n = numel(candidateIds);

representedScaleCount = zeros(n,1);
minimumScaleRows = nan(n,1);
medianPairCorr = nan(n,1);
minimumPairCorr = nan(n,1);
maximumPairCorr = nan(n,1);
for c = 1:n
    rows = find(labels==candidateIds(c));
    scales = unique(T.scale_index(rows));
    counts = zeros(numel(scales),1);
    profiles = nan(numel(scales),size(Z,2));
    for j = 1:numel(scales)
        scaleRows = rows(T.scale_index(rows)==scales(j));
        counts(j) = numel(scaleRows);
        profiles(j,:) = mean(Z(scaleRows,:),1);
    end
    keep = counts>=params.minimum_nodes_per_scale;
    representedScaleCount(c) = nnz(keep);
    if any(keep)
        minimumScaleRows(c) = min(counts(keep));
    end
    profiles = profiles(keep,:);
    pair = nan(0,1);
    for i = 1:size(profiles,1)-1
        for j = i+1:size(profiles,1)
            pair(end+1,1) = i_corr(profiles(i,:),profiles(j,:)); %#ok<AGROW>
        end
    end
    pair = pair(isfinite(pair));
    if ~isempty(pair)
        medianPairCorr(c) = median(pair);
        minimumPairCorr(c) = min(pair);
        maximumPairCorr(c) = max(pair);
    end
end

sufficient = representedScaleCount>=params.minimum_scales_for_audit & ...
    isfinite(medianPairCorr);
status = repmat("insufficient_scale_coverage",n,1);
status(sufficient & medianPairCorr>= ...
    params.cross_scale_consistency_correlation) = ...
    "cross_scale_consistent";
status(sufficient & medianPairCorr< ...
    params.cross_scale_consistency_correlation) = ...
    "scale_dependent_documented";
[tf,loc] = ismember(candidateIds,string( ...
    Validation.topology.motif_candidate_id));
assert(all(tf), ...
    'analyze_run10_cross_scale_consistency:TopologyMismatch', ...
    'Every candidate requires topology provenance.');
eligible = logical(Validation.topology. ...
    eligible_for_behavioral_interpretation(loc));

Audit = table(candidateIds,representedScaleCount,minimumScaleRows, ...
    medianPairCorr,minimumPairCorr,maximumPairCorr,sufficient,status, ...
    eligible,repmat("audit_only_not_validation_status_gate",n,1), ...
    repmat("scale_identity_used_only_to_form_postfreeze_profiles",n,1), ...
    repmat(params.expected_membership_sha256,n,1), ...
    repmat("none",n,1), ...
    'VariableNames', {'motif_candidate_id', ...
    'sufficiently_represented_scale_count','minimum_rows_per_retained_scale', ...
    'median_pairwise_scale_profile_correlation', ...
    'minimum_pairwise_scale_profile_correlation', ...
    'maximum_pairwise_scale_profile_correlation', ...
    'sufficient_scale_evidence','cross_scale_status', ...
    'run09_graph_interpretation_eligible','status_role', ...
    'scale_usage','frozen_membership_sha256', ...
    'experimental_labels_used'});
end

function value = i_corr(a,b)
ok = isfinite(a)&isfinite(b);
if nnz(ok)<3 || std(a(ok))<=eps || std(b(ok))<=eps
    value = NaN;
else
    C = corrcoef(a(ok),b(ok));
    value = C(1,2);
end
end
