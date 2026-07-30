function [Summary, Agreement] = analyze_run10_rater_agreement( ...
        Ratings, Validation, params)
%ANALYZE_RUN10_RATER_AGREEMENT Summarize frozen blinded ratings.

candidateIds = Validation.candidate_ids;
[tf,loc] = ismember(candidateIds,string( ...
    Validation.topology.motif_candidate_id));
assert(all(tf), ...
    'analyze_run10_rater_agreement:TopologyMismatch', ...
    'Every candidate requires topology provenance.');
eligible = logical(Validation.topology. ...
    eligible_for_behavioral_interpretation(loc));
n = numel(candidateIds);

if isempty(Ratings)
    reviewStatus = repmat("not_evaluated_run09_residual",n,1);
    reviewStatus(eligible) = "awaiting_ratings";
    Summary = table(candidateIds,eligible,zeros(n,1),zeros(n,1), ...
        nan(n,1),nan(n,1),nan(n,1),nan(n,1),false(n,1), ...
        reviewStatus,repmat("ratings_absent_not_negative_evidence",n,1), ...
        repmat(params.expected_membership_sha256,n,1), ...
        repmat("none",n,1), ...
        'VariableNames', {'motif_candidate_id', ...
        'run09_graph_interpretation_eligible','rater_count', ...
        'tile_rating_count','mean_candidate_coherence_score', ...
        'mean_tile_ambiguity_score','boundary_minus_core_ambiguity_score', ...
        'interactive_tile_fraction','blinded_coherence_gate_pass', ...
        'blinded_review_status','missing_rating_policy', ...
        'frozen_membership_sha256','experimental_labels_used'});
    Agreement = table("__global__","candidate_overall_coherence_score", ...
        0,0,NaN,NaN,NaN,false,"awaiting_ratings","none", ...
        'VariableNames', {'agreement_scope','rating_field','rater_count', ...
        'item_count','exact_agreement_fraction','within_one_agreement_fraction', ...
        'krippendorff_alpha_ordinal','agreement_gate_pass', ...
        'agreement_status','experimental_labels_used'});
    return
end

ratingCandidates = unique(string(Ratings.motif_candidate_id));
assert(all(ismember(ratingCandidates,Validation.eligible_candidate_ids)), ...
    'analyze_run10_rater_agreement:ResidualRating', ...
    'Only graph-eligible candidates belong in the primary rating packet.');
tile = string(Ratings.rating_level)=="tile";
candidateLevel = string(Ratings.rating_level)=="candidate";
coherence = i_numeric(Ratings.candidate_overall_coherence_score);
ambiguity = i_numeric(Ratings.ambiguity_score);

alpha = i_ordinal_alpha(Ratings(candidateLevel,:),coherence(candidateLevel));
globalRaters = numel(unique(string(Ratings.rater_id)));
globalItems = numel(unique(string(Ratings.review_id(candidateLevel))));
globalPass = globalRaters>=params.human_minimum_raters && ...
    isfinite(alpha) && alpha>=params.human_agreement_threshold;

raterCount = zeros(n,1);
tileCount = zeros(n,1);
meanCoherence = nan(n,1);
meanAmbiguity = nan(n,1);
boundaryMinusCore = nan(n,1);
interactiveFraction = nan(n,1);
coherencePass = false(n,1);
reviewStatus = repmat("not_evaluated_run09_residual",n,1);
agreementExact = nan(n,1);
agreementWithinOne = nan(n,1);
for c = 1:n
    rows = string(Ratings.motif_candidate_id)==candidateIds(c);
    if ~any(rows)
        if eligible(c)
            reviewStatus(c) = "insufficient_ratings";
        end
        continue
    end
    raterCount(c) = numel(unique(string(Ratings.rater_id(rows))));
    tileRows = rows & tile;
    candidateRows = rows & candidateLevel;
    tileCount(c) = nnz(tileRows);
    meanCoherence(c) = mean(coherence(candidateRows),'omitnan');
    meanAmbiguity(c) = mean(ambiguity(tileRows),'omitnan');
    boundary = tileRows & logical(Ratings.hidden_boundary_ambiguity);
    core = tileRows & ~logical(Ratings.hidden_boundary_ambiguity);
    if any(boundary) && any(core)
        boundaryMinusCore(c) = mean(ambiguity(boundary),'omitnan') - ...
            mean(ambiguity(core),'omitnan');
    end
    interactiveFraction(c) = mean( ...
        string(Ratings.interactive_impression(tileRows))=="interactive");
    scores = coherence(candidateRows);
    [agreementExact(c),agreementWithinOne(c)] = i_pair_agreement(scores);
    coherencePass(c) = raterCount(c)>=params.human_minimum_raters && ...
        isfinite(meanCoherence(c)) && ...
        meanCoherence(c)>=params.human_candidate_coherence_threshold;
    if raterCount(c)<params.human_minimum_raters
        reviewStatus(c) = "insufficient_ratings";
    elseif coherencePass(c) && globalPass
        reviewStatus(c) = "blinded_review_supported";
    elseif coherencePass(c)
        reviewStatus(c) = "coherent_but_agreement_below_rule";
    else
        reviewStatus(c) = "blinded_review_not_supported";
    end
end

Summary = table(candidateIds,eligible,raterCount,tileCount,meanCoherence, ...
    meanAmbiguity,boundaryMinusCore,interactiveFraction,coherencePass, ...
    reviewStatus,repmat("missing_values_remain_missing",n,1), ...
    repmat(params.expected_membership_sha256,n,1), ...
    repmat("none",n,1), ...
    'VariableNames', {'motif_candidate_id', ...
    'run09_graph_interpretation_eligible','rater_count', ...
    'tile_rating_count','mean_candidate_coherence_score', ...
    'mean_tile_ambiguity_score','boundary_minus_core_ambiguity_score', ...
    'interactive_tile_fraction','blinded_coherence_gate_pass', ...
    'blinded_review_status','missing_rating_policy', ...
    'frozen_membership_sha256','experimental_labels_used'});

candidateRows = table(candidateIds, ...
    repmat("candidate_overall_coherence_score",n,1),raterCount, ...
    ones(n,1),agreementExact,agreementWithinOne,repmat(alpha,n,1), ...
    repmat(globalPass,n,1),reviewStatus,repmat("none",n,1), ...
    'VariableNames', {'agreement_scope','rating_field','rater_count', ...
    'item_count','exact_agreement_fraction', ...
    'within_one_agreement_fraction','krippendorff_alpha_ordinal', ...
    'agreement_gate_pass','agreement_status','experimental_labels_used'});
globalRow = table("__global__","candidate_overall_coherence_score", ...
    globalRaters,globalItems,mean(agreementExact,'omitnan'), ...
    mean(agreementWithinOne,'omitnan'),alpha,globalPass, ...
    i_if(globalPass,"agreement_gate_passed","agreement_gate_not_passed"), ...
    "none",'VariableNames',candidateRows.Properties.VariableNames);
Agreement = [candidateRows;globalRow];
end

function x = i_numeric(xIn)
if isnumeric(xIn)
    x = double(xIn);
else
    x = str2double(string(xIn));
end
end

function [exact,withinOne] = i_pair_agreement(scores)
scores = scores(isfinite(scores));
if numel(scores)<2
    exact = NaN;
    withinOne = NaN;
    return
end
d = abs(scores(:)-scores(:)');
mask = triu(true(size(d)),1);
exact = mean(d(mask)==0);
withinOne = mean(d(mask)<=1);
end

function alpha = i_ordinal_alpha(T,scores)
review = string(T.review_id);
observed = nan(0,1);
for id = unique(review)'
    x = scores(review==id);
    x = x(isfinite(x));
    if numel(x)>=2
        d = (x(:)-x(:)').^2;
        observed = [observed;d(triu(true(size(d)),1))]; %#ok<AGROW>
    end
end
allScores = scores(isfinite(scores));
if isempty(observed)||numel(allScores)<2
    alpha = NaN;
    return
end
expected = (allScores(:)-allScores(:)').^2;
expected = expected(triu(true(size(expected)),1));
de = mean(expected);
if de<=eps
    alpha = NaN;
else
    alpha = 1-mean(observed)./de;
end
end

function value = i_if(condition,a,b)
if condition
    value = string(a);
else
    value = string(b);
end
end
