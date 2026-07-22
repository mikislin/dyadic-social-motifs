function Audit = compute_graph_event_prevalence_fold_audit(Coverage)
%COMPUTE_GRAPH_EVENT_PREVALENCE_FOLD_AUDIT Compare enriched and baseline prevalence.
%
% The comparison is post-fit and does not alter graph nodes or edges. A
% Haldane-Anscombe half-count is used so all reported prevalence-fold values
% remain finite when an event is absent from one comparison set.

required = {'event_id','event_rule','scale_index','chunk_sec', ...
    'n_graph_nodes','n_graph_event_nodes','baseline_n_graph_nodes', ...
    'baseline_n_graph_event_nodes'};
assert(istable(Coverage) && all(ismember(required, Coverage.Properties.VariableNames)), ...
    'compute_graph_event_prevalence_fold_audit:BadCoverageTable', ...
    'Coverage table is missing required baseline/current count columns.');

Audit = table();
for i = 1:height(Coverage)
    Audit = [Audit; local_row(Coverage.event_id(i), Coverage.event_rule(i), ...
        "per_scale", double(Coverage.scale_index(i)), double(Coverage.chunk_sec(i)), ...
        double(Coverage.n_graph_event_nodes(i)), double(Coverage.n_graph_nodes(i)), ...
        double(Coverage.baseline_n_graph_event_nodes(i)), ...
        double(Coverage.baseline_n_graph_nodes(i)))]; %#ok<AGROW>
end

events = unique(string(Coverage.event_id), 'stable');
for event = events'
    idx = string(Coverage.event_id) == event;
    rule = string(Coverage.event_rule(find(idx, 1)));
    Audit = [Audit; local_row(event, rule, "all_scales", NaN, NaN, ...
        sum(double(Coverage.n_graph_event_nodes(idx)), 'omitnan'), ...
        sum(double(Coverage.n_graph_nodes(idx)), 'omitnan'), ...
        sum(double(Coverage.baseline_n_graph_event_nodes(idx)), 'omitnan'), ...
        sum(double(Coverage.baseline_n_graph_nodes(idx)), 'omitnan'))]; %#ok<AGROW>
end
end

function one = local_row(eventId, eventRule, scope, scaleIndex, chunkSec, ...
    currentEvents, currentNodes, baselineEvents, baselineNodes)
one = table();
one.event_id = string(eventId);
one.event_rule = string(eventRule);
one.aggregation_scope = string(scope);
one.scale_index = scaleIndex;
one.chunk_sec = chunkSec;
one.current_event_nodes = currentEvents;
one.current_total_nodes = currentNodes;
one.baseline_event_nodes = baselineEvents;
one.baseline_total_nodes = baselineNodes;
one.current_prevalence = currentEvents ./ max(currentNodes, 1);
one.baseline_prevalence = baselineEvents ./ max(baselineNodes, 1);
currentSmooth = (currentEvents + 0.5) ./ max(currentNodes + 1, 1);
baselineSmooth = (baselineEvents + 0.5) ./ max(baselineNodes + 1, 1);
one.prevalence_fold_enriched_over_baseline = currentSmooth ./ max(baselineSmooth, eps);
one.log2_prevalence_fold = log2(one.prevalence_fold_enriched_over_baseline);
one.half_count_smoothing_used = true;
one.baseline_zero_event_flag = baselineEvents == 0;
one.current_zero_event_flag = currentEvents == 0;
one.comparison_role = "postfit_sampling_prevalence_audit_not_population_prevalence";
one.labels_used_for_prevalence_audit = "none";
one.arena_used_for_prevalence_audit = false;
one.condition_used_for_prevalence_audit = false;
end
