%DEMO_RUN08_CONDITION_BLIND_GRAPH_PRINCIPLES Compatibility entry point.
%
% The implementation lives in demo/demo04_run08_condition_blind_graph_principles.m.
% Required audit principles:
%   build_condition_blind_knn_graph
%   event_flags_used_for_demo_graph = false
%   labels_used_for_graph = "none"
%   demo_run08_event_coverage_by_scale.csv

repoRoot = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(repoRoot, 'demo', 'demo04_run08_condition_blind_graph_principles.m'));
