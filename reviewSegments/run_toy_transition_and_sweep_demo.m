%% Toy demo for motif transition statistics
% This demo does not require the full dyadic-social-motifs pipeline.
% It creates a toy Segments struct and computes transition statistics.

clear; clc;

motif_id = [1 1 2 2 3 2 1 4 4 2 3 3 1]';
session_index = ones(numel(motif_id),1);
start_time_s = (0:numel(motif_id)-1)' * 2;
stop_time_s = start_time_s + 1.5;
start_frame = round(start_time_s * 80);
stop_frame = round(stop_time_s * 80);
duration_s = stop_time_s - start_time_s;
n_anchors = round(duration_s / 0.1);
median_posterior = 0.9 + 0.1 * rand(numel(motif_id),1);

Segments = struct();
Segments.table = table(session_index, motif_id, start_frame, stop_frame, ...
    start_time_s, stop_time_s, duration_s, n_anchors, median_posterior);
Segments.summary = table(height(Segments.table), sum(n_anchors), median(duration_s), median(n_anchors), 0, ...
    'VariableNames', {'n_segments','n_anchors','median_segment_duration_s','median_segment_anchors','short_duration_fraction'});

Transitions = compute_motif_transition_stats(Segments, 'NumMotifs', 4);
disp(Transitions.table);
plot_motif_transition_stats(Transitions);
