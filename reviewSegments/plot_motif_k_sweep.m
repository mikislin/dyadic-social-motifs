function Fig = plot_motif_k_sweep(Sweep)
%PLOT_MOTIF_K_SWEEP Plot segment-aware K-sweep metrics.

T = Sweep.table;
Fig = figure('Name', 'Motif K sweep diagnostics', 'Color', 'w');
tiledlayout(2,3, 'TileSpacing','compact', 'Padding','compact');

nexttile;
plot(T.K, T.meanARI, '-o');
xlabel('K'); ylabel('Bootstrap ARI'); title('Stability'); grid on;

nexttile;
plot(T.K, T.median_max_posterior, '-o');
xlabel('K'); ylabel('Median max posterior'); title('Confidence'); ylim([0 1]); grid on;

nexttile;
plot(T.K, T.median_segment_duration_s, '-o');
xlabel('K'); ylabel('Median segment duration (s)'); title('Bout duration'); grid on;

nexttile;
plot(T.K, T.short_duration_fraction, '-o');
xlabel('K'); ylabel('Short segment fraction'); title('Flicker after smoothing'); grid on;

nexttile;
plot(T.K, T.min_occupancy_fraction, '-o'); hold on;
plot(T.K, T.max_occupancy_fraction, '-o');
xlabel('K'); ylabel('Occupancy fraction'); title('Occupancy range');
legend({'min','max'}, 'Location','best'); grid on;

nexttile;
plot(T.K, T.transition_entropy_bits, '-o'); hold on;
plot(T.K, T.occupancy_entropy_bits, '-o');
xlabel('K'); ylabel('Entropy (bits)'); title('Entropy');
legend({'transition','occupancy'}, 'Location','best'); grid on;
end
