function Fig = plot_motif_transition_stats(Transitions, varargin)
%PLOT_MOTIF_TRANSITION_STATS Plot transition matrix and motif entropy.

p = inputParser;
p.addParameter('Title', 'Motif transition statistics', @(x)ischar(x) || isstring(x));
p.parse(varargin{:});
P = p.Results;

C = Transitions.counts;
Pr = Transitions.prob;
K = size(C,1);

Fig = figure('Name', char(P.Title), 'Color', 'w');
tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact');

nexttile;
imagesc(Pr);
axis image;
colorbar;
xlabel('To motif'); ylabel('From motif');
title('Transition probability');
xticks(1:K); yticks(1:K);

nexttile;
imagesc(C);
axis image;
colorbar;
xlabel('To motif'); ylabel('From motif');
title('Transition counts');
xticks(1:K); yticks(1:K);

nexttile;
bar(1:K, Transitions.outEntropy(:));
xlabel('Motif'); ylabel('Outgoing entropy (bits)');
title(sprintf('Global entropy %.2f bits', Transitions.globalEntropy));
box off;
end
