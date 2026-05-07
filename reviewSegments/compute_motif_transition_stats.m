function Transitions = compute_motif_transition_stats(Segments, varargin)
%COMPUTE_MOTIF_TRANSITION_STATS Segment-level motif transition statistics.
%
% Transitions = compute_motif_transition_stats(Segments, ...)
%
% Inputs
%   Segments : output of build_motif_segments. Requires Segments.table with
%              session_index and motif_id/cluster plus start/stop fields.
%
% Name-value options
%   NumMotifs             []      inferred from segment labels
%   RemoveSelfTransitions true    ignore consecutive same-motif transitions
%   Normalize             'row'   'row', 'joint', or 'none'
%   Verbose               true
%
% Output fields
%   Transitions.counts          K x K transition-count matrix
%   Transitions.prob           normalized transition matrix
%   Transitions.outEntropy     per-motif outgoing entropy, bits
%   Transitions.inEntropy      per-motif incoming entropy, bits
%   Transitions.globalEntropy  entropy of all nonzero transition events, bits
%   Transitions.table          edge table with counts/probabilities
%   Transitions.sessionTable   per-session transition summary
%   Transitions.params         options used

p = inputParser;
p.addParameter('NumMotifs', [], @(x) isempty(x) || (isscalar(x) && x >= 1));
p.addParameter('RemoveSelfTransitions', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Normalize', 'row', @(x) any(strcmpi(string(x), ["row","joint","none"])));
p.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
P = p.Results;

T = local_get_segment_table(Segments);
labelName = local_label_name(T);

if isempty(T)
    error('compute_motif_transition_stats:EmptySegments', 'Segments.table is empty.');
end

T = local_sort_segments(T);
labels = double(T.(labelName));
if isempty(P.NumMotifs)
    K = max(labels, [], 'omitnan');
else
    K = P.NumMotifs;
end
K = max(K, 1);

counts = zeros(K, K);
sessions = unique(T.session_index, 'stable');
sessionRows = struct([]);
r = 0;

for si = 1:numel(sessions)
    s = sessions(si);
    idx = find(T.session_index == s);
    z = double(T.(labelName)(idx));
    z = z(isfinite(z) & z >= 1 & z <= K);
    C = zeros(K,K);
    nTrans = 0;
    nSelf = 0;
    for i = 1:(numel(z)-1)
        a = z(i); b = z(i+1);
        if a == b
            nSelf = nSelf + 1;
            if logical(P.RemoveSelfTransitions)
                continue
            end
        end
        C(a,b) = C(a,b) + 1;
        nTrans = nTrans + 1;
    end
    counts = counts + C;
    r = r + 1;
    sessionRows(r).session_index = s; %#ok<AGROW>
    sessionRows(r).n_segments = numel(z);
    sessionRows(r).n_transitions = nTrans;
    sessionRows(r).n_self_transitions_removed = nSelf;
    sessionRows(r).transition_entropy_bits = local_entropy(C(:));
end

prob = local_normalize_counts(counts, char(P.Normalize));
outEntropy = nan(K,1);
inEntropy = nan(K,1);
for k = 1:K
    outEntropy(k) = local_entropy(prob(k,:));
    if strcmpi(P.Normalize, 'row')
        cin = counts(:,k);
        inEntropy(k) = local_entropy(cin ./ max(sum(cin), eps));
    else
        inEntropy(k) = local_entropy(prob(:,k));
    end
end

globalEntropy = local_entropy(counts(:) ./ max(sum(counts(:)), eps));
edgeTable = local_edge_table(counts, prob);
if isempty(sessionRows)
    sessionTable = table();
else
    sessionTable = struct2table(sessionRows);
end

Transitions = struct();
Transitions.counts = counts;
Transitions.prob = prob;
Transitions.outEntropy = outEntropy;
Transitions.inEntropy = inEntropy;
Transitions.globalEntropy = globalEntropy;
Transitions.table = edgeTable;
Transitions.sessionTable = sessionTable;
Transitions.params = P;

if logical(P.Verbose)
    fprintf('compute_motif_transition_stats | motifs = %d | transitions = %d | global entropy = %.3f bits\n', ...
        K, sum(counts(:)), globalEntropy);
end
end

function T = local_get_segment_table(Segments)
if istable(Segments)
    T = Segments;
elseif isstruct(Segments) && isfield(Segments, 'table') && istable(Segments.table)
    T = Segments.table;
else
    error('compute_motif_transition_stats:BadInput', 'Input must be Segments struct with Segments.table or a segment table.');
end
end

function name = local_label_name(T)
vars = string(T.Properties.VariableNames);
if any(vars == "motif_id")
    name = 'motif_id';
elseif any(vars == "cluster")
    name = 'cluster';
elseif any(vars == "cluster_id")
    name = 'cluster_id';
else
    error('compute_motif_transition_stats:NoLabel', 'Segment table needs motif_id, cluster, or cluster_id.');
end
end

function T = local_sort_segments(T)
vars = string(T.Properties.VariableNames);
if any(vars == "start_time_s")
    T = sortrows(T, {'session_index','start_time_s'});
elseif any(vars == "start_frame")
    T = sortrows(T, {'session_index','start_frame'});
else
    T = sortrows(T, {'session_index'});
end
end

function P = local_normalize_counts(C, mode)
switch lower(mode)
    case 'row'
        P = C ./ max(sum(C,2), eps);
    case 'joint'
        P = C ./ max(sum(C(:)), eps);
    otherwise
        P = C;
end
end

function H = local_entropy(p)
p = p(:);
p = p(isfinite(p) & p > 0);
if isempty(p)
    H = 0;
else
    p = p ./ sum(p);
    H = -sum(p .* log2(p));
end
end

function E = local_edge_table(counts, prob)
K = size(counts,1);
from_motif = [];
to_motif = [];
count = [];
probability = [];
for i = 1:K
    for j = 1:K
        if counts(i,j) > 0
            from_motif(end+1,1) = i; %#ok<AGROW>
            to_motif(end+1,1) = j; %#ok<AGROW>
            count(end+1,1) = counts(i,j); %#ok<AGROW>
            probability(end+1,1) = prob(i,j); %#ok<AGROW>
        end
    end
end
E = table(from_motif, to_motif, count, probability);
E = sortrows(E, {'from_motif','probability'}, {'ascend','descend'});
end
