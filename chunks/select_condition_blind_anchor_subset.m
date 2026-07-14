function selected = select_condition_blind_anchor_subset(candidates, maxAnchors, opts)
%SELECT_CONDITION_BLIND_ANCHOR_SUBSET Deterministically materialize anchors.
%
% The selector uses only session_index, raw_index, and anchor_time_s. It first
% gives each session time-even coverage when capacity allows, then fills the
% remaining budget by a global time-even pass. Provenance labels may be present
% in the input table, but they are not read by this function.

arguments
    candidates table
    maxAnchors (1,1) double {mustBePositive}
    opts.minAnchorsPerSession (1,1) double {mustBeGreaterThanOrEqual(opts.minAnchorsPerSession,0)} = 4
end

if isempty(candidates)
    selected = candidates;
    return
end

required = ["session_index", "raw_index", "anchor_time_s"];
missing = setdiff(required, string(candidates.Properties.VariableNames));
assert(isempty(missing), 'select_condition_blind_anchor_subset:MissingColumn', ...
    'Candidate table missing required columns: %s', strjoin(missing, ', '));

candidates = sortrows(candidates, {'raw_index', 'session_index', 'anchor_time_s'});
candidates.candidate_anchor_id = (1:height(candidates))';
maxAnchors = min(height(candidates), floor(maxAnchors));

selectedMask = false(height(candidates), 1);
sessions = unique(candidates.session_index, 'stable')';
if opts.minAnchorsPerSession > 0 && maxAnchors >= numel(sessions)
    for sess = sessions
        idx = find(candidates.session_index == sess);
        nWant = min([numel(idx), opts.minAnchorsPerSession, maxAnchors - nnz(selectedMask)]);
        if nWant <= 0
            break
        end
        pick = local_even_positions(idx, nWant);
        selectedMask(pick) = true;
    end
end

remaining = maxAnchors - nnz(selectedMask);
if remaining > 0
    idx = find(~selectedMask);
    pick = local_even_positions(idx, remaining);
    selectedMask(pick) = true;
end

selected = candidates(selectedMask, :);
selected = sortrows(selected, {'raw_index', 'session_index', 'anchor_time_s'});
selected.anchor_id = (1:height(selected))';
selected.anchor_selection_rank = selected.anchor_id;
selected.is_materialized_anchor = true(height(selected), 1);
selected.anchor_selection_rule = repmat( ...
    "deterministic_time_even_by_session_then_global_no_condition_labels", ...
    height(selected), 1);
end

function pick = local_even_positions(idx, nWant)
idx = idx(:);
nWant = min(nWant, numel(idx));
if nWant <= 0
    pick = zeros(0, 1);
elseif nWant >= numel(idx)
    pick = idx;
else
    pos = unique(round(linspace(1, numel(idx), nWant)));
    while numel(pos) < nWant
        missing = setdiff(1:numel(idx), pos, 'stable');
        pos(end + 1) = missing(1); %#ok<AGROW>
    end
    pick = idx(sort(pos(:)));
end
end
