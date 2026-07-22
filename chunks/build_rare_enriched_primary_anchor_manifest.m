function outputs = build_rare_enriched_primary_anchor_manifest(repoRoot, params)
%BUILD_RARE_ENRICHED_PRIMARY_ANCHOR_MANIFEST Build stage-2 anchor bank.
%
% The locked primary run_06/run_07/run_08 outputs remain inputs. Candidate
% eligibility uses only frame masks, transformed feature availability, time,
% session identifiers, run_06 event patterns, and baseline graph rare-stratum
% seeds. Experimental and arena labels are joined only after selection as
% provenance copied from the locked primary manifest.

if nargin < 1 || strlength(string(repoRoot)) == 0
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
end
if nargin < 2 || isempty(params)
    params = load_multiscale_chunking_config();
end
assert(string(params.anchor_manifest_mode) == "rare_enriched", ...
    'build_rare_enriched_primary_anchor_manifest:BadMode', ...
    'This builder requires anchor_manifest_mode=rare_enriched.');

repoRoot = string(repoRoot);
baselineChunkRoot = resolve_repo_path(repoRoot, params.baseline_chunk_input_dir);
baselineGraphRoot = resolve_repo_path(repoRoot, params.baseline_graph_input_dir);
outRoot = resolve_repo_path(repoRoot, params.output_dir);
assert(string(outRoot) ~= string(baselineChunkRoot) && string(outRoot) ~= string(baselineGraphRoot), ...
    'build_rare_enriched_primary_anchor_manifest:BaselineOverwriteBlocked', ...
    'Rare-enriched output root must differ from locked baseline roots.');
local_ensure_dir(outRoot);

basePath = fullfile(baselineChunkRoot, 'primary_scale_specific_anchor_manifest.csv');
scalePath = fullfile(baselineChunkRoot, 'primary_operational_scales.csv');
transformPath = fullfile(baselineChunkRoot, 'chunk_feature_transform_audit.csv');
eventPath = fullfile(baselineChunkRoot, 'primary_chunk_event_summary_audit.csv');
definitionPath = fullfile(baselineGraphRoot, char(params.rare_strata_definition_file));
seedPath = fullfile(baselineGraphRoot, char(params.rare_strata_seed_manifest_file));
required = [string(basePath), string(scalePath), string(transformPath), ...
    string(eventPath), string(definitionPath), string(seedPath)];
for i = 1:numel(required)
    assert(isfile(required(i)), 'build_rare_enriched_primary_anchor_manifest:MissingInput', ...
        'Missing locked baseline input: %s', required(i));
end

Base = local_read_csv(basePath);
PrimaryScales = sortrows(local_read_csv(scalePath), 'chunk_sec');
Transform = local_read_csv(transformPath);
BaselineEvent = local_read_csv(eventPath);
Definition = local_read_csv(definitionPath);
Seed = local_read_csv(seedPath);
local_assert_label_free_inputs(Base, BaselineEvent, Definition, Seed);

[Session, Base] = local_session_table_and_base(Base, params);
stats = local_stats_from_transform_audit(Transform);
targetPerScale = local_target_per_scale(params);
strata = local_strata_order();
quotaFraction = local_quota_fractions(params);

nScale = height(PrimaryScales);
poolCells = cell(nScale, 1);
countRows = table();
eligibilityRows = table();

fprintf('Rare-enriched run_06 companion path | sessions=%d | scales=%d | target/scale=%d\n', ...
    height(Session), nScale, targetPerScale);
for r = 1:height(Session)
    dyad = local_load_dyad(repoRoot, string(Session.feature_output_file(r)));
    Seq = prepare_dyad_timeseries(dyad, 'stats', stats);
    Seq.rowFiniteFrac = mean(isfinite(Seq.X), 2);
    Seq.rowFiniteFrac(~isfinite(Seq.rowFiniteFrac)) = 0;
    Seq.validCumsum = [0; cumsum(double(Seq.validMask(:)))];
    Seq.finiteCumsum = [0; cumsum(double(Seq.rowFiniteFrac(:)))];
    assert(abs(double(Seq.fps) - double(params.fps)) <= 1e-10, ...
        'build_rare_enriched_primary_anchor_manifest:FpsMismatch', ...
        'Feature-session FPS does not match run_06 config.');

    for s = 1:nScale
        scaleIndex = double(PrimaryScales.scale_index(s));
        scaleSec = double(PrimaryScales.chunk_sec(s));
        A = find_condition_blind_chunk_anchors(Seq, scaleSec, ...
            'strideSec', params.stride_sec, ...
            'anchorMode', lower(string(params.anchor_mode)), ...
            'minValidFrac', params.min_chunk_valid_frac, ...
            'minFeatureFiniteFrac', params.min_chunk_feature_finite_frac, ...
            'requireAnchorFrameValid', logical(params.require_anchor_frame_valid));

        oneCount = table(scaleIndex, double(Session.session_index(r)), ...
            double(Session.raw_index(r)), string(Session.session_id(r)), height(A), ...
            'VariableNames', {'scale_index','session_index','raw_index','session_id','n_valid_candidates'});
        countRows = [countRows; oneCount]; %#ok<AGROW>
        if isempty(A)
            continue
        end

        C = local_candidate_rows(A, Seq, dyad, Session(r, :), PrimaryScales(s, :), ...
            Definition, Seed, Base, params);
        for j = 1:numel(strata)
            flagName = "flag_" + strata(j);
            one = table(scaleIndex, double(Session.session_index(r)), ...
                string(strata(j)), nnz(C.(flagName)), ...
                'VariableNames', {'scale_index','session_index','rare_stratum_id','n_eligible_candidates'});
            eligibilityRows = [eligibilityRows; one]; %#ok<AGROW>
        end

        essential = false(height(C), 1);
        for j = 1:(numel(strata) - 1)
            essential = essential | C.("flag_" + strata(j));
        end
        under = C.flag_undercovered_scale_session & ~essential;
        keepUnder = local_even_mask(C.anchor_time_s, under, ...
            params.rare_enriched_undercovered_candidates_per_cell);
        keep = essential | keepUnder;
        if any(keep)
            if isempty(poolCells{s})
                poolCells{s} = C(keep, :);
            else
                poolCells{s} = [poolCells{s}; C(keep, :)]; %#ok<AGROW>
            end
        end
    end
    if r == 1 || mod(r, 20) == 0 || r == height(Session)
        fprintf('  candidate sessions %d/%d\n', r, height(Session));
    end
end

Base.anchor_stage = repmat("base_time_even", height(Base), 1);
Base.anchor_source = repmat("locked_primary_scale_specific_anchor_manifest", height(Base), 1);
Base.anchor_manifest_mode = repmat("rare_enriched", height(Base), 1);
Base.source_primary_scale_specific_anchor_id = double(Base.primary_scale_specific_anchor_id);
Base.rare_stratum_id = repmat("none", height(Base), 1);
Base.rare_stratum_rule = repmat("not_rare_stratum_member", height(Base), 1);
Base.rare_stratum_score = zeros(height(Base), 1);
Base.rare_strata_membership_ids = repmat("none", height(Base), 1);
Base.rare_strata_membership_count = zeros(height(Base), 1);
Base.selection_phase = repmat("retained_base", height(Base), 1);
Base.quota_requested_stratum_id = repmat("none", height(Base), 1);
Base.final_assigned_rare_stratum_id = repmat("none", height(Base), 1);
Base.selection_composite_score = nan(height(Base), 1);
Base.fill_composite_score = nan(height(Base), 1);
Base.fill_reason = repmat("not_applicable_retained_base", height(Base), 1);
Base.duplicate_resolution_rule = repmat("unique_scale_session_anchor_frame_base_precedence", height(Base), 1);

RareManifest = Base(zeros(0, 1), :);
samplingTargets = table();
for s = 1:nScale
    C = poolCells{s};
    if isempty(C)
        C = local_empty_candidate_pool();
    end
    C.candidate_pool_row_id = (1:height(C))';
    baseScale = Base(double(Base.scale_index) == double(PrimaryScales.scale_index(s)), :);
    [selected, targetAudit] = local_select_rare_candidates(C, baseScale, targetPerScale, ...
        strata, quotaFraction);
    samplingTargets = [samplingTargets; targetAudit]; %#ok<AGROW>
    if isempty(selected)
        continue
    end
    rareOne = local_candidates_to_manifest(selected, Base);
    RareManifest = [RareManifest; rareOne]; %#ok<AGROW>
end

% Annotate retained baseline rows from the regenerated condition-blind pool.
nonemptyPool = ~cellfun(@isempty, poolCells);
if any(nonemptyPool)
    Pool = vertcat(poolCells{nonemptyPool});
else
    Pool = table();
end
Base = local_annotate_base_rare_membership(Base, Pool, strata);
Expanded = [Base; RareManifest];
Expanded = sortrows(Expanded, {'scale_index','raw_index','anchor_frame','anchor_stage'});
key = local_anchor_key(Expanded.scale_index, Expanded.raw_index, Expanded.anchor_frame);
assert(numel(unique(key)) == height(Expanded), ...
    'build_rare_enriched_primary_anchor_manifest:DuplicateAnchor', ...
    'Expanded bank contains duplicate scale/session/anchor-frame keys.');
Expanded.anchor_id = zeros(height(Expanded), 1);
Expanded.anchor_selection_rank = zeros(height(Expanded), 1);
for s = unique(double(Expanded.scale_index), 'stable')'
    idx = find(double(Expanded.scale_index) == s);
    Expanded.anchor_id(idx) = (1:numel(idx))';
    Expanded.anchor_selection_rank(idx) = (1:numel(idx))';
end
Expanded.expanded_scale_specific_anchor_id = Expanded.anchor_id;
Expanded.expanded_anchor_global_id = (1:height(Expanded))';
Expanded.labels_used_for_anchor_selection = repmat("none", height(Expanded), 1);
Expanded.arena_used_for_anchor_selection = false(height(Expanded), 1);
Expanded.condition_used_for_anchor_selection = false(height(Expanded), 1);
Expanded.expanded_anchor_rule = repmat( ...
    "retained_time_even_base_plus_condition_blind_event_graph_seed_and_coverage_strata", ...
    height(Expanded), 1);

[Expanded, WeightAudit] = local_inclusion_probabilities(Expanded, countRows, eligibilityRows);
EventSummary = summarize_primary_chunk_events(repoRoot, Session, Expanded, ...
    'contactThreshold', params.event_contact_threshold, ...
    'stateThreshold', params.event_state_threshold, ...
    'turnThresholdDeg', params.event_turn_threshold_deg);
local_assert_assigned_event_strata_match(Expanded, EventSummary);
BankSummary = local_bank_summary(Expanded, EventSummary, PrimaryScales, targetPerScale, ...
    baselineChunkRoot, baselineGraphRoot);
SamplingSummary = local_sampling_summary(Expanded, eligibilityRows, samplingTargets, ...
    PrimaryScales, strata);

writetable(Expanded, fullfile(outRoot, char(params.expanded_anchor_manifest_file)));
writetable(BankSummary, fullfile(outRoot, char(params.expanded_chunk_bank_summary_file)));
writetable(WeightAudit, fullfile(outRoot, char(params.expanded_weight_audit_file)));
writetable(SamplingSummary, fullfile(outRoot, char(params.rare_strata_sampling_summary_file)));
writetable(EventSummary, fullfile(outRoot, char(params.expanded_event_summary_file)));
writetable(countRows, fullfile(outRoot, 'expanded_candidate_scale_session_audit.csv'));
local_copy_baseline_audits(baselineChunkRoot, outRoot);
local_write_provenance_audit(outRoot, baselineChunkRoot, baselineGraphRoot, ...
    basePath, definitionPath, seedPath);

outputs = struct();
outputs.output_root = string(outRoot);
outputs.baseline_chunk_root = string(baselineChunkRoot);
outputs.baseline_graph_root = string(baselineGraphRoot);
outputs.expanded_anchor_manifest_path = string(fullfile(outRoot, char(params.expanded_anchor_manifest_file)));
outputs.expanded_bank_summary_path = string(fullfile(outRoot, char(params.expanded_chunk_bank_summary_file)));
outputs.expanded_weight_audit_path = string(fullfile(outRoot, char(params.expanded_weight_audit_file)));
outputs.rare_strata_sampling_summary_path = string(fullfile(outRoot, char(params.rare_strata_sampling_summary_file)));
outputs.expanded_event_summary_path = string(fullfile(outRoot, char(params.expanded_event_summary_file)));
outputs.n_expanded_anchors = height(Expanded);
outputs.n_base_anchors = nnz(Expanded.anchor_stage == "base_time_even");
outputs.n_rare_enriched_anchors = nnz(Expanded.anchor_stage == "rare_strata_enriched");
end

function local_assert_assigned_event_strata_match(T, E)
checks = ["contact_present","contact_dwell_fraction"; ...
    "contact_transition","contact_transition_count"; ...
    "large_heading_turn","heading_large_turn_count"];
for i = 1:size(checks, 1)
    idx = string(T.rare_stratum_id) == checks(i, 1);
    value = E.(checks(i, 2));
    assert(all(double(value(idx)) > 0), ...
        'build_rare_enriched_primary_anchor_manifest:EventStratumAuditMismatch', ...
        'Assigned %s anchors must be positive in saved event column %s.', ...
        checks(i, 1), checks(i, 2));
end
end

function [Session, Base] = local_session_table_and_base(Base, params)
[~, first] = unique(double(Base.feature_row_index), 'stable');
Session = sortrows(Base(first, :), {'raw_index','session_index'});
if string(params.run_mode) == "smoke"
    nSmoke = min(height(Session), floor(params.smoke_max_sessions));
    Session = Session(local_even_indices(height(Session), nSmoke), :);
end
Base = Base(ismember(string(Base.session_id), string(Session.session_id)), :);
Session.source_feature_row_index = double(Session.feature_row_index);
Session.feature_row_index = (1:height(Session))';
[tf, loc] = ismember(string(Base.session_id), string(Session.session_id));
assert(all(tf), 'build_rare_enriched_primary_anchor_manifest:SessionAlignment', ...
    'Baseline anchors could not be aligned to selected feature sessions.');
Base.feature_row_index = double(Session.feature_row_index(loc));
end

function idx = local_even_indices(nAvailable, nWant)
nWant = min(max(floor(nWant), 0), nAvailable);
if nWant == 0
    idx = zeros(0, 1);
elseif nWant == nAvailable
    idx = (1:nAvailable)';
else
    idx = unique(round(linspace(1, nAvailable, nWant))', 'stable');
    if numel(idx) < nWant
        fill = setdiff((1:nAvailable)', idx, 'stable');
        idx = sort([idx; fill(1:(nWant - numel(idx)))]);
    end
end
end

function target = local_target_per_scale(params)
target = floor(params.rare_enriched_target_anchors_per_scale);
if string(params.run_mode) == "smoke"
    target = min(target, floor(params.smoke_scale_specific_max_anchors_per_scale));
end
end

function C = local_candidate_rows(A, Seq, dyad, sessionRow, scaleRow, Definition, Seed, Base, params)
n = height(A);
L = max(1, round(double(scaleRow.chunk_sec) .* double(Seq.fps)));
switch lower(string(params.anchor_mode))
    case "past"
        left = L - 1; right = 0;
    otherwise
        left = floor((L - 1) / 2); right = ceil((L - 1) / 2);
end
st = double(A.anchor_frame) - left;
en = double(A.anchor_frame) + right;
M = local_candidate_event_metrics(dyad, st, en, params);

C = table();
C.scale_index = repmat(double(scaleRow.scale_index), n, 1);
C.primary_scale_rank = repmat(local_table_double(scaleRow, 'primary_scale_rank', 1, NaN), n, 1);
C.chunk_sec = repmat(double(scaleRow.chunk_sec), n, 1);
C.chunk_frames = repmat(L, n, 1);
C.feature_row_index = repmat(double(sessionRow.feature_row_index), n, 1);
C.session_index = repmat(double(sessionRow.session_index), n, 1);
C.raw_index = repmat(double(sessionRow.raw_index), n, 1);
C.anchor_frame = double(A.anchor_frame);
C.anchor_time_s = double(A.anchor_time_s);
C.start_frame = st;
C.stop_frame = en;
C.start_time_s = double(Seq.time(st));
C.stop_time_s = double(Seq.time(en));
C.valid_frac = double(A.min_scale_valid_frac);
C.feature_finite_frac = double(A.min_scale_feature_finite_frac);
C.anchor_frame_valid = logical(A.anchor_frame_valid);
C.chunk_is_valid = true(n, 1);
C.contact_dwell_fraction = M.contact_dwell_fraction;
C.contact_transition_count = M.contact_transition_count;
C.heading_large_turn_count = M.heading_large_turn_count;
C.radial_speed_mean_mm_s = M.radial_speed_mean_mm_s;
C.source_seed_primary_anchor_id = nan(n, 1);

strata = local_strata_order();
for j = 1:numel(strata)
    C.("flag_" + strata(j)) = false(n, 1);
    C.("score_" + strata(j)) = zeros(n, 1);
end
C.flag_contact_present = isfinite(M.contact_dwell_fraction) & M.contact_dwell_fraction > 0;
C.score_contact_present = max(M.contact_dwell_fraction, 0);
C.flag_contact_transition = isfinite(M.contact_transition_count) & M.contact_transition_count > 0;
C.score_contact_transition = max(M.contact_transition_count, 0);
C.flag_large_heading_turn = isfinite(M.heading_large_turn_count) & M.heading_large_turn_count > 0;
C.score_large_heading_turn = max(M.heading_large_turn_count, 0);

scaleIndex = double(scaleRow.scale_index);
radialThreshold = local_definition_threshold(Definition, "high_radial_speed", scaleIndex, Inf);
C.flag_high_radial_speed = isfinite(M.radial_speed_mean_mm_s) & M.radial_speed_mean_mm_s >= radialThreshold;
C.score_high_radial_speed = M.radial_speed_mean_mm_s ./ max(abs(radialThreshold), eps);

baseCellCount = nnz(double(Base.scale_index) == scaleIndex & ...
    double(Base.session_index) == double(sessionRow.session_index));
underThreshold = local_definition_threshold(Definition, "undercovered_scale_session", scaleIndex, -Inf);
C.flag_undercovered_scale_session(:) = baseCellCount <= underThreshold;
C.score_undercovered_scale_session(:) = (underThreshold + 1) ./ (baseCellCount + 1);

seedStrata = ["low_density_tail","graph_periphery","short_motif_instability"];
for j = 1:numel(seedStrata)
    idxSeed = double(Seed.scale_index) == scaleIndex & ...
        string(Seed.source_session_id) == string(sessionRow.session_id) & ...
        string(Seed.rare_stratum_id) == seedStrata(j);
    rows = find(idxSeed);
    for rr = rows(:)'
        d = abs(C.anchor_time_s - double(Seed.anchor_time_s(rr)));
        hit = d <= double(params.rare_enriched_seed_radius_sec);
        proximity = max(0, 1 - d ./ max(double(params.rare_enriched_seed_radius_sec), eps));
        candidateScore = proximity .* max(double(Seed.rare_stratum_score(rr)), eps);
        scoreName = "score_" + seedStrata(j);
        flagName = "flag_" + seedStrata(j);
        improve = hit & candidateScore >= C.(scoreName);
        C.(flagName) = C.(flagName) | hit;
        C.(scoreName)(improve) = candidateScore(improve);
        C.source_seed_primary_anchor_id(improve) = double(Seed.source_primary_scale_specific_anchor_id(rr));
    end
end

C.rare_composite_score = zeros(n, 1);
for j = 1:numel(strata)
    x = C.("score_" + strata(j));
    x(~isfinite(x) | x < 0) = 0;
    C.rare_composite_score = C.rare_composite_score + log1p(x) + double(C.("flag_" + strata(j)));
end
end

function M = local_candidate_event_metrics(dyad, st, en, params)
mask = local_frame_mask(dyad);
contactRaw = local_feature(dyad, 'in_contact');
heading = local_feature(dyad, 'heading_diff_deg');
radial = local_feature(dyad, 'radial_speed_12');

contactValid = mask & isfinite(contactRaw);
contact = contactRaw > params.event_contact_threshold & contactValid;
validCount = local_window_sum(double(mask), st, en);
contactCount = local_window_sum(double(contact), st, en);

% Match summarize_primary_chunk_events exactly: binary states are compressed
% over frameMask-valid observations before transitions are counted. This also
% prevents candidate eligibility and the saved post-selection audit from
% silently using different transition definitions.
contactTransition = local_compressed_transition_count(double(contact), mask, ...
    st, en, "binary", params.event_turn_threshold_deg);

headingValid = mask & isfinite(heading);
headingTurn = local_compressed_transition_count(heading, headingValid, ...
    st, en, "heading", params.event_turn_threshold_deg);

radialValid = mask & isfinite(radial);
radialSafe = radial;
radialSafe(~radialValid) = 0;
radialSum = local_window_sum(radialSafe, st, en);
radialCount = local_window_sum(double(radialValid), st, en);

M = struct();
M.contact_dwell_fraction = contactCount ./ max(validCount, 1);
M.contact_transition_count = contactTransition;
M.heading_large_turn_count = headingTurn;
M.radial_speed_mean_mm_s = radialSum ./ max(radialCount, 1);
M.radial_speed_mean_mm_s(radialCount == 0) = NaN;
end

function sumValue = local_window_sum(x, st, en)
cs = [0; cumsum(double(x(:)))];
sumValue = cs(en + 1) - cs(st);
end

function count = local_compressed_transition_count(x, valid, st, en, kind, threshold)
x = double(x(:));
valid = logical(valid(:));
validIdx = find(valid);
event = false(numel(valid), 1);
if numel(validIdx) >= 2
    switch string(kind)
        case "binary"
            changed = x(validIdx(2:end)) ~= x(validIdx(1:end-1));
        case "heading"
            h = unwrap(deg2rad(x(validIdx)));
            changed = abs(rad2deg(diff(h))) >= threshold;
        otherwise
            error('build_rare_enriched_primary_anchor_manifest:BadTransitionKind', ...
                'Unknown compressed transition kind: %s', kind);
    end
    event(validIdx(2:end)) = changed;
end
count = local_window_sum(double(event), st, en);
if isempty(validIdx)
    return
elseif numel(validIdx) == 1
    firstValid = repmat(validIdx, numel(st), 1);
    firstValid(double(st) > validIdx) = NaN;
else
    firstValid = interp1(double(validIdx), double(validIdx), double(st), 'next', NaN);
    firstValid(double(st) <= double(validIdx(1))) = double(validIdx(1));
end
hasFirst = isfinite(firstValid) & firstValid <= double(en);
firstValidSafe = firstValid;
firstValidSafe(~hasFirst) = 1;
count(hasFirst) = count(hasFirst) - double(event(firstValidSafe(hasFirst)));
end

function [Selected, TargetAudit] = local_select_rare_candidates(C, Base, targetPerScale, strata, quotaFraction)
Selected = C(zeros(0, 1), :);
TargetAudit = table();
nAdd = max(0, targetPerScale - height(Base));
target = floor(nAdd .* quotaFraction(:));
if ~isempty(target)
    target(end) = target(end) + nAdd - sum(target);
end
for j = 1:numel(strata)
    one = table();
    one.scale_index = double(Base.scale_index(1));
    one.rare_stratum_id = strata(j);
    one.target_rare_additions = target(j);
    one.n_selectable_pool_before_base_exclusion = 0;
    one.n_selectable_after_base_exclusion = 0;
    one.n_selectable_after_prior_quota_assignments = 0;
    one.n_excluded_as_locked_base = 0;
    one.n_depleted_by_prior_quota_assignments = 0;
    one.n_selected_quota_stage = 0;
    one.n_selected_fill_stage = 0;
    one.n_selected_exclusive_assignment = 0;
    one.n_selected_any_membership = 0;
    one.n_selected_with_multiple_memberships = 0;
    one.quota_shortfall_after_quota_stage = target(j);
    if target(j) > 0
        one.quota_shortfall_reason = "no_selectable_candidates";
    else
        one.quota_shortfall_reason = "none";
    end
    TargetAudit = [TargetAudit; one]; %#ok<AGROW>
end
if nAdd == 0 || isempty(C)
    return
end

baseKey = local_anchor_key(Base.scale_index, Base.raw_index, Base.anchor_frame);
poolKey = local_anchor_key(C.scale_index, C.raw_index, C.anchor_frame);
available = ~ismember(poolKey, baseKey);
chosen = false(height(C), 1);
assigned = strings(height(C), 1);
selectionPhase = repmat("not_selected", height(C), 1);
quotaRequested = repmat("none", height(C), 1);
fillReason = repmat("not_fill", height(C), 1);
fillCompositeScore = nan(height(C), 1);
for j = 1:numel(strata)
    poolBeforeBase = C.("flag_" + strata(j));
    afterBase = available & poolBeforeBase;
    eligible = available & ~chosen & C.("flag_" + strata(j));
    pick = local_balanced_pick(C, eligible, target(j), C.("score_" + strata(j)));
    chosen(pick) = true;
    assigned(pick) = strata(j);
    selectionPhase(pick) = "quota";
    quotaRequested(pick) = strata(j);

    TargetAudit.n_selectable_pool_before_base_exclusion(j) = nnz(poolBeforeBase);
    TargetAudit.n_selectable_after_base_exclusion(j) = nnz(afterBase);
    TargetAudit.n_selectable_after_prior_quota_assignments(j) = nnz(eligible);
    TargetAudit.n_excluded_as_locked_base(j) = nnz(poolBeforeBase & ~available);
    TargetAudit.n_depleted_by_prior_quota_assignments(j) = nnz(afterBase & chosen) - numel(pick);
    TargetAudit.n_selected_quota_stage(j) = numel(pick);
    TargetAudit.quota_shortfall_after_quota_stage(j) = target(j) - numel(pick);
    if TargetAudit.quota_shortfall_after_quota_stage(j) == 0
        TargetAudit.quota_shortfall_reason(j) = "none";
    elseif nnz(poolBeforeBase) < target(j)
        TargetAudit.quota_shortfall_reason(j) = "insufficient_selectable_pool";
    elseif nnz(afterBase) < target(j)
        TargetAudit.quota_shortfall_reason(j) = "locked_base_exclusion_depletion";
    elseif nnz(eligible) < target(j)
        TargetAudit.quota_shortfall_reason(j) = "priority_overlap_depletion";
    else
        TargetAudit.quota_shortfall_reason(j) = "unexpected_selection_shortfall";
    end
end
remaining = nAdd - nnz(chosen);
fillChosen = false(height(C), 1);
if remaining > 0
    pick = local_balanced_pick(C, available & ~chosen, remaining, C.rare_composite_score);
    chosen(pick) = true;
    fillChosen(pick) = true;
    selectionPhase(pick) = "fill";
    fillReason(pick) = "quota_shortfall_composite_score_fill";
    fillCompositeScore(pick) = C.rare_composite_score(pick);
    for r = pick(:)'
        assigned(r) = local_primary_stratum(C(r, :), strata);
    end
end
assert(nnz(chosen) == nAdd, ...
    'build_rare_enriched_primary_anchor_manifest:InsufficientRareCandidates', ...
    'Only %d of %d requested rare additions were available at scale_index=%g.', ...
    nnz(chosen), nAdd, double(Base.scale_index(1)));
membershipCount = zeros(height(C), 1);
for j = 1:numel(strata)
    membershipCount = membershipCount + double(C.("flag_" + strata(j)));
    TargetAudit.n_selected_fill_stage(j) = nnz(fillChosen & assigned == strata(j));
    TargetAudit.n_selected_exclusive_assignment(j) = nnz(chosen & assigned == strata(j));
    TargetAudit.n_selected_any_membership(j) = nnz(chosen & C.("flag_" + strata(j)));
end
for j = 1:numel(strata)
    TargetAudit.n_selected_with_multiple_memberships(j) = ...
        nnz(chosen & C.("flag_" + strata(j)) & membershipCount > 1);
end
assert(sum(TargetAudit.n_selected_quota_stage) + sum(TargetAudit.n_selected_fill_stage) == nAdd, ...
    'build_rare_enriched_primary_anchor_manifest:SelectionPhaseCountMismatch', ...
    'Quota-stage and fill-stage counts must account for every rare addition.');
Selected = C(chosen, :);
Selected.rare_stratum_id = assigned(chosen);
Selected.rare_stratum_rule = strings(height(Selected), 1);
Selected.rare_stratum_score = zeros(height(Selected), 1);
Selected.rare_strata_membership_ids = strings(height(Selected), 1);
for i = 1:height(Selected)
    Selected.rare_stratum_rule(i) = local_stratum_rule(Selected.rare_stratum_id(i));
    Selected.rare_stratum_score(i) = Selected.("score_" + Selected.rare_stratum_id(i))(i);
    hit = strings(0, 1);
    for j = 1:numel(strata)
        if Selected.("flag_" + strata(j))(i)
            hit(end + 1, 1) = strata(j); %#ok<AGROW>
        end
    end
    Selected.rare_strata_membership_ids(i) = strjoin(hit, ";");
end
Selected.rare_strata_membership_count = membershipCount(chosen);
Selected.selection_phase = selectionPhase(chosen);
Selected.quota_requested_stratum_id = quotaRequested(chosen);
Selected.final_assigned_rare_stratum_id = assigned(chosen);
Selected.selection_composite_score = C.rare_composite_score(chosen);
Selected.fill_composite_score = fillCompositeScore(chosen);
Selected.fill_reason = fillReason(chosen);
assert(all(Selected.rare_strata_membership_count >= 1), ...
    'build_rare_enriched_primary_anchor_manifest:MissingSelectedMembership', ...
    'Every rare addition must retain at least one condition-blind rare-stratum membership.');
end

function pick = local_balanced_pick(C, mask, nWant, score)
idx = find(mask);
nWant = min(floor(nWant), numel(idx));
if nWant <= 0
    pick = zeros(0, 1);
    return
end
T = table(idx, double(C.session_index(idx)), double(C.raw_index(idx)), ...
    double(C.anchor_time_s(idx)), double(score(idx)), ...
    'VariableNames', {'idx','session_index','raw_index','anchor_time_s','score'});
T.score(~isfinite(T.score)) = -Inf;
T = sortrows(T, {'session_index','score','anchor_time_s'}, {'ascend','descend','ascend'});
T.within_session_rank = zeros(height(T), 1);
for sess = unique(T.session_index, 'stable')'
    rows = find(T.session_index == sess);
    T.within_session_rank(rows) = (1:numel(rows))';
end
T = sortrows(T, {'within_session_rank','score','raw_index','anchor_time_s'}, ...
    {'ascend','descend','ascend','ascend'});
pick = T.idx(1:nWant);
end

function T = local_candidates_to_manifest(C, Base)
templateKey = local_cell_key(Base.scale_index, Base.session_index);
[u, first] = unique(templateKey, 'stable');
[tf, loc] = ismember(local_cell_key(C.scale_index, C.session_index), u);
assert(all(tf), 'build_rare_enriched_primary_anchor_manifest:MissingTemplate', ...
    'Rare candidates lack matching baseline scale-session provenance templates.');
T = Base(first(loc), :);
copyNames = ["scale_index","primary_scale_rank","chunk_sec","chunk_frames", ...
    "feature_row_index","session_index","raw_index","anchor_frame","anchor_time_s", ...
    "start_frame","stop_frame","start_time_s","stop_time_s","valid_frac", ...
    "feature_finite_frac","anchor_frame_valid","chunk_is_valid"];
for name = copyNames
    if ismember(name, C.Properties.VariableNames) && ismember(name, T.Properties.VariableNames)
        value = C.(name);
        if ~(isnumeric(value) && all(isnan(value)))
            T.(name) = value;
        end
    end
end
T.candidate_anchor_id = double(C.candidate_pool_row_id);
T.anchor_id = nan(height(T), 1);
T.anchor_selection_rank = nan(height(T), 1);
T.is_materialized_anchor = true(height(T), 1);
T.anchor_selection_rule = repmat("condition_blind_rare_strata_quota_then_session_balanced_fill", height(T), 1);
T.primary_scale_specific_anchor_id = nan(height(T), 1);
T.primary_anchor_global_id = nan(height(T), 1);
T.anchor_selection_scope = repmat("rare_enriched_primary_bank", height(T), 1);
T.provenance_role = repmat("provenance_joined_after_condition_blind_rare_selection", height(T), 1);
T.anchor_stage = repmat("rare_strata_enriched", height(T), 1);
T.anchor_source = repmat("condition_blind_event_and_graph_seed_candidate_pool", height(T), 1);
T.anchor_manifest_mode = repmat("rare_enriched", height(T), 1);
T.source_primary_scale_specific_anchor_id = double(C.source_seed_primary_anchor_id);
T.rare_stratum_id = string(C.rare_stratum_id);
T.rare_stratum_rule = string(C.rare_stratum_rule);
T.rare_stratum_score = double(C.rare_stratum_score);
T.rare_strata_membership_ids = string(C.rare_strata_membership_ids);
T.rare_strata_membership_count = double(C.rare_strata_membership_count);
T.selection_phase = string(C.selection_phase);
T.quota_requested_stratum_id = string(C.quota_requested_stratum_id);
T.final_assigned_rare_stratum_id = string(C.final_assigned_rare_stratum_id);
T.selection_composite_score = double(C.selection_composite_score);
T.fill_composite_score = double(C.fill_composite_score);
T.fill_reason = string(C.fill_reason);
T.duplicate_resolution_rule = repmat("unique_scale_session_anchor_frame_base_precedence", height(T), 1);
end

function Base = local_annotate_base_rare_membership(Base, Pool, strata)
if isempty(Pool)
    return
end
[tf, loc] = ismember(local_anchor_key(Base.scale_index, Base.raw_index, Base.anchor_frame), ...
    local_anchor_key(Pool.scale_index, Pool.raw_index, Pool.anchor_frame));
for i = find(tf(:))'
    row = Pool(loc(i), :);
    id = local_primary_stratum(row, strata);
    Base.rare_stratum_id(i) = id;
    Base.rare_stratum_rule(i) = local_stratum_rule(id);
    Base.rare_stratum_score(i) = row.("score_" + id);
    hit = strings(0, 1);
    for j = 1:numel(strata)
        if row.("flag_" + strata(j))
            hit(end + 1, 1) = strata(j); %#ok<AGROW>
        end
    end
    Base.rare_strata_membership_ids(i) = strjoin(hit, ";");
    Base.rare_strata_membership_count(i) = numel(hit);
    Base.final_assigned_rare_stratum_id(i) = id;
end
end

function [T, Audit] = local_inclusion_probabilities(T, Count, Eligibility)
n = height(T);
baseP = nan(n, 1);
rareP = zeros(n, 1);
cellKey = local_cell_key(T.scale_index, T.session_index);
countKey = local_cell_key(Count.scale_index, Count.session_index);
for key = unique(cellKey, 'stable')'
    idx = cellKey == key;
    row = find(countKey == key, 1);
    assert(~isempty(row) && Count.n_valid_candidates(row) >= 1, ...
        'build_rare_enriched_primary_anchor_manifest:MissingCandidateCount', ...
        'Missing candidate count for expanded anchor cell %s.', key);
    nBase = nnz(idx & T.anchor_stage == "base_time_even");
    baseP(idx) = min(1, nBase ./ double(Count.n_valid_candidates(row)));
end

rareRows = find(T.anchor_stage == "rare_strata_enriched");
for rr = rareRows(:)'
    same = double(T.scale_index) == double(T.scale_index(rr)) & ...
        double(T.session_index) == double(T.session_index(rr)) & ...
        string(T.rare_stratum_id) == string(T.rare_stratum_id(rr)) & ...
        T.anchor_stage == "rare_strata_enriched";
    erow = find(double(Eligibility.scale_index) == double(T.scale_index(rr)) & ...
        double(Eligibility.session_index) == double(T.session_index(rr)) & ...
        string(Eligibility.rare_stratum_id) == string(T.rare_stratum_id(rr)), 1);
    assert(~isempty(erow) && Eligibility.n_eligible_candidates(erow) >= nnz(same), ...
        'build_rare_enriched_primary_anchor_manifest:BadRareProbabilityDenominator', ...
        'Rare-stratum eligibility denominator is missing or smaller than selected count.');
    rareP(rr) = min(1, nnz(same) ./ double(Eligibility.n_eligible_candidates(erow)));
end

finalP = 1 - (1 - baseP) .* (1 - rareP);
assert(all(isfinite(finalP) & finalP > 0 & finalP <= 1), ...
    'build_rare_enriched_primary_anchor_manifest:InvalidInclusionProbability', ...
    'Final inclusion probabilities must be finite and in (0,1].');
T.base_inclusion_probability = baseP;
T.rare_inclusion_probability = rareP;
T.final_inclusion_probability = finalP;
T.pca_training_weight = ones(n, 1);
T.graph_training_weight = ones(n, 1);
T.audit_inverse_probability_weight = 1 ./ finalP;
T.audit_weight_normalized_within_scale = nan(n, 1);
for s = unique(double(T.scale_index), 'stable')'
    idx = double(T.scale_index) == s;
    T.audit_weight_normalized_within_scale(idx) = T.audit_inverse_probability_weight(idx) ./ ...
        mean(T.audit_inverse_probability_weight(idx), 'omitnan');
end
T.inclusion_probability_rule = repmat( ...
    "scale_session_base_rate_union_assigned_stratum_scale_session_rate", n, 1);
T.weights_used_for_pca = false(n, 1);
T.weights_used_for_graph_edges = false(n, 1);
T.audit_weight_interpretation = repmat( ...
    "assigned_stratum_enrichment_sensitivity_not_exact_overlap_adjusted_probability", n, 1);

keep = ["expanded_anchor_global_id","scale_index","chunk_sec","session_index", ...
    "raw_index","anchor_frame","anchor_stage","selection_phase", ...
    "quota_requested_stratum_id","rare_stratum_id","final_assigned_rare_stratum_id", ...
    "rare_strata_membership_ids","rare_strata_membership_count", ...
    "selection_composite_score","fill_composite_score","fill_reason", ...
    "base_inclusion_probability","rare_inclusion_probability", ...
    "final_inclusion_probability","pca_training_weight","graph_training_weight", ...
    "audit_inverse_probability_weight","audit_weight_normalized_within_scale", ...
    "inclusion_probability_rule","audit_weight_interpretation", ...
    "labels_used_for_anchor_selection", ...
    "arena_used_for_anchor_selection","condition_used_for_anchor_selection"];
Audit = T(:, cellstr(keep));
end

function Summary = local_bank_summary(T, E, PrimaryScales, target, baselineChunkRoot, baselineGraphRoot)
Summary = table();
for s = 1:height(PrimaryScales)
    idx = double(T.scale_index) == double(PrimaryScales.scale_index(s));
    one = table();
    one.scale_index = double(PrimaryScales.scale_index(s));
    one.chunk_sec = double(PrimaryScales.chunk_sec(s));
    one.hierarchical_role = local_table_string(PrimaryScales, 'hierarchical_role', s, "");
    one.target_anchors_per_scale = target;
    one.n_base_anchors = nnz(idx & T.anchor_stage == "base_time_even");
    one.n_rare_enriched_anchors = nnz(idx & T.anchor_stage == "rare_strata_enriched");
    one.n_expanded_anchors = nnz(idx);
    one.n_sessions = numel(unique(T.session_index(idx)));
    one.n_contact_present = nnz(idx & E.contact_dwell_fraction > 0);
    one.n_contact_transition = nnz(idx & E.contact_transition_count > 0);
    one.n_large_heading_turn = nnz(idx & E.heading_large_turn_count > 0);
    one.mean_final_inclusion_probability = mean(T.final_inclusion_probability(idx), 'omitnan');
    one.effective_sample_size_inverse_probability = local_effective_sample_size(T.audit_inverse_probability_weight(idx));
    one.baseline_run06_root = string(baselineChunkRoot);
    one.baseline_run08_root = string(baselineGraphRoot);
    one.labels_used_for_expanded_anchor_bank = "none";
    one.arena_used_for_expanded_anchor_bank = false;
    one.condition_used_for_expanded_anchor_bank = false;
    Summary = [Summary; one]; %#ok<AGROW>
end
end

function Summary = local_sampling_summary(T, Eligibility, Targets, PrimaryScales, strata)
Summary = table();
for s = 1:height(PrimaryScales)
    scaleIndex = double(PrimaryScales.scale_index(s));
    for j = 1:numel(strata)
        idxE = double(Eligibility.scale_index) == scaleIndex & string(Eligibility.rare_stratum_id) == strata(j);
        idxT = double(T.scale_index) == scaleIndex & T.anchor_stage == "rare_strata_enriched" & ...
            string(T.rare_stratum_id) == strata(j);
        idxTarget = double(Targets.scale_index) == scaleIndex & string(Targets.rare_stratum_id) == strata(j);
        one = table();
        one.scale_index = scaleIndex;
        one.chunk_sec = double(PrimaryScales.chunk_sec(s));
        one.rare_stratum_id = strata(j);
        one.rare_stratum_rule = local_stratum_rule(strata(j));
        one.target_rare_additions = sum(Targets.target_rare_additions(idxTarget));
        one.n_eligible_candidates = sum(Eligibility.n_eligible_candidates(idxE));
        one.n_selected_rare_anchors = nnz(idxT);
        one.n_selectable_pool_before_base_exclusion = ...
            sum(Targets.n_selectable_pool_before_base_exclusion(idxTarget));
        one.n_selectable_after_base_exclusion = ...
            sum(Targets.n_selectable_after_base_exclusion(idxTarget));
        one.n_selectable_after_prior_quota_assignments = ...
            sum(Targets.n_selectable_after_prior_quota_assignments(idxTarget));
        one.n_excluded_as_locked_base = sum(Targets.n_excluded_as_locked_base(idxTarget));
        one.n_depleted_by_prior_quota_assignments = ...
            sum(Targets.n_depleted_by_prior_quota_assignments(idxTarget));
        one.n_selected_quota_stage = sum(Targets.n_selected_quota_stage(idxTarget));
        one.n_selected_fill_stage = sum(Targets.n_selected_fill_stage(idxTarget));
        one.n_selected_exclusive_assignment = ...
            sum(Targets.n_selected_exclusive_assignment(idxTarget));
        one.n_selected_any_membership = sum(Targets.n_selected_any_membership(idxTarget));
        one.n_selected_with_multiple_memberships = ...
            sum(Targets.n_selected_with_multiple_memberships(idxTarget));
        one.quota_shortfall_after_quota_stage = ...
            sum(Targets.quota_shortfall_after_quota_stage(idxTarget));
        one.quota_shortfall_reason = strjoin(unique( ...
            string(Targets.quota_shortfall_reason(idxTarget)), 'stable'), ";");
        one.exclusive_assignment_minus_target = ...
            one.n_selected_exclusive_assignment - one.target_rare_additions;
        one.any_membership_minus_target = ...
            one.n_selected_any_membership - one.target_rare_additions;
        one.overlap_fraction_among_any_members = one.n_selected_with_multiple_memberships ./ ...
            max(one.n_selected_any_membership, 1);
        one.n_sessions_with_eligible_candidates = nnz(Eligibility.n_eligible_candidates(idxE) > 0);
        one.n_sessions_with_selected_anchors = numel(unique(T.session_index(idxT)));
        one.mean_rare_inclusion_probability = mean(T.rare_inclusion_probability(idxT), 'omitnan');
        one.labels_used_for_rare_sampling = "none";
        one.arena_used_for_rare_sampling = false;
        one.condition_used_for_rare_sampling = false;
        one.quota_audit_interpretation = ...
            "exclusive_assignment_is_priority_ordered;any_membership_is_overlap_aware;fill_is_composite_score_ranked";
        assert(one.n_selected_rare_anchors == one.n_selected_exclusive_assignment, ...
            'build_rare_enriched_primary_anchor_manifest:ExclusiveAssignmentAuditMismatch', ...
            'Exclusive assignment count does not match selected rare-anchor count.');
        assert(one.n_selected_any_membership >= one.n_selected_exclusive_assignment, ...
            'build_rare_enriched_primary_anchor_manifest:AnyMembershipAuditMismatch', ...
            'Any-membership count must be at least the exclusive assignment count.');
        Summary = [Summary; one]; %#ok<AGROW>
    end
end
end

function local_copy_baseline_audits(sourceRoot, outRoot)
files = ["primary_operational_scales.csv","embedding_dimension_audit.csv", ...
    "pca_loading_stability.csv","arena_sensitivity_audit.csv", ...
    "chunk_feature_transform_audit.csv"];
for file = files
    source = fullfile(sourceRoot, file);
    assert(isfile(source), 'build_rare_enriched_primary_anchor_manifest:MissingBaselineAudit', ...
        'Missing baseline audit required by enriched run_07: %s', source);
    copyfile(source, fullfile(outRoot, file), 'f');
end
end

function local_write_provenance_audit(outRoot, chunkRoot, graphRoot, basePath, definitionPath, seedPath)
T = table();
T.anchor_manifest_mode = "rare_enriched";
T.baseline_run06_root = string(chunkRoot);
T.baseline_run07_root = string(fullfile(fileparts(chunkRoot), 'embedding_motif_discovery'));
T.baseline_run08_root = string(graphRoot);
T.baseline_primary_anchor_manifest = string(basePath);
T.baseline_graph_node_manifest = string(fullfile(graphRoot, 'graph_node_manifest.csv'));
T.baseline_graph_edge_list = string(fullfile(graphRoot, 'graph_edge_list.csv'));
T.rare_strata_definition_csv = string(definitionPath);
T.rare_strata_seed_anchor_manifest_csv = string(seedPath);
T.labels_used_for_expanded_anchor_selection = "none";
T.arena_used_for_expanded_anchor_selection = false;
T.condition_used_for_expanded_anchor_selection = false;
writetable(T, fullfile(outRoot, 'enriched_input_provenance_audit.csv'));
end

function local_assert_label_free_inputs(Base, Event, Definition, Seed)
local_assert_none(Base, 'labels_used_for_primary_anchor_selection');
local_assert_false(Base, 'arena_used_for_primary_anchor_selection');
local_assert_false(Base, 'condition_used_for_primary_anchor_selection');
local_assert_none(Event, 'labels_used_for_event_summary');
local_assert_false(Event, 'arena_used_for_event_summary');
local_assert_false(Event, 'condition_used_for_event_summary');
local_assert_none(Definition, 'labels_used_for_rare_strata');
local_assert_false(Definition, 'arena_used_for_rare_strata');
local_assert_false(Definition, 'condition_used_for_rare_strata');
local_assert_none(Seed, 'labels_used_for_rare_strata');
local_assert_false(Seed, 'arena_used_for_rare_strata');
local_assert_false(Seed, 'condition_used_for_rare_strata');
end

function local_assert_none(T, name)
if ismember(name, T.Properties.VariableNames)
    assert(all(string(T.(name)) == "none"), ...
        'build_rare_enriched_primary_anchor_manifest:LabelLeak', '%s must equal none.', name);
end
end

function local_assert_false(T, name)
if ismember(name, T.Properties.VariableNames)
    assert(~any(logical(T.(name))), ...
        'build_rare_enriched_primary_anchor_manifest:LabelLeak', '%s must be false.', name);
end
end

function stats = local_stats_from_transform_audit(T)
if ismember('channel_index', T.Properties.VariableNames)
    T = sortrows(T, 'channel_index');
end
stats.center = double(T.robust_center(:))';
stats.scale = double(T.robust_scale(:))';
stats.impute = double(T.impute_value(:))';
stats.n_fit_rows = local_table_double(T, 'n_fit_rows_total', 1, NaN);
stats.fit_scope = local_table_string(T, 'stats_fit_scope', 1, "run06_condition_blind_transform_audit");
end

function dyad = local_load_dyad(repoRoot, featureOutputFile)
pathText = resolve_repo_path(repoRoot, featureOutputFile);
S = load(pathText, 'dyad', 'status');
assert(isfield(S, 'dyad'), 'build_rare_enriched_primary_anchor_manifest:MissingDyad', ...
    'Missing dyad in feature file: %s', pathText);
if isfield(S, 'status')
    assert(string(S.status) == "success", ...
        'build_rare_enriched_primary_anchor_manifest:BadFeatureStatus', ...
        'Feature file does not report success: %s', pathText);
end
dyad = S.dyad;
end

function mask = local_frame_mask(dyad)
if isfield(dyad, 'frameMask') && ~isempty(dyad.frameMask)
    mask = logical(dyad.frameMask(:));
else
    mask = true(size(dyad.X, 1), 1);
end
end

function x = local_feature(dyad, name)
name = string(name);
if isfield(dyad, 'raw') && isfield(dyad.raw, char(name))
    x = double(dyad.raw.(char(name))(:));
elseif isfield(dyad, char(name))
    x = double(dyad.(char(name))(:));
else
    idx = find(string(dyad.featureNames(:)) == name, 1);
    assert(~isempty(idx), 'build_rare_enriched_primary_anchor_manifest:MissingFeature', ...
        'Missing condition-blind event feature %s.', name);
    x = double(dyad.X(:, idx));
end
end

function threshold = local_definition_threshold(T, id, scaleIndex, fallback)
idx = string(T.rare_stratum_id) == string(id) & double(T.scale_index) == double(scaleIndex);
if any(idx)
    threshold = double(T.rare_stratum_threshold(find(idx, 1)));
else
    threshold = fallback;
end
end

function mask = local_even_mask(time, eligible, nWant)
mask = false(numel(eligible), 1);
idx = find(eligible);
nWant = min(numel(idx), floor(nWant));
if nWant <= 0
    return
end
if nWant == numel(idx)
    mask(idx) = true;
    return
end
[~, ord] = sort(double(time(idx)));
idx = idx(ord);
pos = unique(round(linspace(1, numel(idx), nWant)));
while numel(pos) < nWant
    missing = setdiff(1:numel(idx), pos, 'stable');
    pos(end + 1) = missing(1); %#ok<AGROW>
end
mask(idx(pos)) = true;
end

function strata = local_strata_order()
strata = ["contact_transition","contact_present","large_heading_turn", ...
    "low_density_tail","graph_periphery","high_radial_speed", ...
    "short_motif_instability","undercovered_scale_session"];
end

function q = local_quota_fractions(params)
q = [params.rare_quota_contact_transition; params.rare_quota_contact_present; ...
    params.rare_quota_large_heading_turn; params.rare_quota_low_density_tail; ...
    params.rare_quota_graph_periphery; params.rare_quota_high_radial_speed; ...
    params.rare_quota_short_motif_instability; params.rare_quota_undercovered_scale_session];
end

function id = local_primary_stratum(row, strata)
id = "undercovered_scale_session";
for j = 1:numel(strata)
    if row.("flag_" + strata(j))
        id = strata(j);
        return
    end
end
end

function rule = local_stratum_rule(id)
switch string(id)
    case "contact_transition"
        rule = "contact_transition_count_gt_zero";
    case "contact_present"
        rule = "contact_dwell_fraction_gt_zero";
    case "large_heading_turn"
        rule = "heading_large_turn_count_gt_zero";
    case "low_density_tail"
        rule = "temporal_neighborhood_of_baseline_high_knn_radius_seed";
    case "graph_periphery"
        rule = "temporal_neighborhood_of_baseline_low_in_degree_seed";
    case "high_radial_speed"
        rule = "radial_speed_mean_ge_baseline_within_scale_quantile";
    case "short_motif_instability"
        rule = "temporal_neighborhood_of_baseline_short_motif_instability_seed";
    case "undercovered_scale_session"
        rule = "baseline_scale_session_support_le_within_scale_quantile";
    otherwise
        rule = "not_rare_stratum_member";
end
end

function n = local_effective_sample_size(w)
w = double(w(:));
w = w(isfinite(w) & w > 0);
if isempty(w)
    n = 0;
else
    n = sum(w).^2 ./ sum(w.^2);
end
end

function key = local_anchor_key(scaleIndex, rawIndex, anchorFrame)
key = string(round(double(scaleIndex(:)))) + "_" + ...
    string(round(double(rawIndex(:)))) + "_" + string(round(double(anchorFrame(:))));
end

function key = local_cell_key(scaleIndex, sessionIndex)
key = string(round(double(scaleIndex(:)))) + "_" + string(round(double(sessionIndex(:))));
end

function C = local_empty_candidate_pool()
C = table();
end

function value = local_table_double(T, name, row, defaultValue)
if ismember(name, T.Properties.VariableNames)
    value = double(T.(name)(row));
else
    value = defaultValue;
end
end

function value = local_table_string(T, name, row, defaultValue)
if ismember(name, T.Properties.VariableNames)
    value = string(T.(name)(row));
else
    value = string(defaultValue);
end
end

function T = local_read_csv(pathText)
opts = detectImportOptions(pathText, 'FileType', 'text', 'Delimiter', ',', 'TextType', 'string');
T = readtable(pathText, opts);
end

function local_ensure_dir(pathText)
if ~exist(pathText, 'dir')
    mkdir(pathText);
end
end
