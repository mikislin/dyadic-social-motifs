function Audit = persist_run08_rare_definition_stages(repoRoot, outRoot, params, ...
    postfitDefinition, postfitMembership, postfitSeed)
%PERSIST_RUN08_RARE_DEFINITION_STAGES Separate selection-time and post-fit strata.
%
% The locked baseline definition used to select enriched anchors is preserved
% verbatim in explicit selection-time files. Current-graph re-evaluation is
% written to separate post-fit files and cannot be mistaken for the selection
% rule that generated the enriched bank.

repoRoot = string(repoRoot);
outRoot = string(outRoot);
lockedRoot = resolve_repo_path(repoRoot, params.baseline_graph_input_dir);
lockedDefinitionPath = fullfile(lockedRoot, 'rare_strata_definition.csv');
lockedMembershipPath = fullfile(lockedRoot, 'rare_strata_node_membership.csv');
lockedSeedPath = fullfile(lockedRoot, 'rare_strata_seed_anchor_manifest.csv');

LockedDefinition = local_optional_read(lockedDefinitionPath);
LockedMembership = local_optional_read(lockedMembershipPath);
LockedSeed = local_optional_read(lockedSeedPath);

PostDefinition = postfitDefinition;
PostDefinition.definition_stage = repmat("postfit_current_graph_re_evaluation", height(PostDefinition), 1);
PostDefinition.definition_source_graph_root = repmat(outRoot, height(PostDefinition), 1);
PostDefinition.active_for_anchor_selection = false(height(PostDefinition), 1);
PostDefinition.recomputed_after_enriched_graph_fit = ...
    repmat(string(params.anchor_manifest_mode) == "rare_enriched", height(PostDefinition), 1);
PostMembership = postfitMembership;
PostMembership.definition_stage = repmat("postfit_current_graph_re_evaluation", height(PostMembership), 1);
PostMembership.active_for_anchor_selection = false(height(PostMembership), 1);
PostSeed = postfitSeed;
PostSeed.definition_stage = repmat("postfit_current_graph_re_evaluation", height(PostSeed), 1);
PostSeed.active_for_anchor_selection = false(height(PostSeed), 1);

writetable(PostDefinition, fullfile(outRoot, 'rare_strata_postfit_definition.csv'));
writetable(PostMembership, fullfile(outRoot, 'rare_strata_postfit_node_membership.csv'));
writetable(PostSeed, fullfile(outRoot, 'rare_strata_postfit_seed_anchor_manifest.csv'));

% Backward-compatible names remain explicitly marked as post-fit artifacts.
writetable(PostDefinition, fullfile(outRoot, 'rare_strata_definition.csv'));
writetable(PostMembership, fullfile(outRoot, 'rare_strata_node_membership.csv'));
writetable(PostSeed, fullfile(outRoot, 'rare_strata_seed_anchor_manifest.csv'));

if ~isempty(LockedDefinition)
    LockedDefinition.definition_stage = repmat("locked_selection_time_baseline", height(LockedDefinition), 1);
    LockedDefinition.definition_source_graph_root = repmat(lockedRoot, height(LockedDefinition), 1);
    LockedDefinition.active_for_anchor_selection = true(height(LockedDefinition), 1);
    LockedDefinition.recomputed_after_enriched_graph_fit = false(height(LockedDefinition), 1);
    writetable(LockedDefinition, fullfile(outRoot, 'rare_strata_selection_definition_locked.csv'));
end
if ~isempty(LockedMembership)
    LockedMembership.definition_stage = repmat("locked_selection_time_baseline", height(LockedMembership), 1);
    LockedMembership.active_for_anchor_selection = true(height(LockedMembership), 1);
    writetable(LockedMembership, fullfile(outRoot, 'rare_strata_selection_node_membership_locked.csv'));
end
if ~isempty(LockedSeed)
    LockedSeed.definition_stage = repmat("locked_selection_time_baseline", height(LockedSeed), 1);
    LockedSeed.active_for_anchor_selection = true(height(LockedSeed), 1);
    writetable(LockedSeed, fullfile(outRoot, 'rare_strata_selection_seed_anchor_manifest_locked.csv'));
end

Provenance = table();
Provenance = [Provenance; local_provenance_row( ...
    "locked_selection_time_baseline", lockedRoot, lockedDefinitionPath, ...
    height(LockedDefinition), true, false, params)];
Provenance = [Provenance; local_provenance_row( ...
    "postfit_current_graph_re_evaluation", outRoot, ...
    fullfile(outRoot, 'rare_strata_postfit_definition.csv'), ...
    height(PostDefinition), false, true, params)];
writetable(Provenance, fullfile(outRoot, 'rare_strata_definition_provenance_audit.csv'));

Comparison = local_definition_comparison(LockedDefinition, PostDefinition);
writetable(Comparison, fullfile(outRoot, 'rare_strata_definition_stage_comparison_audit.csv'));

Audit = struct();
Audit.provenance = Provenance;
Audit.comparison = Comparison;
Audit.lockedDefinition = LockedDefinition;
Audit.postfitDefinition = PostDefinition;
end

function T = local_optional_read(pathText)
if isfile(pathText)
    T = readtable(pathText, 'TextType', 'string');
else
    T = table();
end
end

function row = local_provenance_row(stage, root, pathText, nRows, active, recomputed, params)
row = table();
row.definition_stage = string(stage);
row.definition_source_graph_root = string(root);
row.definition_csv = string(pathText);
row.n_definition_rows = nRows;
row.active_for_anchor_selection = logical(active);
row.recomputed_after_current_graph_fit = logical(recomputed);
row.anchor_manifest_mode = string(params.anchor_manifest_mode);
row.stage_interpretation = local_stage_interpretation(stage);
row.labels_used_for_definition = "none";
row.arena_used_for_definition = false;
row.condition_used_for_definition = false;
end

function text = local_stage_interpretation(stage)
if string(stage) == "locked_selection_time_baseline"
    text = "definition_available_before_enriched_anchor_selection";
else
    text = "diagnostic_re_evaluation_not_a_retrospective_selection_rule";
end
end

function C = local_definition_comparison(Locked, Post)
C = table();
if isempty(Post)
    return
end
for i = 1:height(Post)
    one = table();
    one.rare_stratum_id = string(Post.rare_stratum_id(i));
    one.scale_index = double(Post.scale_index(i));
    one.chunk_sec = double(Post.chunk_sec(i));
    one.postfit_threshold = double(Post.rare_stratum_threshold(i));
    one.postfit_secondary_threshold = double(Post.rare_stratum_secondary_threshold(i));
    one.postfit_tertiary_threshold = double(Post.rare_stratum_tertiary_threshold(i));
    one.postfit_n_member_nodes = double(Post.n_member_nodes(i));
    one.locked_threshold = NaN;
    one.locked_secondary_threshold = NaN;
    one.locked_tertiary_threshold = NaN;
    one.locked_n_member_nodes = NaN;
    if ~isempty(Locked) && all(ismember({'rare_stratum_id','scale_index'}, Locked.Properties.VariableNames))
        idx = string(Locked.rare_stratum_id) == one.rare_stratum_id & ...
            double(Locked.scale_index) == one.scale_index;
        if any(idx)
            j = find(idx, 1);
            one.locked_threshold = double(Locked.rare_stratum_threshold(j));
            one.locked_secondary_threshold = double(Locked.rare_stratum_secondary_threshold(j));
            one.locked_tertiary_threshold = double(Locked.rare_stratum_tertiary_threshold(j));
            one.locked_n_member_nodes = double(Locked.n_member_nodes(j));
        end
    end
    one.threshold_delta_postfit_minus_locked = one.postfit_threshold - one.locked_threshold;
    one.member_count_delta_postfit_minus_locked = one.postfit_n_member_nodes - one.locked_n_member_nodes;
    one.comparison_role = "selection_time_definition_vs_postfit_re_evaluation_audit";
    one.active_definition_for_enriched_selection = "locked_selection_time_baseline";
    one.labels_used_for_comparison = "none";
    one.arena_used_for_comparison = false;
    one.condition_used_for_comparison = false;
    C = [C; one]; %#ok<AGROW>
end
end
