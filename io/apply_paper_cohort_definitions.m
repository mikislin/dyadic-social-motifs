function manifest = apply_paper_cohort_definitions(manifest)
%APPLY_PAPER_COHORT_DEFINITIONS Add cohort and contrast annotations.
%
% These columns prevent reuse of bare condition labels across cohorts or
% arena sizes. 

assert(istable(manifest), 'apply_paper_cohort_definitions:InputNotTable', ...
    'Input must be a table.');

n = height(manifest);
vars = string(manifest.Properties.VariableNames);

if ~ismember("raw_index", vars)
    manifest.raw_index = (1:n)';
end
if ~ismember("condition", vars)
    manifest.condition = strings(n,1);
end
if ~ismember("arena", vars)
    manifest.arena = strings(n,1);
end
if ~ismember("arena_label", vars)
    manifest.arena_label = strings(n,1);
end
if ~ismember("session_file", vars)
    manifest.session_file = strings(n,1);
end

rawIndex = double(manifest.raw_index);
condition = string(manifest.condition);
arena = string(manifest.arena);
arenaLabel = string(manifest.arena_label);
missingArenaLabel = arenaLabel == "";
arenaLabel(missingArenaLabel) = local_arena_label(arena(missingArenaLabel));

sessionFile = string(manifest.session_file);
sessionRelPath = sessionFile;
if ismember("session_rel_path", vars)
    hasRel = strlength(string(manifest.session_rel_path)) > 0;
    sessionRelPath(hasRel) = string(manifest.session_rel_path(hasRel));
end

manifest.session_path = sessionRelPath;
manifest.session_rel_path = sessionRelPath;
manifest.session_source_dir = repmat("files2run", n, 1);
manifest.raw_database_file = repmat("ALLPAIRS-Feb2025.mat", n, 1);
manifest.manifest_schema_version = repmat("foundation_v2_label_safe", n, 1);
manifest.manifest_generated_by = repmat("paper/run_01_manifest.m", n, 1);

manifest.condition_id = condition;
manifest.condition_raw = condition;
manifest.dyad_condition = condition;
manifest.condition_label = local_condition_label(condition);
manifest.dyad_condition_label = manifest.condition_label;
manifest.condition_family_id = local_condition_family_id(condition);
manifest.condition_family = manifest.condition_family_id;
manifest.condition_family_label = local_condition_family_label(manifest.condition_family_id);
manifest.mouse_pair_label = local_mouse_pair_label(manifest, condition);
manifest.arena_display_label = local_arena_display_label(arenaLabel);
manifest.social_context_label = local_code_label(local_column_or_default(manifest, ...
    "social_context", strings(n,1)));
manifest.drug_label = local_code_label(local_column_or_default(manifest, ...
    "drug", strings(n,1)));
manifest.dreadd_target_label = local_code_label(local_column_or_default(manifest, ...
    "dreadd_target", strings(n,1)));

manifest.cohort_id = repmat("unassigned", n, 1);
manifest.cohort_label = repmat("unassigned", n, 1);
manifest.cohort_short_label = repmat("unassigned", n, 1);
manifest.cohort_condition = strings(n,1);
manifest.cohort_condition_label = strings(n,1);
manifest.cohort_definition_source = repmat("user_2026_06_25", n, 1);
manifest.cohort_definition_note = strings(n,1);
manifest.defined_contrast_set = repmat("unassigned", n, 1);
manifest.defined_contrast_label = repmat("Unassigned", n, 1);

manifest.include_defined_primary_big_contrast = false(n,1);
manifest.include_defined_primary_small_contrast = false(n,1);
manifest.include_defined_wtwt_arena_contrast = false(n,1);
manifest.include_defined_anesthetized_context = false(n,1);
manifest.include_defined_pfc_chemogenetic = false(n,1);
manifest.include_defined_dcn_chemogenetic = false(n,1);

defs = local_definitions();
for d = 1:numel(defs)
    mask = local_range_mask(rawIndex, defs(d).ranges);
    manifest.cohort_id(mask) = defs(d).cohort_id;
    manifest.cohort_label(mask) = defs(d).cohort_label;
    manifest.cohort_short_label(mask) = defs(d).cohort_short_label;
    manifest.defined_contrast_set(mask) = defs(d).contrast_set;
    manifest.defined_contrast_label(mask) = defs(d).contrast_label;
    manifest.cohort_definition_note(mask) = defs(d).note;

    if defs(d).expected_condition ~= ""
        mismatch = mask & condition ~= defs(d).expected_condition;
        manifest.cohort_definition_note(mismatch) = ...
            manifest.cohort_definition_note(mismatch) + ...
            ";condition_mismatch_expected_" + defs(d).expected_condition;
    end
end

manifest.cohort_condition = manifest.cohort_id + ":" + condition;
manifest.cohort_condition_label = manifest.cohort_short_label + " " + ...
    manifest.condition_label;
manifest.arena_condition = arenaLabel + ":" + condition;
manifest.arena_condition_label = manifest.arena_display_label + " " + ...
    manifest.condition_label;
manifest.analysis_group_id = manifest.cohort_id + "|" + arenaLabel + "|" + condition;
manifest.analysis_group_label = manifest.cohort_label;
manifest.plot_group_label = local_plot_group_label(manifest.cohort_short_label, ...
    arenaLabel, manifest.condition_label);
manifest.label_text_safe_for_matlab = local_label_text_is_safe(manifest);

primary = ismember(condition, ["WT_WT","M_WT","M_M"]);
manifest.include_defined_primary_big_contrast = ...
    manifest.cohort_id == "cohort1" & arenaLabel == "big" & primary;
manifest.include_defined_primary_small_contrast = ...
    manifest.cohort_id == "cohort2" & arenaLabel == "small" & primary;
manifest.include_defined_wtwt_arena_contrast = ...
    ((manifest.cohort_id == "cohort1" & arenaLabel == "big") | ...
     (manifest.cohort_id == "cohort2" & arenaLabel == "small")) & ...
    condition == "WT_WT";
manifest.include_defined_anesthetized_context = ...
    manifest.cohort_id == "cohort1" & ismember(condition, ["M_A","WT_A"]);
manifest.include_defined_pfc_chemogenetic = ...
    manifest.cohort_id == "cohort3" & contains(condition, "PFC");
manifest.include_defined_dcn_chemogenetic = ...
    manifest.cohort_id == "cohort3_dcn_extension";

manifest.arena_label = arenaLabel;
end

function defs = local_definitions()
defs = struct('cohort_id', {}, 'cohort_label', {}, 'expected_condition', {}, ...
    'ranges', {}, 'contrast_set', {}, 'note', {}, ...
    'cohort_short_label', {}, 'contrast_label', {});

defs(end+1) = local_def("cohort1_single_mouse", "Cohort 1 big arena WT single baseline", ...
    "WT", [1 20], "single_mouse_baseline", ...
    "single-mouse baseline retained for egocentric/social-context analyses", ...
    "C1", "Single mouse baseline");

defs(end+1) = local_def("cohort1", "Cohort 1 big arena WT Pair", ...
    "WT_WT", [21 31], "block1_big_free_dyad_primary", "", ...
    "C1", "Block 1 big arena primary dyads");
defs(end+1) = local_def("cohort1", "Cohort 1 big arena Mixed Pair", ...
    "M_WT", [32 41], "block1_big_free_dyad_primary", "", ...
    "C1", "Block 1 big arena primary dyads");
defs(end+1) = local_def("cohort1", "Cohort 1 big arena Mut Pair", ...
    "M_M", [42 52], "block1_big_free_dyad_primary", "", ...
    "C1", "Block 1 big arena primary dyads");
defs(end+1) = local_def("cohort1", "Cohort 1 big arena Mut with anesthetized center", ...
    "M_A", [53 66], "anesthetized_center_context", ...
    "reserved for special egocentric/social-context analyses", ...
    "C1", "Anesthetized center context");
defs(end+1) = local_def("cohort1", "Cohort 1 big arena WT with anesthetized center", ...
    "WT_A", [67 75], "anesthetized_center_context", ...
    "reserved for special egocentric/social-context analyses", ...
    "C1", "Anesthetized center context");

defs(end+1) = local_def("cohort2", "Cohort 2 small arena WT Pair", ...
    "WT_WT", [76 93], "block1_small_free_dyad_primary", "", ...
    "C2", "Block 1 small arena primary dyads");
defs(end+1) = local_def("cohort2", "Cohort 2 small arena Mixed Pair", ...
    "M_WT", [94 117], "block1_small_free_dyad_primary", "", ...
    "C2", "Block 1 small arena primary dyads");
defs(end+1) = local_def("cohort2", "Cohort 2 small arena Mut Pair", ...
    "M_M", [118 129], "block1_small_free_dyad_primary", "", ...
    "C2", "Block 1 small arena primary dyads");

defs(end+1) = local_def("cohort3", "Cohort 3 small arena WT DCZ PFC", ...
    "WT_DCZPFC", [130 139], "cohort3_pfc_chemogenetic", "", ...
    "C3", "Cohort 3 PFC chemogenetic contrast");
defs(end+1) = local_def("cohort3", "Cohort 3 small arena WT Pair reference early", ...
    "WT_WT", [140 149], "cohort3_pfc_chemogenetic", "", ...
    "C3", "Cohort 3 PFC chemogenetic contrast");
defs(end+1) = local_def("cohort3", "Cohort 3 small arena WT Pair reference late", ...
    "WT_WT", [160 189], "cohort3_pfc_chemogenetic", "", ...
    "C3", "Cohort 3 PFC chemogenetic contrast");
defs(end+1) = local_def("cohort3", "Cohort 3 small arena Mut DCZ PFC", ...
    "M_DCZPFC", [190 199], "cohort3_pfc_chemogenetic", "", ...
    "C3", "Cohort 3 PFC chemogenetic contrast");
defs(end+1) = local_def("cohort3", "Cohort 3 small arena Mut DCZ control PFC", ...
    "M_DCZcPFC", [200 204], "cohort3_pfc_chemogenetic", ...
    "control partner PFC chemogenetic group", ...
    "C3", "Cohort 3 PFC chemogenetic contrast");
defs(end+1) = local_def("cohort3", "Cohort 3 small arena Mut SAL PFC", ...
    "M_SALPFC", [219 228], "cohort3_pfc_chemogenetic", "", ...
    "C3", "Cohort 3 PFC chemogenetic contrast");
defs(end+1) = local_def("cohort3", "Cohort 3 small arena Mut SAL control PFC", ...
    "M_SALcPFC", [229 233], "cohort3_pfc_chemogenetic", ...
    "control partner PFC chemogenetic group", ...
    "C3", "Cohort 3 PFC chemogenetic contrast");

defs(end+1) = local_def("cohort3_dcn_extension", "Cohort 3 DCN extension small arena WT DCZ DCN", ...
    "WT_DCZDCN", [150 159], "cohort3_dcn_chemogenetic", ...
    "DCN chemogenetic group present in raw data; documented outside PFC contrast", ...
    "C3 DCN", "Cohort 3 DCN chemogenetic extension");
defs(end+1) = local_def("cohort3_dcn_extension", "Cohort 3 DCN extension small arena Mut DCZ control DCN", ...
    "M_DCZcDCN", [205 208], "cohort3_dcn_chemogenetic", ...
    "DCN chemogenetic group present in raw data; documented outside PFC contrast", ...
    "C3 DCN", "Cohort 3 DCN chemogenetic extension");
defs(end+1) = local_def("cohort3_dcn_extension", "Cohort 3 DCN extension small arena Mut DCZ DCN", ...
    "M_DCZDCN", [209 218], "cohort3_dcn_chemogenetic", ...
    "DCN chemogenetic group present in raw data; documented outside PFC contrast", ...
    "C3 DCN", "Cohort 3 DCN chemogenetic extension");
defs(end+1) = local_def("cohort3_dcn_extension", "Cohort 3 DCN extension small arena Mut SAL control DCN", ...
    "M_SALcDCN", [234 237], "cohort3_dcn_chemogenetic", ...
    "DCN chemogenetic group present in raw data; documented outside PFC contrast", ...
    "C3 DCN", "Cohort 3 DCN chemogenetic extension");
defs(end+1) = local_def("cohort3_dcn_extension", "Cohort 3 DCN extension small arena Mut SAL DCN", ...
    "M_SALDCN", [238 247], "cohort3_dcn_chemogenetic", ...
    "DCN chemogenetic group present in raw data; documented outside PFC contrast", ...
    "C3 DCN", "Cohort 3 DCN chemogenetic extension");
end

function def = local_def(cohortId, cohortLabel, expectedCondition, ranges, ...
    contrastSet, note, cohortShortLabel, contrastLabel)
def = struct();
def.cohort_id = cohortId;
def.cohort_label = cohortLabel;
def.expected_condition = expectedCondition;
def.ranges = ranges;
def.contrast_set = contrastSet;
def.note = note;
def.cohort_short_label = cohortShortLabel;
def.contrast_label = contrastLabel;
end

function mask = local_range_mask(rawIndex, ranges)
mask = false(size(rawIndex));
for r = 1:size(ranges,1)
    mask = mask | (rawIndex >= ranges(r,1) & rawIndex <= ranges(r,2));
end
end

function label = local_arena_label(code)
code = upper(strtrim(string(code)));
label = strings(size(code));
label(code == "B") = "big";
label(code == "S") = "small";
end

function label = local_condition_label(condition)
condition = string(condition);
label = replace(condition, "_", " ");
for i = 1:numel(condition)
    switch condition(i)
        case "WT"
            label(i) = "WT single";
        case "WT_WT"
            label(i) = "WT Pair";
        case "M_WT"
            label(i) = "Mixed Pair";
        case "M_M"
            label(i) = "Mut Pair";
        case "M_A"
            label(i) = "Mut + anesthetized center";
        case "WT_A"
            label(i) = "WT + anesthetized center";
        case "WT_DCZPFC"
            label(i) = "WT DCZ PFC";
        case "M_DCZPFC"
            label(i) = "Mut DCZ PFC";
        case "M_SALPFC"
            label(i) = "Mut SAL PFC";
        case "M_SALcPFC"
            label(i) = "Mut SAL control PFC";
        case "M_DCZcPFC"
            label(i) = "Mut DCZ control PFC";
        case "WT_DCZDCN"
            label(i) = "WT DCZ DCN";
        case "M_DCZDCN"
            label(i) = "Mut DCZ DCN";
        case "M_SALDCN"
            label(i) = "Mut SAL DCN";
        case "M_SALcDCN"
            label(i) = "Mut SAL control DCN";
        case "M_DCZcDCN"
            label(i) = "Mut DCZ control DCN";
    end
end
end

function family = local_condition_family_id(condition)
condition = string(condition);
family = repmat("other", size(condition));
family(condition == "WT") = "single_mouse";
family(ismember(condition, ["WT_WT","M_WT","M_M"])) = "freely_moving_dyad";
family(ismember(condition, ["M_A","WT_A"])) = "anesthetized_center";
family(contains(condition, "PFC")) = "pfc_chemogenetic";
family(contains(condition, "DCN")) = "dcn_chemogenetic";
end

function label = local_condition_family_label(familyId)
familyId = string(familyId);
label = local_code_label(familyId);
label(familyId == "single_mouse") = "Single mouse";
label(familyId == "freely_moving_dyad") = "Freely moving dyad";
label(familyId == "anesthetized_center") = "Anesthetized center";
label(familyId == "pfc_chemogenetic") = "PFC chemogenetic";
label(familyId == "dcn_chemogenetic") = "DCN chemogenetic";
label(familyId == "other") = "Other";
end

function label = local_mouse_pair_label(manifest, condition)
condition = string(condition);
label = local_condition_label(condition);

vars = string(manifest.Properties.VariableNames);
if all(ismember(["mouse_type_1","mouse_type_2"], vars))
    m1 = string(manifest.mouse_type_1);
    m2 = string(manifest.mouse_type_2);
    for i = 1:numel(condition)
        if m1(i) ~= "" && m2(i) ~= ""
            label(i) = local_mouse_type_label(m1(i)) + " + " + ...
                local_mouse_type_label(m2(i));
        elseif m1(i) ~= ""
            label(i) = local_mouse_type_label(m1(i));
        end
    end
end

label(condition == "WT_WT") = "WT Pair";
label(condition == "M_WT") = "Mixed Pair";
label(condition == "M_M") = "Mut Pair";
end

function label = local_mouse_type_label(mouseType)
mouseType = string(mouseType);
label = local_code_label(mouseType);
label(mouseType == "WT") = "WT";
label(mouseType == "MUT") = "Mut";
label(mouseType == "ANES") = "Anesthetized center";
label(mouseType == "WT_DREADD_PFC") = "WT DREADD PFC";
label(mouseType == "WT_DREADD_DCN") = "WT DREADD DCN";
end

function label = local_arena_display_label(arenaLabel)
arenaLabel = string(arenaLabel);
label = local_code_label(arenaLabel);
label(arenaLabel == "big") = "Big arena";
label(arenaLabel == "small") = "Small arena";
end

function label = local_plot_group_label(cohortShortLabel, arenaLabel, conditionLabel)
label = strtrim(cohortShortLabel + " " + arenaLabel + " " + conditionLabel);
label = local_sanitize_label(label);
end

function value = local_column_or_default(T, name, defaultValue)
if ismember(name, string(T.Properties.VariableNames))
    value = string(T.(name));
else
    value = defaultValue;
end
end

function label = local_code_label(value)
value = string(value);
label = replace(value, "_", " ");
label = strtrim(label);
for i = 1:numel(label)
    if strlength(label(i)) == 0
        continue
    end
    words = split(label(i));
    for w = 1:numel(words)
        token = words(w);
        if all(isstrprop(char(token), 'upper')) || any(isstrprop(char(token), 'digit'))
            words(w) = token;
        else
            words(w) = upper(extractBetween(token, 1, 1)) + lower(extractAfter(token, 1));
        end
    end
    label(i) = strjoin(words, " ");
end
label(value == "none") = "None";
label(value == "single") = "Single";
label(value == "freely_moving_dyad") = "Freely moving dyad";
label(value == "anesthetized_partner") = "Anesthetized partner";
label(value == "chemogenetic_context") = "Chemogenetic context";
end

function label = local_sanitize_label(label)
label = replace(string(label), "_", " ");
label = regexprep(label, '\s+', ' ');
label = strtrim(label);
end

function ok = local_label_text_is_safe(manifest)
labelVars = ["condition_label","dyad_condition_label","condition_family_label", ...
    "mouse_pair_label","arena_display_label","social_context_label", ...
    "drug_label","dreadd_target_label","cohort_label", ...
    "cohort_short_label","cohort_condition_label","defined_contrast_label", ...
    "arena_condition_label","analysis_group_label","plot_group_label"];
ok = true(height(manifest),1);
vars = string(manifest.Properties.VariableNames);
for v = labelVars(ismember(labelVars, vars))
    ok = ok & ~contains(string(manifest.(v)), "_");
end
end
