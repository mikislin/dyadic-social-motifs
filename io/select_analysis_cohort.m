function mask = select_analysis_cohort(manifest, cohort)
%SELECT_ANALYSIS_COHORT Return named manifest inclusion masks.
%
% Centralizing mask, keeps paper scripts from re-implementing cohort
% logic, especially the exclusion of anesthetized-center sessions from motif
% discovery while retaining them for egocentric social-context analyses.

cohort = lower(strtrim(string(cohort)));
if istable(manifest) && ismember('raw_index', manifest.Properties.VariableNames) && ...
        ~ismember('cohort_id', manifest.Properties.VariableNames)
    manifest = apply_paper_cohort_definitions(manifest);
end

switch cohort
    case {"motif", "motifs", "motif_discovery", "block1", "block1_dyadic"}
        if ismember('include_motif_discovery', manifest.Properties.VariableNames)
            mask = manifest.include_motif_discovery == 1;
        else
            mask = manifest.include_block1_dyadic == 1 & ...
                string(manifest.social_context) ~= "anesthetized_partner";
        end

    case {"block2", "block2_egocentric", "egocentric", "egocentric_context"}
        mask = manifest.include_block2_egocentric == 1;

    case {"paper", "paper_include", "preprocess", "preprocessing"}
        mask = manifest.paper_include == 1;

    case {"anesthetized", "anesthetized_partner", "anesthetized_context"}
        mask = string(manifest.social_context) == "anesthetized_partner";
        if ismember('mouse_type_1', manifest.Properties.VariableNames)
            mask = mask | string(manifest.mouse_type_1) == "ANES";
        end
        if ismember('mouse_type_2', manifest.Properties.VariableNames)
            mask = mask | string(manifest.mouse_type_2) == "ANES";
        end
        if ismember('mouse_type_3', manifest.Properties.VariableNames)
            mask = mask | string(manifest.mouse_type_3) == "ANES";
        end

    case {"cohort1", "big_primary", "block1_big_primary"}
        mask = string(manifest.cohort_id) == "cohort1";
        if any(cohort == ["big_primary","block1_big_primary"])
            mask = mask & manifest.include_defined_primary_big_contrast == 1;
        end

    case {"cohort2", "small_primary", "block1_small_primary"}
        mask = string(manifest.cohort_id) == "cohort2";
        if any(cohort == ["small_primary","block1_small_primary"])
            mask = mask & manifest.include_defined_primary_small_contrast == 1;
        end

    case {"primary", "primary_free_dyads", "block1_primary"}
        mask = manifest.include_defined_primary_big_contrast == 1 | ...
            manifest.include_defined_primary_small_contrast == 1;

    case {"wtwt_arena", "wtwt_arena_reference", "arena_wtwt"}
        mask = manifest.include_defined_wtwt_arena_contrast == 1;

    case {"cohort3", "pfc_chemogenetic", "cohort3_pfc"}
        mask = string(manifest.cohort_id) == "cohort3";
        if any(cohort == ["pfc_chemogenetic","cohort3_pfc"])
            mask = manifest.include_defined_pfc_chemogenetic == 1;
        end

    case {"dcn_chemogenetic", "cohort3_dcn"}
        mask = manifest.include_defined_dcn_chemogenetic == 1;

    otherwise
        error('select_analysis_cohort:UnknownCohort', ...
            'Unknown cohort "%s".', cohort);
end

mask = logical(mask(:));
end
