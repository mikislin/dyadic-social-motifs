function analysisClass = local_analysis_class(T)
analysisClass = repmat("excluded_or_unsupported", height(T), 1);
analysisClass(T.include_block2_egocentric == 1) = "block2_egocentric_context";
analysisClass(T.effective_n_animals == 1 & T.include_block2_egocentric == 1) = ...
    "block2_single_mouse";
analysisClass(select_analysis_cohort(T, "anesthetized_context")) = ...
    "block2_anesthetized_center";
analysisClass(select_analysis_cohort(T, "motif_discovery")) = ...
    "block1_motif_discovery";
end