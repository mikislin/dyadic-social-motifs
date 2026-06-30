function notes = local_qc_notes(T, motifMask, warnThresh, maxThresh)
notes = strings(height(T), 1);
for i = 1:height(T)
    parts = strings(0,1);
    if ~T.preprocess_success(i)
        parts(end+1,1) = "preprocessing_missing_or_failed"; 
    end
    if motifMask(i) && T.badframe_fraction(i) > maxThresh
        parts(end+1,1) = "exclude_motif_high_badframe_fraction";
    elseif motifMask(i) && T.badframe_fraction(i) > warnThresh
        parts(end+1,1) = "review_motif_badframe_fraction";
    end
    if string(T.motif_exclusion_reason(i)) == "anesthetized_center_special_egocentric_context"
        parts(end+1,1) = "reserved_for_block2_anesthetized_context"; 
    end
    notes(i) = strjoin(parts, ';');
end
end