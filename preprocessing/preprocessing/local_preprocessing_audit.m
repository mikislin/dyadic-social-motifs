function audit = local_preprocessing_audit(sessionPreproc)
if ~isfield(sessionPreproc.qc, 'sessionStats') || ...
        ~isfield(sessionPreproc.qc.sessionStats, 'nPredictionIssueAnimalFrames')
    sessionPreproc.qc.sessionStats = summarize_preprocessing_qc(sessionPreproc);
end

stats = sessionPreproc.qc.sessionStats;
audit = local_empty_audit();
audit.badframe_fraction = local_badframe_fraction(sessionPreproc, stats);
audit.n_prediction_issue_animal_frames = struct_field_or_nan(stats, 'nPredictionIssueAnimalFrames');
audit.n_interpolated_prediction_issue_animal_frames = struct_field_or_nan(stats, 'nInterpolatedPredictionIssueAnimalFrames');
audit.n_usable_prediction_issue_animal_frames = struct_field_or_nan(stats, 'nUsablePredictionIssueAnimalFrames');
audit.n_repaired_prediction_issue_animal_frames = struct_field_or_nan(stats, 'nRepairedPredictionIssueAnimalFrames');
audit.n_unresolved_prediction_issue_animal_frames = struct_field_or_nan(stats, 'nUnresolvedPredictionIssueAnimalFrames');
audit.prediction_issue_fraction = struct_field_or_nan(stats, 'predictionIssueFraction');
audit.interpolated_prediction_issue_fraction = struct_field_or_nan(stats, 'interpolatedPredictionIssueFraction');
audit.usable_prediction_issue_fraction = struct_field_or_nan(stats, 'usablePredictionIssueFraction');
audit.repaired_prediction_issue_fraction = struct_field_or_nan(stats, 'repairedPredictionIssueFraction');
audit.unresolved_prediction_issue_fraction = struct_field_or_nan(stats, 'unresolvedPredictionIssueFraction');
audit.prediction_issue_usable_rate = struct_field_or_nan(stats, 'predictionIssueUsableRate');
audit.prediction_issue_repair_rate = struct_field_or_nan(stats, 'predictionIssueRepairRate');
audit.prediction_issue_unresolved_rate = struct_field_or_nan(stats, 'predictionIssueUnresolvedRate');
audit.n_prediction_issue_frames_any_animal = struct_field_or_nan(stats, 'nPredictionIssueFramesAnyAnimal');
audit.n_repaired_prediction_issue_frames_any_animal = struct_field_or_nan(stats, 'nRepairedPredictionIssueFramesAnyAnimal');
audit.n_unresolved_prediction_issue_frames_any_animal = struct_field_or_nan(stats, 'nUnresolvedPredictionIssueFramesAnyAnimal');
audit.prediction_issue_any_animal_fraction = struct_field_or_nan(stats, 'predictionIssueAnyAnimalFraction');
audit.repaired_prediction_issue_any_animal_fraction = struct_field_or_nan(stats, 'repairedPredictionIssueAnyAnimalFraction');
audit.unresolved_prediction_issue_any_animal_fraction = struct_field_or_nan(stats, 'unresolvedPredictionIssueAnyAnimalFraction');
end
