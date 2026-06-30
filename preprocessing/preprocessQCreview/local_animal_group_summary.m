function summary = local_animal_group_summary(animalTable)
if isempty(animalTable) || height(animalTable) == 0
    summary = table();
    return
end
metrics = {'pctBadframes','pctPredictionIssueFrames','pctPredictionIssueRepaired', ...
    'pctPredictionIssueUnresolved','pctInterpFrames','pctLowConfFrames', ...
    'pctJumpFrames','pctGeomFrames','pctArenaFrames'};
metrics = metrics(ismember(metrics, animalTable.Properties.VariableNames));
summary = groupsummary(animalTable, {'analysis_class','arena_condition','condition'}, ...
    {'mean','median','max'}, metrics);
end