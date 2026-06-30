function [animalRows, pointRows] = local_append_session_detail_rows(animalRows, pointRows, T, rowIdx, P, cfg)
nFrames = size(P.clean.tracks, 1);
nNodes = size(P.clean.tracks, 2);
nAnimals = size(P.clean.tracks, 4);
badframes = P.qc.badframes;
if isvector(badframes)
    badframes = badframes(:);
end
stats = P.qc.sessionStats;

for animal = 1:nAnimals
    base = local_base_row(T, rowIdx);
    base.animal = animal;
    base.n_frames = nFrames;
    base.n_nodes = nNodes;

    if isfield(stats, 'animal') && numel(stats.animal) >= animal
        st = stats.animal(animal);
    else
        st = struct();
    end

    animalRow = base;
    animalRow.pctBadframes = struct_field_or_nan(st, 'pctBadframes');
    animalRow.pctInterpFrames = struct_field_or_nan(st, 'pctInterpFrames');
    animalRow.pctJumpFrames = struct_field_or_nan(st, 'pctJumpFrames');
    animalRow.pctGeomFrames = struct_field_or_nan(st, 'pctGeomFrames');
    animalRow.pctArenaFrames = struct_field_or_nan(st, 'pctArenaFrames');
    animalRow.pctLowConfFrames = struct_field_or_nan(st, 'pctLowConfFrames');
    animalRow.pctJumpSamples = struct_field_or_nan(st, 'pctJumpSamples');
    animalRow.pctInterpSamples = struct_field_or_nan(st, 'pctInterpSamples');
    animalRow.nPredictionIssueFrames = struct_field_or_nan(st, 'nPredictionIssueFrames');
    animalRow.nInterpolatedPredictionIssueFrames = struct_field_or_nan(st, 'nInterpolatedPredictionIssueFrames');
    animalRow.nUsablePredictionIssueFrames = struct_field_or_nan(st, 'nUsablePredictionIssueFrames');
    animalRow.nRepairedPredictionIssueFrames = struct_field_or_nan(st, 'nRepairedPredictionIssueFrames');
    animalRow.nUnresolvedPredictionIssueFrames = struct_field_or_nan(st, 'nUnresolvedPredictionIssueFrames');
    animalRow.pctPredictionIssueFrames = struct_field_or_nan(st, 'pctPredictionIssueFrames');
    animalRow.pctPredictionIssueUsable = struct_field_or_nan(st, 'pctPredictionIssueUsable');
    animalRow.pctPredictionIssueRepaired = struct_field_or_nan(st, 'pctPredictionIssueRepaired');
    animalRow.pctPredictionIssueUnresolved = struct_field_or_nan(st, 'pctPredictionIssueUnresolved');
    animalRows = [animalRows; animalRow]; %#ok<AGROW>

    qa = P.qc.animals(animal);
    masks = local_full_qc_masks(qa, nFrames, nNodes);
    animalBad = badframes(:, min(animal, size(badframes, 2)));
    for node = 1:nNodes
        issue = masks.lowConf(:,node) | masks.jump(:,node) | masks.geom(:,node);
        repaired = issue & masks.interp(:,node) & ~animalBad & ~masks.finalNan(:,node);
        unresolved = issue & (animalBad | masks.finalNan(:,node));

        pointRow = base;
        pointRow.node = node;
        pointRow.node_label = local_node_label(cfg, node);
        pointRow.n_prediction_issue_samples = nnz(issue);
        pointRow.n_interpolated_prediction_issue_samples = nnz(issue & masks.interp(:,node));
        pointRow.n_repaired_prediction_issue_samples = nnz(repaired);
        pointRow.n_unresolved_prediction_issue_samples = nnz(unresolved);
        pointRow.fraction_low_confidence_samples = mean(masks.lowConf(:,node));
        pointRow.fraction_jump_samples = mean(masks.jump(:,node));
        pointRow.fraction_geometry_samples = mean(masks.geom(:,node));
        pointRow.fraction_arena_samples = mean(masks.arena(:,node));
        pointRow.fraction_interpolated_samples = mean(masks.interp(:,node));
        pointRow.fraction_final_nan_samples = mean(masks.finalNan(:,node));
        pointRow.fraction_prediction_issue_samples = mean(issue);
        pointRow.fraction_repaired_prediction_issue_samples = mean(repaired);
        pointRow.fraction_unresolved_prediction_issue_samples = mean(unresolved);
        pointRow.prediction_issue_repair_rate = local_fraction(nnz(repaired), nnz(issue));
        pointRow.prediction_issue_unresolved_rate = local_fraction(nnz(unresolved), nnz(issue));
        pointRows = [pointRows; pointRow]; %#ok<AGROW>
    end
end
end
