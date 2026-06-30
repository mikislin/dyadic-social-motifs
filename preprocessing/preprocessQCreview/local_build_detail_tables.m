function [animalTable, pointTable, detailAudit] = local_build_detail_tables(T, repoRoot, preprocRoot, cfg)
animalRows = struct([]);
pointRows = struct([]);
loaded = false(height(T), 1);
loadMessage = strings(height(T), 1);

for i = 1:height(T)
    if ~T.preprocess_success(i)
        loadMessage(i) = "preprocessing_not_successful";
        continue
    end

    outFile = local_resolve_output_file(T, i, repoRoot, preprocRoot, cfg);
    if ~isfile(outFile)
        loadMessage(i) = "missing_preprocessed_mat";
        continue
    end

    try
        S = load(outFile, 'sessionPreproc');
        assert(isfield(S, 'sessionPreproc'), 'Missing sessionPreproc');
        P = S.sessionPreproc;
        if ~isfield(P.qc, 'sessionStats') || ...
                ~isfield(P.qc.sessionStats, 'nPredictionIssueAnimalFrames')
            P.qc.sessionStats = summarize_preprocessing_qc(P);
        end
        [animalRows, pointRows] = local_append_session_detail_rows( ...
            animalRows, pointRows, T, i, P, cfg);
        loaded(i) = true;
        loadMessage(i) = "ok";
    catch ME
        loadMessage(i) = string(ME.message);
    end
end

animalTable = local_struct_to_table(animalRows);
pointTable = local_struct_to_table(pointRows);
detailAudit = table(T.raw_index, T.session_id, loaded, loadMessage, ...
    'VariableNames', {'raw_index','session_id','loaded_preprocessed_mat','load_message'});
end