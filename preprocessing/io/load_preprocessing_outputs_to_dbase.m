function dbase = load_preprocessing_outputs_to_dbase(dbase, outDir)
%LOAD_PREPROCESSING_OUTPUTS_TO_DBASE Merge preprocessed outputs back into dbase.
%
% Legacy helper for older dbase-struct workflows. The paper pipeline stores
% preprocessing metadata in derived/preprocessed summaries instead.

for i = 1:numel(dbase)
    inFile = fullfile(outDir, sprintf('session_%04d_preproc.mat', i));
    assert(isfile(inFile), 'Missing preprocessing output: %s', inFile);
    S = load(inFile);
    if isfield(S, 'out')
        sessionPreproc = S.out;
    elseif isfield(S, 'sessionPreproc')
        sessionPreproc = S.sessionPreproc;
    else
        error('load_preprocessing_outputs_to_dbase:MissingOutput', ...
            'Missing out or sessionPreproc variable in %s.', inFile);
    end
    dbase(i).preprocessed_tracks = sessionPreproc.clean.tracks;
    dbase(i).preprocessing_badframes = sessionPreproc.qc.badframes;
    dbase(i).tracks = sessionPreproc.clean.tracks;
    dbase(i).badframes = sessionPreproc.qc.badframes;
    dbase(i).preproc_qc = sessionPreproc.qc.sessionStats;
    dbase(i).preproc_debug = sessionPreproc.debug;
    dbase(i).preproc_params = sessionPreproc.params;
end
end
