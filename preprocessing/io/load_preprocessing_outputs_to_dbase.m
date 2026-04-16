function dbase = load_preprocessing_outputs_to_dbase(dbase, outDir)
%LOAD_PREPROCESSING_OUTPUTS_TO_DBASE Merge preprocessed outputs back into dbase.

for i = 1:numel(dbase)
    inFile = fullfile(outDir, sprintf('session_%04d_preproc.mat', i));
    assert(isfile(inFile), 'Missing preprocessing output: %s', inFile);
    S = load(inFile, 'out');
    dbase(i).tracks = S.out.clean.tracks;
    dbase(i).badframes = S.out.qc.badframes;
    dbase(i).preproc_qc = S.out.qc.sessionStats;
    dbase(i).preproc_debug = S.out.debug;
    dbase(i).preproc_params = S.out.params;
end
end
