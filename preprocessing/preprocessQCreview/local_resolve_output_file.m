function pathOut = local_resolve_output_file(T, rowIdx, repoRoot, preprocRoot, cfg)
candidate = string(T.preprocess_output_file(rowIdx));
if strlength(candidate) > 0 && isfile(candidate)
    pathOut = char(candidate);
    return
end

sessionDir = fullfile(preprocRoot, char(cfg.paths.preprocess_session_subdir));
pathOut = fullfile(sessionDir, sprintf('session_%04d_preprocessed.mat', T.raw_index(rowIdx)));
if isfile(pathOut)
    return
end

if strlength(candidate) > 0
pathOut = resolve_repo_path(repoRoot, candidate);
end
end
