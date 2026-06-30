function relPath = local_repo_relative_path(repoRoot, absPath)
repoRoot = char(repoRoot);
absPath = char(absPath);
if startsWith(absPath, repoRoot)
    relPath = string(extractAfter(string(absPath), strlength(string(repoRoot)) + 1));
else
    relPath = string(absPath);
end
relPath = replace(relPath, '\', '/');
end
