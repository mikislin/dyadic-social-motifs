function pathOut = resolve_repo_path(repoRoot, pathText)
%RESOLVE_REPO_PATH Return an absolute path for a repo-relative config value.

pathText = char(string(pathText));
if startsWith(pathText, filesep) || ~isempty(regexp(pathText, '^[A-Za-z]:[\\/]', 'once'))
    pathOut = pathText;
else
    pathOut = fullfile(char(string(repoRoot)), pathText);
end
end
