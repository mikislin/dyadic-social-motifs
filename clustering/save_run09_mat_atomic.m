function save_run09_mat_atomic(outputPath, variables, useV73)
%SAVE_RUN09_MAT_ATOMIC Save named values from a struct and replace atomically.

if nargin < 3
    useV73 = false;
end
outputPath = string(outputPath);
assert(isstruct(variables) && isscalar(variables), ...
    'save_run09_mat_atomic:BadVariables', ...
    'variables must be a scalar struct.');
parent = string(fileparts(outputPath));
if parent ~= "" && ~isfolder(parent)
    mkdir(parent);
end
temporaryPath = outputPath + ".tmp-" + ...
    string(char(javaMethod('randomUUID', 'java.util.UUID'))) + ".mat";
cleanup = onCleanup(@() i_delete_if_present(temporaryPath)); %#ok<NASGU>
if useV73
    save(temporaryPath, '-struct', 'variables', '-v7.3');
else
    save(temporaryPath, '-struct', 'variables', '-v7');
end
[moved, message] = movefile(temporaryPath, outputPath, 'f');
assert(moved, 'save_run09_mat_atomic:ReplaceFailed', ...
    'Could not atomically replace %s: %s', outputPath, message);
end

function i_delete_if_present(pathText)
if isfile(pathText)
    delete(pathText);
end
end
