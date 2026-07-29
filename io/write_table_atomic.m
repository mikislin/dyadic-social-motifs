function write_table_atomic(T, outputPath)
%WRITE_TABLE_ATOMIC Write a table then atomically replace its target.

outputPath = string(outputPath);
parent = string(fileparts(outputPath));
if parent ~= "" && ~isfolder(parent)
    mkdir(parent);
end
temporaryPath = outputPath + ".tmp-" + ...
    string(char(javaMethod('randomUUID', 'java.util.UUID'))) + ".csv";
cleanup = onCleanup(@() i_delete_if_present(temporaryPath)); %#ok<NASGU>
writetable(T, temporaryPath);
[moved, message] = movefile(temporaryPath, outputPath, 'f');
assert(moved, 'write_table_atomic:ReplaceFailed', ...
    'Could not atomically replace %s: %s', outputPath, message);
end

function i_delete_if_present(pathText)
if isfile(pathText)
    delete(pathText);
end
end
