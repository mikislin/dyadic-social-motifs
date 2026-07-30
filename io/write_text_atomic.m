function write_text_atomic(lines, outputPath)
%WRITE_TEXT_ATOMIC Write UTF-8 text then atomically replace its target.

outputPath = string(outputPath);
parent = string(fileparts(outputPath));
if parent ~= "" && ~isfolder(parent)
    mkdir(parent);
end
temporaryPath = outputPath + ".tmp-" + ...
    string(char(javaMethod('randomUUID','java.util.UUID')));
cleanup = onCleanup(@() i_delete_if_present(temporaryPath)); %#ok<NASGU>
[fid,message] = fopen(temporaryPath,'w','n','UTF-8');
assert(fid>=0,'write_text_atomic:OpenFailed', ...
    'Could not open %s: %s',temporaryPath,message);
fileCleanup = onCleanup(@() fclose(fid));
lines = string(lines);
for i = 1:numel(lines)
    fprintf(fid,'%s\n',lines(i));
end
clear fileCleanup
[moved,message] = movefile(temporaryPath,outputPath,'f');
assert(moved,'write_text_atomic:ReplaceFailed', ...
    'Could not atomically replace %s: %s',outputPath,message);
end

function i_delete_if_present(pathText)
if isfile(pathText)
    delete(pathText);
end
end
