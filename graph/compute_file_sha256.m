function [value, byteCount, lineFeedCount, lastByte] = compute_file_sha256(pathText)
%COMPUTE_FILE_SHA256 SHA-256 of the exact bytes stored in a file.
%
% Bytes are read by MATLAB and passed to MessageDigest. Do not use
% FileInputStream.read with a MATLAB numeric array: Java mutates a converted
% copy rather than the MATLAB buffer, which can silently hash zero-filled
% bytes instead of file contents.

pathText = string(pathText);
assert(isscalar(pathText) && isfile(pathText), ...
    'compute_file_sha256:MissingFile', ...
    'File does not exist: %s', pathText);

[fid, message] = fopen(pathText, 'rb');
assert(fid >= 0, 'compute_file_sha256:OpenFailed', ...
    'Could not open %s: %s', pathText, message);
cleanup = onCleanup(@() fclose(fid));

md = javaMethod('getInstance', 'java.security.MessageDigest', 'SHA-256');
byteCount = 0;
lineFeedCount = 0;
lastByte = NaN;
while true
    bytes = fread(fid, 1024 * 1024, '*uint8');
    if isempty(bytes)
        break
    end
    md.update(typecast(bytes(:), 'int8'));
    byteCount = byteCount + numel(bytes);
    lineFeedCount = lineFeedCount + nnz(bytes == 10);
    lastByte = double(bytes(end));
end

digest = typecast(md.digest(), 'uint8');
value = lower(string(reshape(dec2hex(digest, 2).', 1, [])));
end
