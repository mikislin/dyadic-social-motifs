function value = compute_run09_sha256_text(textValue)
%COMPUTE_RUN09_SHA256_TEXT SHA-256 of an explicitly UTF-8 encoded string.

textValue = string(textValue);
assert(isscalar(textValue) && ~ismissing(textValue), ...
    'compute_run09_sha256_text:BadInput', ...
    'The value to hash must be one nonmissing string scalar.');
bytes = unicode2native(char(textValue), 'UTF-8');
md = javaMethod('getInstance', 'java.security.MessageDigest', 'SHA-256');
md.update(typecast(uint8(bytes(:)), 'int8'));
digest = typecast(md.digest(), 'uint8');
value = lower(string(reshape(dec2hex(digest, 2).', 1, [])));
end
