function value = compute_run09_membership_sha256(membership)
%COMPUTE_RUN09_MEMBERSHIP_SHA256 Cross-platform hash of canonical labels.

membership = uint32(canonicalize_run09_membership(membership));
bytes = zeros(4 * numel(membership), 1, 'uint8');
bytes(1:4:end) = uint8(bitshift(membership, -24));
bytes(2:4:end) = uint8(bitand(bitshift(membership, -16), 255));
bytes(3:4:end) = uint8(bitand(bitshift(membership, -8), 255));
bytes(4:4:end) = uint8(bitand(membership, 255));
md = javaMethod('getInstance', 'java.security.MessageDigest', 'SHA-256');
md.update(typecast(bytes, 'int8'));
digest = typecast(md.digest(), 'uint8');
value = lower(string(reshape(dec2hex(digest, 2).', 1, [])));
end
