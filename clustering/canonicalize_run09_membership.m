function membership = canonicalize_run09_membership(membership)
%CANONICALIZE_RUN09_MEMBERSHIP Relabel communities by first node occurrence.

membership = double(membership(:));
assert(all(isfinite(membership)) && all(membership >= 1) && ...
    all(membership == floor(membership)), ...
    'canonicalize_run09_membership:BadMembership', ...
    'Membership labels must be finite positive integers.');
labelMap = zeros(max(membership), 1);
canonical = zeros(size(membership));
nextLabel = 0;
for i = 1:numel(membership)
    old = membership(i);
    if labelMap(old) == 0
        nextLabel = nextLabel + 1;
        labelMap(old) = nextLabel;
    end
    canonical(i) = labelMap(old);
end
membership = canonical;
end
