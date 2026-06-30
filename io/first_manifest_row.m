function idx = first_manifest_row(manifestTable, mask, errorId, description)
%FIRST_MANIFEST_ROW Return the first manifest row matching a logical mask.

if nargin < 3 || strlength(string(errorId)) == 0
    errorId = 'first_manifest_row:MissingRow';
end
if nargin < 4
    description = "requested manifest subset";
end

mask = logical(mask);
assert(height(manifestTable) == numel(mask), ...
    'first_manifest_row:MaskSizeMismatch', ...
    'Mask length must match manifest table height.');

idx = find(mask, 1, 'first');
assert(~isempty(idx), errorId, ...
    'Could not find a manifest row for %s.', description);
end
