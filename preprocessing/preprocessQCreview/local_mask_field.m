function mask = local_mask_field(qa, fieldName, nFrames, nNodes)
if isfield(qa, fieldName) && ~isempty(qa.(fieldName))
    mask = logical(full(qa.(fieldName)));
else
    mask = false(nFrames, nNodes);
end
if ~isequal(size(mask), [nFrames nNodes])
    mask = false(nFrames, nNodes);
end
end