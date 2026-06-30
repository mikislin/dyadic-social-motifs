function value = struct_field_or_nan(S, fieldName)
%STRUCT_FIELD_OR_NAN Read a scalar struct field, or NaN when absent.

if isstruct(S) && isfield(S, fieldName)
    value = S.(fieldName);
else
    value = NaN;
end
end
