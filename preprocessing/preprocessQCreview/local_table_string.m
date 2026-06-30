function v = local_table_string(T, varName, rowIdx, defaultValue)
if ismember(varName, T.Properties.VariableNames)
    v = string(T.(varName)(rowIdx));
else
    v = string(defaultValue);
end
end