function T = local_struct_to_table(rows)
if isempty(rows)
    T = table();
else
    T = struct2table(rows);
end
end
