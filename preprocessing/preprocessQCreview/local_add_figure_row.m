function rows = local_add_figure_row(rows, files, description)
if isempty(files)
    return
end
for i = 1:numel(files)
    row = struct();
    row.figure_file = string(files(i));
    row.description = string(description);
    rows = [rows; row]; %#ok<AGROW>
end
end
