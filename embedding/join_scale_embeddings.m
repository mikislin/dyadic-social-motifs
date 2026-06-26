function chunkTable = join_scale_embeddings(baseChunkTable, allTables, nameGroups)
%JOIN_SCALE_EMBEDDINGS Join per-scale embedding columns onto chunk metadata.

chunkTable = baseChunkTable;
for s = 1:numel(allTables)
    Ts = allTables{s};
    if isempty(Ts)
        continue
    end
    keepNames = ["chunk_id"; string(nameGroups{s}(:))];
    keepNames = keepNames(ismember(keepNames, string(Ts.Properties.VariableNames)));
    TsSmall = Ts(:, cellstr(keepNames));
    chunkTable = outerjoin(chunkTable, TsSmall, ...
        'Keys', 'chunk_id', 'MergeKeys', true, 'Type', 'left');
end
chunkTable = sortrows(chunkTable, 'chunk_id');
end
