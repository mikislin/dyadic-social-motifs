function base = local_base_row(T, rowIdx)
base = struct();
base.raw_index = T.raw_index(rowIdx);
base.session_id = string(T.session_id(rowIdx));
base.condition = string(T.condition(rowIdx));
base.condition_id = local_table_string(T, 'condition_id', rowIdx, string(T.condition(rowIdx)));
base.condition_label = local_table_string(T, 'condition_label', rowIdx, string(T.condition(rowIdx)));
base.arena = string(T.arena(rowIdx));
base.arena_condition = local_table_string(T, 'arena_condition', rowIdx, "");
base.arena_condition_label = local_table_string(T, 'arena_condition_label', rowIdx, "");
base.analysis_group_id = local_table_string(T, 'analysis_group_id', rowIdx, "");
base.analysis_group_label = local_table_string(T, 'analysis_group_label', rowIdx, "");
base.plot_group_label = local_table_string(T, 'plot_group_label', rowIdx, "");
base.analysis_class = string(T.analysis_class(rowIdx));
end
