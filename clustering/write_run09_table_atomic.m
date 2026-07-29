function write_run09_table_atomic(T, outputPath)
%WRITE_RUN09_TABLE_ATOMIC Write a table then atomically replace its target.
%
% Retained as a run_09 compatibility name. New stages should call the
% repository-wide write_table_atomic helper.

write_table_atomic(T, outputPath);
end
