function summary = local_assign_audit_summary(summary, rowIdx, audit)
fields = fieldnames(audit);
for i = 1:numel(fields)
    if ismember(fields{i}, summary.Properties.VariableNames)
        summary.(fields{i})(rowIdx) = audit.(fields{i});
    end
end
end