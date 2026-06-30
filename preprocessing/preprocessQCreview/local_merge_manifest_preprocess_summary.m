function T = local_merge_manifest_preprocess_summary(M, S)
T = M;
T.preprocess_status = repmat("missing", height(T), 1);
T.preprocess_runtime_sec = nan(height(T), 1);
T.badframe_fraction = nan(height(T), 1);
auditVars = local_preprocess_audit_vars();
for i = 1:numel(auditVars)
    T.(char(auditVars(i))) = nan(height(T), 1);
end
T.preprocess_output_file = strings(height(T), 1);
T.preprocess_error_message = strings(height(T), 1);

[tf, loc] = ismember(T.raw_index, S.raw_index);
T.preprocess_status(tf) = S.status(loc(tf));
T.preprocess_runtime_sec(tf) = S.runtime_sec(loc(tf));
T.badframe_fraction(tf) = S.badframe_fraction(loc(tf));
for i = 1:numel(auditVars)
    v = char(auditVars(i));
    if ismember(v, S.Properties.VariableNames)
        T.(v)(tf) = S.(v)(loc(tf));
    end
end
T.preprocess_output_file(tf) = S.output_file(loc(tf));
if ismember('error_message', S.Properties.VariableNames)
    T.preprocess_error_message(tf) = S.error_message(loc(tf));
end

T.preprocess_success = T.preprocess_status == "success" | ...
    T.preprocess_status == "skipped_existing";
end