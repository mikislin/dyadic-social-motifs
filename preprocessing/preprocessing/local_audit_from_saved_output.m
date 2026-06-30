function audit = local_audit_from_saved_output(S)
if isfield(S, 'preprocess_qc_audit')
    audit = S.preprocess_qc_audit;
    if ~isfield(audit, 'badframe_fraction') && isfield(S, 'badframe_fraction')
        audit.badframe_fraction = S.badframe_fraction;
    end
elseif isfield(S, 'sessionPreproc')
    audit = local_preprocessing_audit(S.sessionPreproc);
elseif isfield(S, 'badframe_fraction')
    audit = local_empty_audit();
    audit.badframe_fraction = S.badframe_fraction;
else
    audit = local_empty_audit();
end
end