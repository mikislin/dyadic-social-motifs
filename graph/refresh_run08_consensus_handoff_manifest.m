function Manifest = refresh_run08_consensus_handoff_manifest(outRoot)
%REFRESH_RUN08_CONSENSUS_HANDOFF_MANIFEST Rehash an existing frozen handoff.
%
% This audit-only repair does not rebuild graph edges, node eligibility, PCA,
% neighborhoods, consensus support, or topology. It streams the exact bytes
% of every artifact already declared by the handoff, verifies existing CSV
% row counts, records byte-accurate SHA-256 values, and atomically replaces
% the handoff manifest with a new deterministic freeze ID.

outRoot = string(outRoot);
manifestPath = fullfile(outRoot, 'run08_to_run09_handoff_manifest.csv');
assert(isfile(manifestPath), ...
    'refresh_run08_consensus_handoff_manifest:MissingManifest', ...
    'Missing run_08 handoff manifest: %s', manifestPath);

Manifest = readtable(manifestPath, 'TextType', 'string');
requiredColumns = ["artifact_id","artifact_path","sha256","file_bytes", ...
    "n_rows","required_for_run09","artifact_status","consensus_freeze_id", ...
    "artifact_role","permitted_run09_use"];
assert(all(ismember(requiredColumns, string(Manifest.Properties.VariableNames))), ...
    'refresh_run08_consensus_handoff_manifest:BadManifestSchema', ...
    'The existing handoff manifest is missing required columns.');

sourcePaths = [
    string(mfilename('fullpath')) + ".m"
    fullfile(fileparts(mfilename('fullpath')), 'compute_file_sha256.m')
    ];
for pathText = sourcePaths'
    [~, stem, ext] = fileparts(pathText);
    artifactId = string(stem) + string(ext);
    sourceMask = Manifest.artifact_id == artifactId;
    if any(sourceMask)
        Manifest.artifact_path(sourceMask) = pathText;
        Manifest.artifact_role(sourceMask) = "source_code_or_config_provenance";
        Manifest.permitted_run09_use(sourceMask) = ...
            "audit_only_exact_implementation_trace";
    else
        one = Manifest(1, :);
        one.artifact_id = artifactId;
        one.artifact_path = pathText;
        one.sha256 = "";
        one.file_bytes = 0;
        one.n_rows = NaN;
        one.required_for_run09 = false;
        one.artifact_status = "pending_rehash";
        one.artifact_role = "source_code_or_config_provenance";
        one.permitted_run09_use = "audit_only_exact_implementation_trace";
        Manifest = [Manifest; one]; %#ok<AGROW>
    end
end

oldRows = double(Manifest.n_rows);
for i = 1:height(Manifest)
    pathText = Manifest.artifact_path(i);
    assert(isfile(pathText), ...
        'refresh_run08_consensus_handoff_manifest:MissingArtifact', ...
        'Declared handoff artifact is missing: %s', pathText);
    [digest, byteCount, lineFeedCount, lastByte] = compute_file_sha256(pathText);
    Manifest.sha256(i) = digest;
    Manifest.file_bytes(i) = byteCount;
    Manifest.artifact_status(i) = "present_hashed_byte_exact";

    if endsWith(lower(pathText), ".csv")
        physicalLineCount = lineFeedCount + double(byteCount > 0 && lastByte ~= 10);
        dataRows = max(physicalLineCount - 1, 0);
        if isfinite(oldRows(i))
            assert(dataRows == oldRows(i), ...
                'refresh_run08_consensus_handoff_manifest:RowCountChanged', ...
                'CSV row count changed for %s: manifest=%d, observed=%d.', ...
                pathText, oldRows(i), dataRows);
        end
        Manifest.n_rows(i) = dataRows;
    else
        Manifest.n_rows(i) = NaN;
    end
end

edgeMask = Manifest.artifact_id == "run08_to_run09_edge_list.csv";
configMask = Manifest.artifact_id == "graph_consensus_freeze_config.csv";
assert(nnz(edgeMask) == 1 && nnz(configMask) == 1, ...
    'refresh_run08_consensus_handoff_manifest:MissingFreezeInputs', ...
    'The handoff must contain exactly one frozen edge input and freeze config.');
edgeHash = Manifest.sha256(edgeMask);
configHash = Manifest.sha256(configMask);
freezeId = "run08_consensus_" + extractBefore(edgeHash, 13) + "_" + ...
    extractBefore(configHash, 9);
Manifest.consensus_freeze_id(:) = freezeId;
Manifest.hash_algorithm = repmat("SHA-256", height(Manifest), 1);
Manifest.hash_scope = repmat("exact_file_bytes", height(Manifest), 1);
Manifest.hash_implementation = repmat( ...
    "matlab_fread_uint8_to_java_message_digest", height(Manifest), 1);
Manifest.manifest_schema_version = repmat( ...
    "run08_handoff_manifest_v2_byte_exact", height(Manifest), 1);

requiredMask = logical(Manifest.required_for_run09);
assert(all(Manifest.artifact_status(requiredMask) == "present_hashed_byte_exact") && ...
    all(strlength(Manifest.sha256(requiredMask)) == 64), ...
    'refresh_run08_consensus_handoff_manifest:RequiredArtifactFailure', ...
    'One or more required run_09 artifacts failed byte-exact hashing.');

temporaryPath = manifestPath + ".tmp.csv";
writetable(Manifest, temporaryPath);
[moved, message] = movefile(temporaryPath, manifestPath, 'f');
assert(moved, 'refresh_run08_consensus_handoff_manifest:AtomicReplaceFailed', ...
    'Could not replace handoff manifest: %s', message);
end
