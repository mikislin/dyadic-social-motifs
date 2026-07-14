function outputs = run_multiscale_chunking_and_scale_selection(repoRoot, opts)
%RUN_MULTISCALE_CHUNKING_AND_SCALE_SELECTION Build run-06 chunks and scales.
%
% Condition-blind contract
%   - Session inclusion uses run_05 success, motif-discovery QC pass, dyadic
%     status, and canonical feature availability.
%   - Anchor selection uses frame validity, transformed feature availability,
%     time, session identity, and predefined caps only.
%   - Transform stats and scale scoring use feature values only.
%   - Condition, cohort, arena, drug, genotype, and outcome labels are copied
%     only into audit/provenance tables and audit figures.

if nargin < 1 || strlength(string(repoRoot)) == 0
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
end
if nargin < 2 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'configPath') || strlength(string(opts.configPath)) == 0
    opts.configPath = fullfile(repoRoot, 'config', 'multiscale_chunking_config.csv');
end

repoRoot = string(repoRoot);
cd(repoRoot);
addpath(genpath(repoRoot));

params = load_multiscale_chunking_config(opts.configPath);
outRoot = resolve_repo_path(repoRoot, params.output_dir);
local_assert_output_dir_compatible(repoRoot, outRoot, params);
figDir = fullfile(outRoot, 'figures');
sourceDir = fullfile(outRoot, 'figure_sources');
logDir = fullfile(outRoot, 'logs');
local_ensure_dir(outRoot);
local_ensure_dir(figDir);
local_ensure_dir(sourceDir);
local_ensure_dir(logDir);

diary(fullfile(logDir, 'run_06_build_multiscale_chunks_latest.log'));
cleanup = onCleanup(@() diary('off')); %#ok<NASGU>

paths = local_output_paths(outRoot);
fprintf('run_06_build_multiscale_chunks_and_select_scales\n');
fprintf('Repo root: %s\n', repoRoot);
fprintf('Output root: %s\n', outRoot);
fprintf('Run mode: %s\n', params.run_mode);

writetable(params.config_table, paths.parameterAudit);

[eligibleSessions, featureDict] = local_load_run05_inputs(repoRoot, params);
[featureNames, canonicalMeta] = default_dyad_feature_metadata();
local_validate_feature_dictionary(featureDict, canonicalMeta, featureNames);

selectedSessionRows = local_select_run_sessions(eligibleSessions, params);
sessionTable = eligibleSessions(selectedSessionRows, :);
sessionTable.feature_row_index = (1:height(sessionTable))';
sessionTable.session_index = (1:height(sessionTable))';

scaleTable = local_make_scale_table(params);
writetable(scaleTable, paths.scaleBank);

fprintf('Eligible run_05 sessions: %d\n', height(eligibleSessions));
fprintf('Selected sessions in this run: %d\n', height(sessionTable));
fprintf('Scale bank size: %d\n', height(scaleTable));

[stats, channelMeta, statsAudit] = local_fit_transform_stats(repoRoot, sessionTable, params);
transformAudit = local_feature_transform_audit(channelMeta, canonicalMeta, stats, statsAudit);
writetable(transformAudit, paths.transformAudit);

[candidateAnchors, inputAudit, scaleSessionCoverage] = local_build_anchor_candidates(repoRoot, eligibleSessions, ...
    selectedSessionRows, sessionTable, scaleTable, params);
scaleAnchorCoverage = local_scale_anchor_coverage_summary(scaleSessionCoverage);
selectedAnchors = select_condition_blind_anchor_subset(candidateAnchors, ...
    local_anchor_cap(params), ...
    'minAnchorsPerSession', params.min_materialized_anchors_per_session);
selectedAnchors = local_add_anchor_provenance_role(selectedAnchors);

inputAudit = local_add_materialized_counts(inputAudit, selectedAnchors);
writetable(inputAudit, paths.inputSessionAudit);
writetable(scaleSessionCoverage, paths.scaleSessionAnchorCoverage);
writetable(scaleAnchorCoverage, paths.scaleAnchorCoverage);

[ChunkSet, anchorManifest, shardManifest] = local_build_scale_summary_chunkset(repoRoot, outRoot, ...
    sessionTable, selectedAnchors, scaleTable, stats, canonicalMeta, params);
writetable(anchorManifest, paths.anchorManifest);
writetable(shardManifest, paths.scaleSummaryShardManifest);

chunkReport = validate_multiscale_chunk_bank(ChunkSet, ...
    'makePlots', false, ...
    'verbose', true);

chunkBankSummary = local_chunk_bank_summary(ChunkSet, chunkReport, candidateAnchors, scaleAnchorCoverage);
validitySummary = local_chunk_validity_summary(anchorManifest);
sessionScaleCounts = local_session_scale_counts(anchorManifest);
arenaQcCounts = local_arena_qc_counts(anchorManifest);
writetable(chunkBankSummary, paths.chunkBankSummary);
writetable(validitySummary, paths.validitySummary);
writetable(sessionScaleCounts, fullfile(sourceDir, 'chunk_valid_counts_by_session_scale.csv'));
writetable(arenaQcCounts, fullfile(sourceDir, 'chunk_valid_counts_by_arena_qc.csv'));

scoreFeatures = string(canonicalMeta.Name(canonicalMeta.ClusteringCandidate == 1));
ScaleScore = score_multiscale_chunk_bank(ChunkSet, ...
    'nPCs', params.score_n_pcs, ...
    'varianceTargetPct', params.score_variance_target_pct, ...
    'nClusters', params.score_n_clusters, ...
    'maxChunksPerScale', params.score_max_chunks_per_scale, ...
    'featureNames', scoreFeatures, ...
    'representationMode', params.score_representation_mode, ...
    'summaryTemporalBins', params.summary_temporal_bins, ...
    'summaryDctCoeffs', params.summary_dct_coeffs, ...
    'summaryIncludeDct', params.summary_include_dct, ...
    'summaryUseScaledFeatures', params.summary_use_scaled_features, ...
    'rngSeed', params.rng_seed, ...
    'microMaxSec', params.scale_band_micro_max_sec, ...
    'motifMaxSec', params.scale_band_motif_max_sec, ...
    'verbose', true);
ChunkSet = local_drop_summary_matrices(ChunkSet);
ScaleScore.selectedTable = select_operational_timescales(ScaleScore, ...
    'nMicro', params.selected_n_micro, ...
    'nMotif', params.selected_n_motif, ...
    'nContext', params.selected_n_context, ...
    'minLogGap', params.selected_min_log_gap);
scaleStability = estimate_scale_selection_stability(ScaleScore, ...
    'nBootstraps', local_stability_bootstrap_count(params), ...
    'nMicro', params.selected_n_micro, ...
    'nMotif', params.selected_n_motif, ...
    'nContext', params.selected_n_context, ...
    'minLogGap', params.selected_min_log_gap, ...
    'stabilityThreshold', params.stability_selection_frequency_threshold, ...
    'rngSeed', params.rng_seed + 106);
scaleReport = validate_scale_usefulness_scores(ScaleScore, ChunkSet, ...
    'makePlots', false, ...
    'verbose', true);

scaleScores = ScaleScore.scaleTable;
scaleScores.no_condition_labels_used = true(height(scaleScores), 1);
scaleScores.scoring_feature_rule = repmat("all_34_canonical_clustering_candidates", height(scaleScores), 1);
scaleScores.scoring_representation_rule = repmat(params.score_representation_mode, height(scaleScores), 1);
embeddingDimensionAudit = local_embedding_dimension_audit(ScaleScore.embeddingDimensionAudit, scaleScores, params);
selectedScales = local_annotate_selected_scales(ScaleScore.selectedTable, scaleStability, params);
ScaleScore.selectedTable = selectedScales;
selectionAudit = local_scale_selection_audit(scaleScores, selectedScales, params, scaleStability);
pcaLoadingStability = estimate_pca_loading_stability(ScaleScore, ...
    'nSplits', local_pca_loading_stability_split_count(params), ...
    'nPCs', params.score_n_pcs, ...
    'stabilityThreshold', params.pca_loading_stability_threshold, ...
    'rngSeed', params.rng_seed + 607);
arenaSensitivity = audit_arena_embedding_sensitivity(ScaleScore, canonicalMeta, ...
    'nEmbeddingPCs', params.arena_sensitivity_n_embedding_pcs, ...
    'topNFeatures', params.arena_sensitivity_top_n_features);
primaryScales = local_primary_operational_scales(selectedScales, pcaLoadingStability, params);
[primaryAnchorManifest, primaryChunkBankSummary] = build_primary_scale_specific_anchor_manifest( ...
    repoRoot, sessionTable, primaryScales, stats, ...
    'strideSec', params.stride_sec, ...
    'anchorMode', lower(params.anchor_mode), ...
    'minValidFrac', params.min_chunk_valid_frac, ...
    'minFeatureFiniteFrac', params.min_chunk_feature_finite_frac, ...
    'requireAnchorFrameValid', params.require_anchor_frame_valid, ...
    'maxAnchorsPerScale', local_scale_specific_anchor_cap(params), ...
    'minAnchorsPerSession', params.scale_specific_min_anchors_per_session);
primaryEventSummary = summarize_primary_chunk_events(repoRoot, sessionTable, primaryAnchorManifest, ...
    'contactThreshold', params.event_contact_threshold, ...
    'stateThreshold', params.event_state_threshold, ...
    'turnThresholdDeg', params.event_turn_threshold_deg);
writetable(scaleScores, paths.scaleScores);
writetable(selectedScales, paths.selectedScales);
writetable(primaryScales, paths.primaryScales);
writetable(embeddingDimensionAudit, paths.embeddingDimensionAudit);
writetable(scaleStability, paths.scaleSelectionStability);
writetable(selectionAudit, paths.scaleSelectionAudit);
writetable(pcaLoadingStability, paths.pcaLoadingStability);
writetable(arenaSensitivity, paths.arenaSensitivity);
writetable(primaryAnchorManifest, paths.primaryAnchorManifest);
writetable(primaryChunkBankSummary, paths.primaryChunkBankSummary);
writetable(primaryEventSummary, paths.primaryEventSummary);

coverageSource = chunkBankSummary;
exampleSource = local_example_trace_source(ChunkSet, selectedScales, params);
arenaShiftSource = local_arena_shift_after_transform(ChunkSet);
transformSource = transformAudit;
writetable(coverageSource, fullfile(sourceDir, 'chunk_coverage_by_temporal_scale.csv'));
writetable(scaleAnchorCoverage, fullfile(sourceDir, 'scale_anchor_coverage_source.csv'));
writetable(exampleSource, fullfile(sourceDir, 'chunk_example_trace_source.csv'));
writetable(arenaShiftSource, fullfile(sourceDir, 'chunk_arena_shift_after_transform.csv'));
writetable(transformSource, fullfile(sourceDir, 'chunk_feature_transform_source.csv'));
writetable(primaryScales, fullfile(sourceDir, 'primary_operational_scales_source.csv'));
writetable(primaryChunkBankSummary, fullfile(sourceDir, 'primary_scale_specific_chunk_bank_summary.csv'));
writetable(primaryEventSummary, fullfile(sourceDir, 'primary_chunk_event_summary_source.csv'));
writetable(pcaLoadingStability, fullfile(sourceDir, 'pca_loading_stability_source.csv'));
writetable(arenaSensitivity, fullfile(sourceDir, 'arena_embedding_sensitivity_source.csv'));

figureManifest = make_run06_chunk_qc_figures(outRoot, params);
writetable(figureManifest, paths.figureManifest);

chunkMatPath = fullfile(outRoot, char(params.chunk_mat_file));
chunkMatSaved = false;
if params.save_chunk_mat_artifact
    tmpChunkMatPath = string(chunkMatPath) + ".tmp";
    if isfile(tmpChunkMatPath)
        delete(tmpChunkMatPath);
    end
    save(tmpChunkMatPath, 'ChunkSet', 'ScaleScore', 'chunkReport', 'scaleReport', ...
        'scaleStability', 'embeddingDimensionAudit', 'pcaLoadingStability', ...
        'arenaSensitivity', 'primaryScales', 'primaryChunkBankSummary', ...
        'params', 'scaleTable', 'statsAudit', '-v7.3');
    movefile(tmpChunkMatPath, chunkMatPath, 'f');
    chunkMatSaved = true;
end

fprintf('Wrote parameter audit: %s\n', paths.parameterAudit);
fprintf('Wrote input session audit: %s\n', paths.inputSessionAudit);
fprintf('Wrote scale-session anchor coverage audit: %s\n', paths.scaleSessionAnchorCoverage);
fprintf('Wrote scale anchor coverage audit: %s\n', paths.scaleAnchorCoverage);
fprintf('Wrote transform audit: %s\n', paths.transformAudit);
fprintf('Wrote anchor manifest: %s\n', paths.anchorManifest);
fprintf('Wrote scale summary shard manifest: %s\n', paths.scaleSummaryShardManifest);
fprintf('Wrote chunk bank summary: %s\n', paths.chunkBankSummary);
fprintf('Wrote scale scores: %s\n', paths.scaleScores);
fprintf('Wrote selected scales: %s\n', paths.selectedScales);
fprintf('Wrote primary operational scales: %s\n', paths.primaryScales);
fprintf('Wrote embedding dimension audit: %s\n', paths.embeddingDimensionAudit);
fprintf('Wrote scale selection stability: %s\n', paths.scaleSelectionStability);
fprintf('Wrote PCA loading stability: %s\n', paths.pcaLoadingStability);
fprintf('Wrote arena sensitivity audit: %s\n', paths.arenaSensitivity);
fprintf('Wrote primary scale-specific anchor manifest: %s\n', paths.primaryAnchorManifest);
fprintf('Wrote primary chunk event summary: %s\n', paths.primaryEventSummary);
fprintf('Wrote figure manifest: %s\n', paths.figureManifest);
if chunkMatSaved
    fprintf('Wrote ignored MAT artifact: %s\n', chunkMatPath);
else
    fprintf('Skipped ignored MAT artifact by config: %s\n', chunkMatPath);
end

outputs = struct();
outputs.output_root = string(outRoot);
outputs.chunk_mat_path = string(chunkMatPath);
outputs.chunk_mat_saved = chunkMatSaved;
outputs.parameter_audit_path = string(paths.parameterAudit);
outputs.input_session_audit_path = string(paths.inputSessionAudit);
outputs.scale_session_anchor_coverage_path = string(paths.scaleSessionAnchorCoverage);
outputs.scale_anchor_coverage_path = string(paths.scaleAnchorCoverage);
outputs.feature_transform_audit_path = string(paths.transformAudit);
outputs.anchor_manifest_path = string(paths.anchorManifest);
outputs.scale_summary_shard_manifest_path = string(paths.scaleSummaryShardManifest);
outputs.chunk_bank_summary_path = string(paths.chunkBankSummary);
outputs.chunk_validity_summary_path = string(paths.validitySummary);
outputs.scale_scores_path = string(paths.scaleScores);
outputs.selected_scales_path = string(paths.selectedScales);
outputs.primary_scales_path = string(paths.primaryScales);
outputs.embedding_dimension_audit_path = string(paths.embeddingDimensionAudit);
outputs.scale_selection_stability_path = string(paths.scaleSelectionStability);
outputs.scale_selection_audit_path = string(paths.scaleSelectionAudit);
outputs.pca_loading_stability_path = string(paths.pcaLoadingStability);
outputs.arena_sensitivity_path = string(paths.arenaSensitivity);
outputs.primary_anchor_manifest_path = string(paths.primaryAnchorManifest);
outputs.primary_chunk_bank_summary_path = string(paths.primaryChunkBankSummary);
outputs.primary_event_summary_path = string(paths.primaryEventSummary);
outputs.figure_manifest_path = string(paths.figureManifest);
outputs.n_eligible_sessions = height(eligibleSessions);
outputs.n_selected_sessions = height(sessionTable);
outputs.n_candidate_anchors = height(candidateAnchors);
outputs.n_materialized_anchors = height(selectedAnchors);
outputs.n_chunk_rows = height(anchorManifest);
outputs.n_scales = height(scaleTable);
outputs.n_selected_scales = height(selectedScales);
outputs.n_primary_scales = height(primaryScales);
end

function paths = local_output_paths(outRoot)
paths = struct();
paths.parameterAudit = fullfile(outRoot, 'chunk_parameter_audit.csv');
paths.inputSessionAudit = fullfile(outRoot, 'chunk_input_session_audit.csv');
paths.scaleSessionAnchorCoverage = fullfile(outRoot, 'scale_session_anchor_coverage_audit.csv');
paths.scaleAnchorCoverage = fullfile(outRoot, 'scale_anchor_coverage_audit.csv');
paths.transformAudit = fullfile(outRoot, 'chunk_feature_transform_audit.csv');
paths.anchorManifest = fullfile(outRoot, 'chunk_anchor_manifest.csv');
paths.scaleSummaryShardManifest = fullfile(outRoot, 'scale_summary_shard_manifest.csv');
paths.chunkBankSummary = fullfile(outRoot, 'chunk_bank_summary.csv');
paths.validitySummary = fullfile(outRoot, 'chunk_validity_summary.csv');
paths.scaleScores = fullfile(outRoot, 'scale_usefulness_scores.csv');
paths.selectedScales = fullfile(outRoot, 'selected_operational_scales.csv');
paths.primaryScales = fullfile(outRoot, 'primary_operational_scales.csv');
paths.embeddingDimensionAudit = fullfile(outRoot, 'embedding_dimension_audit.csv');
paths.scaleSelectionStability = fullfile(outRoot, 'scale_selection_stability.csv');
paths.scaleSelectionAudit = fullfile(outRoot, 'scale_selection_audit.csv');
paths.pcaLoadingStability = fullfile(outRoot, 'pca_loading_stability.csv');
paths.arenaSensitivity = fullfile(outRoot, 'arena_sensitivity_audit.csv');
paths.primaryAnchorManifest = fullfile(outRoot, 'primary_scale_specific_anchor_manifest.csv');
paths.primaryChunkBankSummary = fullfile(outRoot, 'primary_scale_specific_chunk_bank_summary.csv');
paths.primaryEventSummary = fullfile(outRoot, 'primary_chunk_event_summary_audit.csv');
paths.figureManifest = fullfile(outRoot, 'chunk_qc_figure_manifest.csv');
paths.scaleBank = fullfile(outRoot, 'scale_bank_configured.csv');
end

function [ChunkSet, anchorManifest, shardManifest] = local_build_scale_summary_chunkset(repoRoot, outRoot, ...
    sessionTable, selectedAnchors, scaleTable, stats, canonicalMeta, params)
assert(params.use_scale_summary_shards, ...
    'run_multiscale_chunking:ScaleSummaryShardsRequired', ...
    'run_06 local mode now requires use_scale_summary_shards=true to avoid dense all-scale tensors.');

shardDir = fullfile(outRoot, char(params.scale_summary_shard_dir));
local_ensure_dir(shardDir);

scoreFeatures = string(canonicalMeta.Name(canonicalMeta.ClusteringCandidate == 1));
nScale = height(scaleTable);
scaleCells = cell(nScale, 1);
anchorCells = cell(nScale, 1);
manifestRows = cell(nScale, 1);

for s = 1:nScale
    scaleRow = scaleTable(s, :);
    shardPath = local_scale_summary_shard_path(shardDir, scaleRow);
    [ScaleShard, shardStatus] = local_load_or_build_scale_summary_shard(repoRoot, sessionTable, ...
        selectedAnchors, scaleRow, stats, params, scoreFeatures, shardPath);

    scaleCells{s} = ScaleShard.scale;
    anchorCells{s} = ScaleShard.scale.meta;
    manifestRows{s} = local_scale_summary_shard_manifest_row(ScaleShard, shardPath, shardStatus);
    fprintf('Scale shard %02d/%02d %.4fs: %s\n', s, nScale, ScaleShard.chunk_sec, shardStatus);
end

Scale = vertcat(scaleCells{:});
anchorManifest = vertcat(anchorCells{:});
anchorManifest = sortrows(anchorManifest, {'scale_index', 'anchor_id'});
shardManifest = vertcat(manifestRows{:});

firstScale = Scale(1);
ChunkSet = struct();
ChunkSet.stats = stats;
ChunkSet.obsNames = firstScale.obsNames;
ChunkSet.channelMeta = firstScale.channelMeta;
ChunkSet.featureNames = firstScale.featureNames;
ChunkSet.featureMeta = firstScale.featureMeta;
ChunkSet.chunkSec = double(scaleTable.chunk_sec(:))';
ChunkSet.strideSec = median(diff(unique(selectedAnchors.anchor_time_s)), 'omitnan');
if ~isfinite(ChunkSet.strideSec)
    ChunkSet.strideSec = NaN;
end
ChunkSet.anchorMode = char(params.anchor_mode);
ChunkSet.scale = Scale;
ChunkSet.chunkTable = sortrows(anchorManifest, 'chunk_id');
ChunkSet.sessions = local_session_stubs_from_table(sessionTable);
ChunkSet.nSessions = height(sessionTable);
ChunkSet.nObs = numel(firstScale.obsNames);
ChunkSet.scaleBankMeta = scaleTable;
ChunkSet.anchorTable = selectedAnchors;
ChunkSet.materializedAnchorCount = height(selectedAnchors);
ChunkSet.manifestChunkRows = height(anchorManifest);
ChunkSet.scaleSummaryShardManifest = shardManifest;
ChunkSet.usesScaleSummaryShards = true;
end

function [ScaleShard, statusText] = local_load_or_build_scale_summary_shard(repoRoot, sessionTable, ...
    selectedAnchors, scaleRow, stats, params, scoreFeatures, shardPath)
statusText = "built";
summaryProfile = local_summary_profile_for_scale(scaleRow, params);
if params.reuse_scale_summary_shards && isfile(shardPath)
    S = load(shardPath, 'ScaleShard');
    if isfield(S, 'ScaleShard') && local_scale_summary_shard_compatible(S.ScaleShard, ...
            selectedAnchors, scaleRow, params, summaryProfile)
        ScaleShard = S.ScaleShard;
        ScaleShard.summary_profile = summaryProfile.profile;
        ScaleShard.summary_temporal_band = summaryProfile.temporal_band;
        statusText = "reused";
        return
    end
    statusText = "rebuilt_incompatible_existing_shard";
end

ScaleShard = build_scale_summary_shard_from_anchor_manifest(repoRoot, sessionTable, ...
    selectedAnchors, scaleRow, stats, ...
    'fps', params.fps, ...
    'anchorMode', lower(params.anchor_mode), ...
    'minValidFrac', params.min_chunk_valid_frac, ...
    'minFeatureFiniteFrac', params.min_chunk_feature_finite_frac, ...
    'featureNames', scoreFeatures, ...
    'summaryTemporalBins', summaryProfile.temporal_bins, ...
    'summaryDctCoeffs', summaryProfile.dct_coeffs, ...
    'summaryIncludeDct', params.summary_include_dct, ...
    'summaryUseScaledFeatures', params.summary_use_scaled_features);
ScaleShard.summary_profile = summaryProfile.profile;
ScaleShard.summary_temporal_band = summaryProfile.temporal_band;

tmpPath = string(shardPath) + ".tmp";
if isfile(tmpPath)
    delete(tmpPath);
end
save(tmpPath, 'ScaleShard', '-v7.3');
movefile(tmpPath, shardPath, 'f');
end

function tf = local_scale_summary_shard_compatible(ScaleShard, selectedAnchors, scaleRow, params, summaryProfile)
tf = isstruct(ScaleShard) && isfield(ScaleShard, 'shard_version') && ...
    string(ScaleShard.shard_version) == "run06_scale_summary_shard_v1" && ...
    isfield(ScaleShard, 'scale') && isfield(ScaleShard.scale, 'Xsummary') && ...
    isfield(ScaleShard.scale, 'meta');
if ~tf
    return
end
tf = tf && ScaleShard.n_anchor_rows == height(selectedAnchors);
tf = tf && abs(double(ScaleShard.chunk_sec) - double(scaleRow.chunk_sec(1))) <= 1e-10;
tf = tf && ScaleShard.chunk_frames == max(1, round(double(scaleRow.chunk_sec(1)) * params.fps));
tf = tf && ScaleShard.summary_temporal_bins == summaryProfile.temporal_bins;
tf = tf && ScaleShard.summary_dct_coeffs == summaryProfile.dct_coeffs;
tf = tf && logical(ScaleShard.summary_include_dct) == logical(params.summary_include_dct);
tf = tf && logical(ScaleShard.summary_use_scaled_features) == logical(params.summary_use_scaled_features);
tf = tf && (~isfield(ScaleShard, 'summary_profile') || ...
    string(ScaleShard.summary_profile) == summaryProfile.profile);
tf = tf && height(ScaleShard.scale.meta) == height(selectedAnchors);
tf = tf && size(ScaleShard.scale.Xsummary, 1) == height(selectedAnchors);
end

function row = local_scale_summary_shard_manifest_row(ScaleShard, shardPath, statusText)
row = table();
row.scale_index = double(ScaleShard.scale_index);
row.chunk_sec = double(ScaleShard.chunk_sec);
row.chunk_frames = double(ScaleShard.chunk_frames);
row.n_anchor_rows = double(ScaleShard.n_anchor_rows);
row.n_valid_chunks = nnz(ScaleShard.scale.meta.chunk_is_valid);
row.n_summary_dims = size(ScaleShard.scale.Xsummary, 2);
row.summary_profile = local_struct_string_field(ScaleShard, 'summary_profile', "legacy_global");
row.summary_temporal_band = local_struct_string_field(ScaleShard, 'summary_temporal_band', "");
row.summary_temporal_bins = double(ScaleShard.summary_temporal_bins);
row.summary_dct_coeffs = double(ScaleShard.summary_dct_coeffs);
row.summary_storage_class = string(class(ScaleShard.scale.Xsummary));
row.shard_status = string(statusText);
row.shard_file = string(shardPath);
if isfile(shardPath)
    info = dir(shardPath);
    row.shard_bytes = double(info.bytes);
else
    row.shard_bytes = NaN;
end
row.shard_version = string(ScaleShard.shard_version);
row.labels_used_for_shard = "none";
row.condition_used_for_shard = false;
row.arena_used_for_shard = false;
end

function shardPath = local_scale_summary_shard_path(shardDir, scaleRow)
scaleIndex = double(scaleRow.scale_index(1));
chunkSec = double(scaleRow.chunk_sec(1));
secText = strrep(sprintf('%0.4f', chunkSec), '.', 'p');
stem = sprintf('run06_scale_%03d_%ss_summary_shard.mat', scaleIndex, secText);
shardPath = fullfile(shardDir, stem);
end

function sessions = local_session_stubs_from_table(sessionTable)
sessions = cell(height(sessionTable), 1);
for i = 1:height(sessionTable)
    S = struct();
    S.sessionIndex = i;
    S.raw_index = local_table_double(sessionTable, 'raw_index', i, i);
    S.session_id = local_table_string(sessionTable, 'session_id', i, "");
    S.session_file = local_table_string(sessionTable, 'session_file', i, "");
    S.feature_output_file = local_table_string(sessionTable, 'feature_output_file', i, "");
    S.arena = local_table_string(sessionTable, 'arena', i, "");
    S.arena_label = local_table_string(sessionTable, 'arena_label', i, "");
    S.cohort_id = local_table_string(sessionTable, 'cohort_id', i, "");
    S.condition_id = local_table_string(sessionTable, 'condition_id', i, "");
    S.qc_set = local_table_string(sessionTable, 'qc_set', i, "");
    S.fps = local_table_double(sessionTable, 'fps', i, NaN);
    S.nFrames = local_table_double(sessionTable, 'n_frames', i, NaN);
    S.validFrameFractionFromSeq = local_table_double(sessionTable, 'valid_frame_fraction', i, NaN);
    sessions{i} = S;
end
end

function ChunkSet = local_drop_summary_matrices(ChunkSet)
if ~isfield(ChunkSet, 'scale')
    return
end
for s = 1:numel(ChunkSet.scale)
    if isfield(ChunkSet.scale(s), 'Xsummary')
        ChunkSet.scale(s).Xsummary = [];
    end
end
end

function [eligible, featureDict] = local_load_run05_inputs(repoRoot, params)
featureRoot = resolve_repo_path(repoRoot, params.feature_output_dir);
sessionAuditPath = fullfile(featureRoot, char(params.feature_session_audit_file));
featureDictPath = fullfile(featureRoot, char(params.feature_dictionary_file));
assert(isfile(sessionAuditPath), 'run_multiscale_chunking:MissingRun05SessionAudit', ...
    'Missing run_05 session audit: %s', sessionAuditPath);
assert(isfile(featureDictPath), 'run_multiscale_chunking:MissingFeatureDictionary', ...
    'Missing run_05 feature dictionary: %s', featureDictPath);

T = readtable(sessionAuditPath, 'TextType', 'string');
featureDict = readtable(featureDictPath, 'TextType', 'string');
required = ["status", "include_motif_discovery", "qc_pass_motif_discovery", ...
    "feature_output_file", "raw_index", "session_id", "n_features"];
missing = required(ismember(required, string(T.Properties.VariableNames)) == 0);
assert(isempty(missing), 'run_multiscale_chunking:BadRun05SessionAudit', ...
    'Run_05 session audit missing required columns: %s', strjoin(missing, ', '));

mask = T.status == "success" & T.include_motif_discovery == 1 & ...
    T.qc_pass_motif_discovery == 1 & T.n_features == height(featureDict);
eligible = sortrows(T(mask, :), 'raw_index');
assert(~isempty(eligible), 'run_multiscale_chunking:NoEligibleSessions', ...
    'No successful run_05 motif-discovery feature sessions were found.');
end

function local_assert_output_dir_compatible(repoRoot, outRoot, params)
productionRoot = resolve_repo_path(repoRoot, "derived/chunks_motif_discovery");
if params.run_mode == "smoke" && ~params.allow_smoke_production_output && ...
        strcmpi(char(outRoot), char(productionRoot))
    error('run_multiscale_chunking:SmokeProductionOutputBlocked', ...
        ['Smoke mode is configured to write to the production chunk output directory: %s\n' ...
        'Use RUN06_CHUNK_OUTPUT_DIR=derived/chunks_motif_discovery_smoke, switch RUN06_CHUNK_RUN_MODE=full, ' ...
        'or set RUN06_ALLOW_SMOKE_PRODUCTION_OUTPUT=true for an intentional audit run.'], char(outRoot));
end
end

function local_validate_feature_dictionary(featureDict, canonicalMeta, featureNames)
assert(height(featureDict) == numel(featureNames), ...
    'run_multiscale_chunking:FeatureCountMismatch', ...
    'Feature dictionary must contain the canonical feature count.');
assert(all(string(featureDict.Name) == string(featureNames(:))), ...
    'run_multiscale_chunking:FeatureNameMismatch', ...
    'Run_05 feature dictionary does not match default_dyad_feature_metadata order.');
assert(all(featureDict.ClusteringCandidate == 1), ...
    'run_multiscale_chunking:FeatureCandidateMismatch', ...
    'All run_06 features must be canonical clustering candidates.');
assert(all(string(featureDict.Unit) == string(canonicalMeta.Unit)), ...
    'run_multiscale_chunking:FeatureUnitMismatch', ...
    'Run_05 feature dictionary units differ from canonical metadata.');
end

function selectedRows = local_select_run_sessions(eligible, params)
if params.run_mode == "smoke"
    n = min(height(eligible), floor(params.smoke_max_sessions));
    selectedRows = local_even_pick((1:height(eligible))', n);
else
    n = height(eligible);
    selectedRows = (1:n)';
end
if isfinite(params.max_sessions)
    n = min(n, floor(params.max_sessions));
    selectedRows = selectedRows(1:n);
end
end

function scaleTable = local_make_scale_table(params)
if params.run_mode == "smoke"
    nScales = params.smoke_n_scales;
else
    nScales = params.n_scales;
end
if params.use_log_scale_bank
    chunkSec = make_logscale_chunk_bank( ...
        'minSec', params.min_scale_sec, ...
        'maxSec', params.max_scale_sec, ...
        'nScales', nScales);
else
    chunkSec = linspace(params.min_scale_sec, params.max_scale_sec, nScales);
end
chunkSec = chunkSec(:);
scale_index = (1:numel(chunkSec))';
temporal_band = strings(numel(chunkSec), 1);
hierarchical_role = strings(numel(chunkSec), 1);
for i = 1:numel(chunkSec)
    temporal_band(i) = local_scale_band_label(chunkSec(i), params);
    hierarchical_role(i) = temporal_band(i);
end
scaleTable = table(scale_index, chunkSec, round(chunkSec .* params.fps), ...
    temporal_band, hierarchical_role, ...
    'VariableNames', {'scale_index', 'chunk_sec', 'chunk_frames_at_config_fps', ...
    'temporal_band', 'hierarchical_role'});
end

function [stats, channelMeta, statsAudit] = local_fit_transform_stats(repoRoot, sessionTable, params)
Xstats = [];
channelMeta = table();
rows = table();
for i = 1:height(sessionTable)
    dyad = local_load_dyad_from_session(repoRoot, sessionTable, i);
    Seq = prepare_dyad_timeseries(dyad);
    local_assert_seq_fps(Seq, params);
    local_assert_canonical_dyad(dyad);
    if isempty(channelMeta)
        channelMeta = Seq.channelMeta;
    end
    ok = Seq.validMask(:) & mean(isfinite(Seq.X), 2) >= params.min_chunk_feature_finite_frac;
    idx = find(ok);
    idx = local_even_pick(idx, min(numel(idx), floor(params.max_scaling_frames_per_session)));
    Xstats = [Xstats; Seq.X(idx, :)]; %#ok<AGROW>

    one = table();
    one.feature_row_index = i;
    one.raw_index = sessionTable.raw_index(i);
    one.session_id = string(sessionTable.session_id(i));
    one.n_valid_frames_available = nnz(ok);
    one.n_frames_sampled_for_stats = numel(idx);
    rows = [rows; one]; %#ok<AGROW>
end
stats = fit_condition_blind_chunk_stats(Xstats, channelMeta);
statsAudit = rows;
stats.n_fit_rows = size(Xstats, 1);
end

function [candidateAnchors, inputAudit, scaleSessionCoverage] = local_build_anchor_candidates(repoRoot, eligibleSessions, selectedRows, sessionTable, scaleTable, params)
candidateAnchors = table();
scaleSessionCoverage = table();
inputAudit = local_initialize_input_audit(eligibleSessions, selectedRows);
scaleSec = scaleTable.chunk_sec(:);

for i = 1:height(sessionTable)
    dyad = local_load_dyad_from_session(repoRoot, sessionTable, i);
    Seq = prepare_dyad_timeseries(dyad);
    local_assert_seq_fps(Seq, params);
    A = find_condition_blind_chunk_anchors(Seq, scaleSec, ...
        'strideSec', params.stride_sec, ...
        'anchorMode', lower(params.anchor_mode), ...
        'minValidFrac', params.min_chunk_valid_frac, ...
        'minFeatureFiniteFrac', params.min_chunk_feature_finite_frac, ...
        'requireAnchorFrameValid', params.require_anchor_frame_valid);
    scaleSessionCoverage = [scaleSessionCoverage; ...
        local_scale_session_anchor_survey(Seq, scaleTable, A, sessionTable, selectedRows, i, params)]; %#ok<AGROW>
    if ~isempty(A)
        n = height(A);
        A.feature_row_index = repmat(i, n, 1);
        A.input_audit_row = repmat(selectedRows(i), n, 1);
        A.session_index = repmat(i, n, 1);
        A.raw_index = repmat(sessionTable.raw_index(i), n, 1);
        A.session_id = repmat(string(sessionTable.session_id(i)), n, 1);
        A.session_file = repmat(local_table_string(sessionTable, 'session_file', i, ""), n, 1);
        A.arena = repmat(local_table_string(sessionTable, 'arena', i, ""), n, 1);
        A.arena_label = repmat(local_table_string(sessionTable, 'arena_label', i, ""), n, 1);
        A.cohort_id = repmat(local_table_string(sessionTable, 'cohort_id', i, ""), n, 1);
        A.cohort_label = repmat(local_table_string(sessionTable, 'cohort_label', i, ""), n, 1);
        A.condition_id = repmat(local_table_string(sessionTable, 'condition_id', i, ""), n, 1);
        A.condition_label = repmat(local_table_string(sessionTable, 'condition_label', i, ""), n, 1);
        A.condition_role = repmat("provenance_only_not_used_for_anchor_selection", n, 1);
        A.qc_set = repmat(local_table_string(sessionTable, 'qc_set', i, ""), n, 1);
        A.qc_pass_motif_discovery = repmat(sessionTable.qc_pass_motif_discovery(i), n, 1);
        A.qc_warn_motif_discovery = repmat(sessionTable.qc_warn_motif_discovery(i), n, 1);
        A.feature_output_file = repmat(string(sessionTable.feature_output_file(i)), n, 1);
        A.valid_frame_fraction = repmat(sessionTable.valid_frame_fraction(i), n, 1);
        A.badframe_mask_fraction = repmat(sessionTable.badframe_mask_fraction(i), n, 1);
        A.session_anchor_candidate_rank = (1:n)';
        candidateAnchors = [candidateAnchors; A]; %#ok<AGROW>
    end

    auditRow = selectedRows(i);
    inputAudit.anchor_candidate_count(auditRow) = height(A);
    inputAudit.n_frames_loaded(auditRow) = numel(Seq.time);
    inputAudit.fps_loaded(auditRow) = Seq.fps;
    inputAudit.valid_frame_fraction_loaded(auditRow) = mean(Seq.validMask);
    inputAudit.anchor_candidate_status(auditRow) = "success";
end

candidateAnchors = sortrows(candidateAnchors, {'raw_index', 'session_index', 'anchor_time_s'});
candidateAnchors.global_anchor_candidate_rank = (1:height(candidateAnchors))';
end

function T = local_scale_session_anchor_survey(Seq, scaleTable, commonAnchors, sessionTable, selectedRows, rowIdx, params)
T = table();
for s = 1:height(scaleTable)
    A = find_condition_blind_chunk_anchors(Seq, scaleTable.chunk_sec(s), ...
        'strideSec', params.stride_sec, ...
        'anchorMode', lower(params.anchor_mode), ...
        'minValidFrac', params.min_chunk_valid_frac, ...
        'minFeatureFiniteFrac', params.min_chunk_feature_finite_frac, ...
        'requireAnchorFrameValid', params.require_anchor_frame_valid);

    one = table();
    one.scale_index = scaleTable.scale_index(s);
    one.chunk_sec = scaleTable.chunk_sec(s);
    one.chunk_frames = scaleTable.chunk_frames_at_config_fps(s);
    one.temporal_band = string(scaleTable.temporal_band(s));
    one.feature_row_index = rowIdx;
    one.input_audit_row = selectedRows(rowIdx);
    one.session_index = rowIdx;
    one.raw_index = sessionTable.raw_index(rowIdx);
    one.session_id = string(sessionTable.session_id(rowIdx));
    one.arena = local_table_string(sessionTable, 'arena', rowIdx, "");
    one.arena_label = local_table_string(sessionTable, 'arena_label', rowIdx, "");
    one.cohort_id = local_table_string(sessionTable, 'cohort_id', rowIdx, "");
    one.condition_id = local_table_string(sessionTable, 'condition_id', rowIdx, "");
    one.qc_set = local_table_string(sessionTable, 'qc_set', rowIdx, "");
    one.n_frames_loaded = numel(Seq.time);
    one.fps_loaded = Seq.fps;
    one.valid_frame_fraction_loaded = mean(Seq.validMask);
    one.n_scale_specific_candidate_anchors = height(A);
    one.n_common_all_scale_candidate_anchors = height(commonAnchors);
    one.common_anchor_retention_fraction = local_fraction(height(commonAnchors), height(A));
    one.first_scale_specific_anchor_time_s = local_anchor_time_bound(A, "first");
    one.last_scale_specific_anchor_time_s = local_anchor_time_bound(A, "last");
    one.scale_specific_anchor_rule = "single_scale_frame_mask_feature_availability_fixed_stride_no_condition_labels";
    one.common_anchor_rule = "all_configured_scales_frame_mask_feature_availability_fixed_stride_no_condition_labels";
    one.labels_used_for_scale_survey = "none";
    one.arena_used_for_scale_survey = false;
    one.condition_used_for_scale_survey = false;
    T = [T; one]; %#ok<AGROW>
end
end

function value = local_anchor_time_bound(A, whichBound)
if isempty(A)
    value = NaN;
elseif whichBound == "first"
    value = min(A.anchor_time_s);
else
    value = max(A.anchor_time_s);
end
end

function value = local_fraction(num, den)
if den <= 0
    value = NaN;
else
    value = double(num) ./ double(den);
end
end

function T = local_scale_anchor_coverage_summary(scaleSessionCoverage)
T = table();
if isempty(scaleSessionCoverage)
    return
end
scales = unique(scaleSessionCoverage(:, {'scale_index', 'chunk_sec', 'chunk_frames', 'temporal_band'}), 'rows');
for i = 1:height(scales)
    idx = scaleSessionCoverage.scale_index == scales.scale_index(i);
    S = scaleSessionCoverage(idx, :);
    one = scales(i, :);
    one.n_sessions_surveyed = height(S);
    one.n_sessions_with_scale_specific_anchors = nnz(S.n_scale_specific_candidate_anchors > 0);
    one.n_sessions_with_common_all_scale_anchors = nnz(S.n_common_all_scale_candidate_anchors > 0);
    one.n_scale_specific_candidate_anchors = sum(S.n_scale_specific_candidate_anchors);
    one.n_common_all_scale_candidate_anchors = sum(S.n_common_all_scale_candidate_anchors);
    one.common_anchor_retention_fraction = local_fraction( ...
        one.n_common_all_scale_candidate_anchors, one.n_scale_specific_candidate_anchors);
    one.mean_session_common_anchor_retention_fraction = mean(S.common_anchor_retention_fraction, 'omitnan');
    one.min_session_common_anchor_retention_fraction = min(S.common_anchor_retention_fraction, [], 'omitnan');
    one.labels_used_for_scale_survey = "none";
    one.arena_used_for_scale_survey = false;
    one.condition_used_for_scale_survey = false;
    T = [T; one]; %#ok<AGROW>
end
end

function inputAudit = local_initialize_input_audit(eligible, selectedRows)
inputAudit = eligible;
inputAudit.selected_for_run06 = false(height(inputAudit), 1);
inputAudit.selected_for_run06(selectedRows) = true;
inputAudit.run06_selection_role = repmat("not_selected_in_this_run_mode", height(inputAudit), 1);
inputAudit.run06_selection_role(selectedRows) = "selected_condition_blind_even_raw_index";
inputAudit.anchor_candidate_count = zeros(height(inputAudit), 1);
inputAudit.materialized_anchor_count = zeros(height(inputAudit), 1);
inputAudit.n_frames_loaded = nan(height(inputAudit), 1);
inputAudit.fps_loaded = nan(height(inputAudit), 1);
inputAudit.valid_frame_fraction_loaded = nan(height(inputAudit), 1);
inputAudit.anchor_candidate_status = repmat("not_run", height(inputAudit), 1);
inputAudit.anchor_selection_rule = repmat( ...
    "frame_mask_feature_availability_fixed_stride_time_even_no_condition_labels", ...
    height(inputAudit), 1);
end

function inputAudit = local_add_materialized_counts(inputAudit, selectedAnchors)
if isempty(selectedAnchors)
    return
end
rows = unique(selectedAnchors.input_audit_row(:))';
for r = rows
    inputAudit.materialized_anchor_count(r) = nnz(selectedAnchors.input_audit_row == r);
end
end

function selectedAnchors = local_add_anchor_provenance_role(selectedAnchors)
if isempty(selectedAnchors)
    return
end
selectedAnchors.provenance_label_policy = repmat( ...
    "condition_cohort_arena_labels_carried_for_audit_only", ...
    height(selectedAnchors), 1);
end

function audit = local_feature_transform_audit(channelMeta, featureMeta, stats, statsAudit)
audit = channelMeta;
audit.channel_index = (1:height(audit))';
audit = movevars(audit, 'channel_index', 'Before', 1);
audit.robust_center = stats.center(:);
audit.robust_scale = stats.scale(:);
audit.impute_value = stats.impute(:);
audit.n_fit_rows_total = repmat(stats.n_fit_rows, height(audit), 1);
audit.stats_fit_scope = repmat(string(stats.fit_scope), height(audit), 1);
audit.condition_blind_transform = true(height(audit), 1);
audit.transform_source = repmat("features/default_dyad_feature_metadata.m plus prepare_dyad_timeseries", height(audit), 1);
audit.n_sessions_contributing = repmat(height(statsAudit), height(audit), 1);
audit.n_frames_sampled_for_stats = repmat(sum(statsAudit.n_frames_sampled_for_stats), height(audit), 1);

family = strings(height(audit), 1);
unit = strings(height(audit), 1);
isCircular = false(height(audit), 1);
isBoolean = false(height(audit), 1);
clusteringCandidate = false(height(audit), 1);
for i = 1:height(audit)
    row = find(featureMeta.Name == string(audit.BaseFeature(i)), 1);
    if ~isempty(row)
        family(i) = featureMeta.Family(row);
        unit(i) = featureMeta.Unit(row);
        isCircular(i) = featureMeta.IsCircular(row);
        isBoolean(i) = featureMeta.IsBoolean(row);
        clusteringCandidate(i) = featureMeta.ClusteringCandidate(row);
    end
end
audit.Family = family;
audit.Unit = unit;
audit.IsCircular = isCircular;
audit.IsBoolean = isBoolean;
audit.ClusteringCandidate = clusteringCandidate;
end

function summary = local_chunk_bank_summary(ChunkSet, Report, candidateAnchors, scaleAnchorCoverage)
summary = Report.scaleSummary;
summary.n_materialized_anchors = repmat(ChunkSet.materializedAnchorCount, height(summary), 1);
summary.n_candidate_anchors_before_cap = repmat(height(candidateAnchors), height(summary), 1);
summary.n_scale_specific_candidate_anchors = nan(height(summary), 1);
summary.n_common_all_scale_candidate_anchors = nan(height(summary), 1);
summary.n_sessions_with_scale_specific_anchors = nan(height(summary), 1);
summary.n_sessions_with_common_all_scale_anchors = nan(height(summary), 1);
summary.common_anchor_retention_fraction = nan(height(summary), 1);
summary.n_valid_chunks = zeros(height(summary), 1);
summary.valid_chunk_fraction = nan(height(summary), 1);
summary.mean_feature_finite_frac = nan(height(summary), 1);
summary.n_sessions_with_valid_chunks = zeros(height(summary), 1);
for s = 1:numel(ChunkSet.scale)
    meta = ChunkSet.scale(s).meta;
    covRow = find(scaleAnchorCoverage.scale_index == summary.scale_index(s), 1);
    if ~isempty(covRow)
        summary.n_scale_specific_candidate_anchors(s) = scaleAnchorCoverage.n_scale_specific_candidate_anchors(covRow);
        summary.n_common_all_scale_candidate_anchors(s) = scaleAnchorCoverage.n_common_all_scale_candidate_anchors(covRow);
        summary.n_sessions_with_scale_specific_anchors(s) = scaleAnchorCoverage.n_sessions_with_scale_specific_anchors(covRow);
        summary.n_sessions_with_common_all_scale_anchors(s) = scaleAnchorCoverage.n_sessions_with_common_all_scale_anchors(covRow);
        summary.common_anchor_retention_fraction(s) = scaleAnchorCoverage.common_anchor_retention_fraction(covRow);
    end
    summary.n_valid_chunks(s) = nnz(meta.chunk_is_valid);
    summary.valid_chunk_fraction(s) = mean(meta.chunk_is_valid);
    summary.mean_feature_finite_frac(s) = mean(meta.feature_finite_frac, 'omitnan');
    summary.n_sessions_with_valid_chunks(s) = numel(unique(meta.session_index(meta.chunk_is_valid)));
end
summary.condition_blind_anchor_rule = repmat("fixed_stride_mask_feature_time_even_no_condition_labels", height(summary), 1);
end

function T = local_chunk_validity_summary(anchorManifest)
scales = unique(anchorManifest.scale_index, 'stable')';
T = table();
for s = scales
    idx = anchorManifest.scale_index == s;
    M = anchorManifest(idx, :);
    one = table();
    one.scale_index = s;
    one.chunk_sec = M.chunk_sec(1);
    one.chunk_frames = M.chunk_frames(1);
    one.n_chunks = height(M);
    one.n_valid_chunks = nnz(M.chunk_is_valid);
    one.valid_chunk_fraction = mean(M.chunk_is_valid);
    one.mean_valid_frac = mean(M.valid_frac, 'omitnan');
    one.min_valid_frac = min(M.valid_frac);
    one.mean_feature_finite_frac = mean(M.feature_finite_frac, 'omitnan');
    one.min_feature_finite_frac = min(M.feature_finite_frac);
    T = [T; one]; %#ok<AGROW>
end
end

function T = local_session_scale_counts(anchorManifest)
groups = unique(anchorManifest(:, {'scale_index', 'chunk_sec', 'raw_index', ...
    'session_id', 'arena_label', 'qc_set'}), 'rows');
T = table();
for i = 1:height(groups)
    idx = anchorManifest.scale_index == groups.scale_index(i) & ...
        anchorManifest.raw_index == groups.raw_index(i);
    one = groups(i, :);
    one.n_chunks = nnz(idx);
    one.n_valid_chunks = nnz(anchorManifest.chunk_is_valid(idx));
    one.mean_valid_frac = mean(anchorManifest.valid_frac(idx), 'omitnan');
    T = [T; one]; %#ok<AGROW>
end
end

function T = local_arena_qc_counts(anchorManifest)
groups = unique(anchorManifest(:, {'scale_index', 'chunk_sec', 'arena_label', 'qc_set'}), 'rows');
T = table();
for i = 1:height(groups)
    idx = anchorManifest.scale_index == groups.scale_index(i) & ...
        anchorManifest.arena_label == groups.arena_label(i) & ...
        anchorManifest.qc_set == groups.qc_set(i);
    one = groups(i, :);
    one.n_chunks = nnz(idx);
    one.n_valid_chunks = nnz(anchorManifest.chunk_is_valid(idx));
    one.mean_valid_frac = mean(anchorManifest.valid_frac(idx), 'omitnan');
    T = [T; one]; %#ok<AGROW>
end
end

function selectedScales = local_annotate_selected_scales(selectedScales, scaleStability, params)
if isempty(selectedScales)
    return
end
selectedScales.bootstrap_selection_frequency = nan(height(selectedScales), 1);
selectedScales.bootstrap_median_role_rank = nan(height(selectedScales), 1);
selectedScales.passes_stability_threshold = false(height(selectedScales), 1);
if istable(scaleStability) && ~isempty(scaleStability)
    [tf, loc] = ismember(selectedScales.scale_index, scaleStability.scale_index);
    selectedScales.bootstrap_selection_frequency(tf) = scaleStability.selection_frequency(loc(tf));
    selectedScales.bootstrap_median_role_rank(tf) = scaleStability.median_role_rank(loc(tf));
    selectedScales.passes_stability_threshold(tf) = scaleStability.passes_stability_threshold(loc(tf));
end
selectedScales.dimension_guard_min_effective_dim = repmat(params.dimension_guard_min_effective_dim, height(selectedScales), 1);
selectedScales.dimension_guard_max_pc1_explained_pct = repmat(params.dimension_guard_max_pc1_explained_pct, height(selectedScales), 1);
selectedScales.dominant_pc_warning = selectedScales.effective_dim < params.dimension_guard_min_effective_dim | ...
    selectedScales.pc1_explained > params.dimension_guard_max_pc1_explained_pct;
selectedScales.passes_dimension_guard = ~selectedScales.dominant_pc_warning;
selectedScales.stable_and_dimension_supported = selectedScales.passes_stability_threshold & selectedScales.passes_dimension_guard;
selectedScales.selection_status_note = repmat("selected_by_predefined_score_then_audited_for_stability_and_dimension", height(selectedScales), 1);
selectedScales.labels_used_for_selection = repmat("none", height(selectedScales), 1);
end

function T = local_scale_selection_audit(scaleScores, selectedScales, params, scaleStability)
T = scaleScores;
T.selected_operational_scale = ismember(T.scale_index, selectedScales.scale_index);
T.selected_role = repmat("", height(T), 1);
T.selection_rank_within_role = nan(height(T), 1);
for i = 1:height(selectedScales)
    row = find(T.scale_index == selectedScales.scale_index(i), 1);
    if ~isempty(row)
        T.selected_role(row) = selectedScales.hierarchical_role(i);
        T.selection_rank_within_role(row) = selectedScales.rank_within_role(i);
    end
end
T.scale_selection_rule = repmat(params.scale_selection_rule, height(T), 1);
T.labels_used_for_selection = repmat("none", height(T), 1);
T.arena_used_for_selection = false(height(T), 1);
T.condition_used_for_selection = false(height(T), 1);
T.dimension_guard_min_effective_dim = repmat(params.dimension_guard_min_effective_dim, height(T), 1);
T.dimension_guard_max_pc1_explained_pct = repmat(params.dimension_guard_max_pc1_explained_pct, height(T), 1);
T.dominant_pc_warning = T.effective_dim < params.dimension_guard_min_effective_dim | ...
    T.pc1_explained > params.dimension_guard_max_pc1_explained_pct;
T.passes_dimension_guard = ~T.dominant_pc_warning;
if nargin >= 4 && istable(scaleStability) && ~isempty(scaleStability)
    [tf, loc] = ismember(T.scale_index, scaleStability.scale_index);
    T.bootstrap_selection_frequency = nan(height(T), 1);
    T.bootstrap_median_role_rank = nan(height(T), 1);
    T.bootstrap_passes_stability_threshold = false(height(T), 1);
    T.bootstrap_selection_frequency(tf) = scaleStability.selection_frequency(loc(tf));
    T.bootstrap_median_role_rank(tf) = scaleStability.median_role_rank(loc(tf));
    T.bootstrap_passes_stability_threshold(tf) = scaleStability.passes_stability_threshold(loc(tf));
end
T.stable_and_dimension_supported = T.selected_operational_scale & ...
    T.bootstrap_passes_stability_threshold & T.passes_dimension_guard;
end

function T = local_embedding_dimension_audit(dimensionAudit, scaleScores, params)
T = dimensionAudit;
if isempty(T)
    return
end
[tf, loc] = ismember(T.scale_index, scaleScores.scale_index);
T.initial_band = strings(height(T), 1);
T.n_pcs_target_pct_from_scale_scores = nan(height(T), 1);
T.cum_scoring_pcs_explained = nan(height(T), 1);
T.initial_band(tf) = scaleScores.initial_band(loc(tf));
T.n_pcs_target_pct_from_scale_scores(tf) = scaleScores.n_pcs_target_pct(loc(tf));
T.cum_scoring_pcs_explained(tf) = scaleScores.cum_embedding_pcs_explained(loc(tf));
T.summary_profile = strings(height(T), 1);
for i = 1:height(T)
    oneScale = table();
    oneScale.chunk_sec = T.chunk_sec(i);
    prof = local_summary_profile_for_scale(oneScale, params);
    T.summary_profile(i) = prof.profile;
end

minPcs = nan(height(T), 1);
maxPcs = nan(height(T), 1);
for i = 1:height(T)
    switch string(T.initial_band(i))
        case "micro"
            minPcs(i) = params.embedding_dim_micro_min_pcs;
            maxPcs(i) = params.embedding_dim_micro_max_pcs;
        case "motif"
            minPcs(i) = params.embedding_dim_motif_min_pcs;
            maxPcs(i) = params.embedding_dim_motif_max_pcs;
        otherwise
            minPcs(i) = params.embedding_dim_context_min_pcs;
            maxPcs(i) = params.embedding_dim_context_max_pcs;
    end
end
targetPcs = T.n_pcs_target_pct_from_scale_scores;
targetPcs(~isfinite(targetPcs)) = T.n_score_pcs_retained(~isfinite(targetPcs));
T.recommended_min_pcs_by_temporal_role = minPcs;
T.recommended_max_pcs_by_temporal_role = maxPcs;
T.recommended_pcs_for_next_embedding = min(max(targetPcs, minPcs), maxPcs);
T.recommendation_reaches_target_variance = targetPcs <= maxPcs;
T.embedding_dimension_rule = repmat( ...
    "condition_blind_multiresolution_pca_target_with_predefined_temporal_role_caps", ...
    height(T), 1);
T.labels_used_for_embedding_dimension = repmat("none", height(T), 1);
end

function T = local_primary_operational_scales(selectedScales, pcaLoadingStability, params)
T = selectedScales;
if isempty(T)
    return
end
T.primary_scale_rule = repmat(params.primary_scale_rule, height(T), 1);
T.pca_loading_median_subspace_similarity = nan(height(T), 1);
T.passes_pca_loading_stability_threshold = false(height(T), 1);
if istable(pcaLoadingStability) && ~isempty(pcaLoadingStability)
    [tf, loc] = ismember(T.scale_index, pcaLoadingStability.scale_index);
    T.pca_loading_median_subspace_similarity(tf) = pcaLoadingStability.median_subspace_similarity(loc(tf));
    T.passes_pca_loading_stability_threshold(tf) = pcaLoadingStability.passes_loading_stability_threshold(loc(tf));
end
T.primary_operational_scale = T.stable_and_dimension_supported;
T.primary_exclusion_reason = strings(height(T), 1);
for i = 1:height(T)
    if T.primary_operational_scale(i)
        T.primary_exclusion_reason(i) = "included_primary_scale";
    elseif ~T.passes_stability_threshold(i)
        T.primary_exclusion_reason(i) = "failed_scale_selection_bootstrap_stability";
    elseif ~T.passes_dimension_guard(i)
        T.primary_exclusion_reason(i) = "failed_pca_dimension_guard";
    else
        T.primary_exclusion_reason(i) = "not_primary_under_stable_and_dimension_supported_rule";
    end
end
T.labels_used_for_primary_scale_rule = repmat("none", height(T), 1);
T.arena_used_for_primary_scale_rule = false(height(T), 1);
T.condition_used_for_primary_scale_rule = false(height(T), 1);
T = T(T.primary_operational_scale, :);
end

function n = local_stability_bootstrap_count(params)
if params.run_mode == "smoke"
    n = params.smoke_stability_bootstraps;
else
    n = params.stability_bootstraps;
end
end

function n = local_pca_loading_stability_split_count(params)
if params.run_mode == "smoke"
    n = params.smoke_pca_loading_stability_splits;
else
    n = params.pca_loading_stability_splits;
end
end

function n = local_scale_specific_anchor_cap(params)
if params.run_mode == "smoke"
    n = params.smoke_scale_specific_max_anchors_per_scale;
else
    n = params.scale_specific_max_anchors_per_scale;
end
end

function T = local_example_trace_source(ChunkSet, selectedScales, params)
features = ["centroid_dist", "mutual_facing", "heading_diff_deg", "in_contact"];
T = table();
if isempty(selectedScales)
    return
end
for r = 1:height(selectedScales)
    [~, s] = min(abs([ChunkSet.scale.chunkSec] - selectedScales.chunk_sec(r)));
    Sc = ChunkSet.scale(s);
    nUse = min(height(Sc.meta), floor(params.example_anchor_count));
    idx = local_even_pick((1:height(Sc.meta))', nUse);
    for ii = 1:numel(idx)
        chunkIdx = idx(ii);
        for f = 1:numel(features)
            if isfield(Sc, 'Xraw') && ~isempty(Sc.Xraw)
                x = extract_chunk_feature_trace(ChunkSet, s, chunkIdx, features(f), ...
                    'invertLog1p', true);
                valid = Sc.valid(chunkIdx, :)';
            else
                [x, valid] = local_extract_chunk_feature_trace_from_meta(ChunkSet, Sc, chunkIdx, features(f), true);
            end
            t = ((0:numel(x)-1)' ./ ChunkSet.sessions{Sc.meta.feature_row_index(chunkIdx)}.fps);
            one = table();
            one.scale_index = repmat(s, numel(x), 1);
            one.chunk_sec = repmat(Sc.chunkSec, numel(x), 1);
            one.hierarchical_role = repmat(string(selectedScales.hierarchical_role(r)), numel(x), 1);
            one.example_rank = repmat(ii, numel(x), 1);
            one.anchor_id = repmat(Sc.meta.anchor_id(chunkIdx), numel(x), 1);
            one.raw_index = repmat(Sc.meta.raw_index(chunkIdx), numel(x), 1);
            one.session_id = repmat(string(Sc.meta.session_id(chunkIdx)), numel(x), 1);
            one.feature_name = repmat(features(f), numel(x), 1);
            one.time_within_chunk_s = t;
            one.value = x;
            one.valid_frame = valid;
            T = [T; one]; %#ok<AGROW>
        end
    end
end
end

function T = local_arena_shift_after_transform(ChunkSet)
T = table();
if isempty(ChunkSet.scale) || isempty(ChunkSet.scale(1).meta)
    return
end
Sc = ChunkSet.scale(1);
if isfield(Sc, 'centerScaled') && ~isempty(Sc.centerScaled)
    X = double(Sc.centerScaled);
else
    L = size(Sc.X, 2);
    if string(ChunkSet.anchorMode) == "past"
        centerIdx = L;
    else
        centerIdx = floor((L - 1) / 2) + 1;
    end
    X = squeeze(double(Sc.X(:, centerIdx, :)));
end
if isvector(X)
    X = X(:);
end
arena = string(Sc.meta.arena_label);
arenaLevels = unique(arena, 'stable');
if numel(arenaLevels) < 2
    return
end
for d = 1:size(X, 2)
    valsA = X(arena == arenaLevels(1), d);
    valsB = X(arena == arenaLevels(2), d);
    one = table();
    one.channel_index = d;
    one.obs_name = string(ChunkSet.channelMeta.ObsName(d));
    one.base_feature = string(ChunkSet.channelMeta.BaseFeature(d));
    one.channel_type = string(ChunkSet.channelMeta.ChannelType(d));
    one.arena_a = arenaLevels(1);
    one.arena_b = arenaLevels(2);
    one.n_a = nnz(isfinite(valsA));
    one.n_b = nnz(isfinite(valsB));
    one.median_a = median(valsA, 'omitnan');
    one.median_b = median(valsB, 'omitnan');
    pooled = iqr([valsA(isfinite(valsA)); valsB(isfinite(valsB))]);
    if ~(isfinite(pooled) && pooled > 0)
        pooled = 1;
    end
    one.robust_median_difference_iqr = (one.median_b - one.median_a) ./ pooled;
    one.audit_only_not_selection = true;
    T = [T; one]; %#ok<AGROW>
end
end

function [x, valid] = local_extract_chunk_feature_trace_from_meta(ChunkSet, Sc, chunkIdx, featureName, invertLog1p)
featureName = string(featureName);
meta = Sc.meta(chunkIdx, :);
sess = meta.feature_row_index(1);
featurePath = resolve_repo_path(string(pwd), string(meta.feature_output_file(1)));
S = load(featurePath, 'dyad', 'status');
assert(isfield(S, 'dyad'), 'run_multiscale_chunking:MissingTraceDyad', ...
    'Missing dyad in feature file: %s', featurePath);
if isfield(S, 'status')
    assert(string(S.status) == "success", ...
        'run_multiscale_chunking:BadTraceDyadStatus', ...
        'Feature file does not report success: %s', featurePath);
end
Seq = prepare_dyad_timeseries(S.dyad, 'stats', ChunkSet.stats);
idx = round(meta.start_frame(1)):round(meta.stop_frame(1));
idx = idx(idx >= 1 & idx <= numel(Seq.time));
valid = logical(Seq.validMask(idx));

fMeta = ChunkSet.featureMeta;
names = string(fMeta.Name);
fRow = find(names == featureName, 1);
assert(~isempty(fRow), 'run_multiscale_chunking:TraceFeatureNotFound', ...
    'Feature not found in ChunkSet.featureMeta: %s', featureName);

base = string(Seq.channelMeta.BaseFeature);
ctype = string(Seq.channelMeta.ChannelType);
isCirc = logical(fMeta.IsCircular(fRow));
isBool = logical(fMeta.IsBoolean(fRow));
hint = string(fMeta.TransformHint(fRow));

if isCirc
    idxCos = find(base == featureName & ctype == "circular_cos", 1);
    idxSin = find(base == featureName & ctype == "circular_sin", 1);
    assert(~isempty(idxCos) && ~isempty(idxSin), ...
        'run_multiscale_chunking:MissingCircularTraceChannels', ...
        'Missing circular cos/sin channels for %s.', featureName);
    x = atan2d(Seq.X(idx, idxSin), Seq.X(idx, idxCos));
else
    ch = find(base == featureName, 1);
    if isempty(ch)
        ch = find(string(Seq.obsNames) == featureName, 1);
    end
    assert(~isempty(ch), 'run_multiscale_chunking:MissingTraceChannel', ...
        'Missing channel for feature %s.', featureName);
    x = Seq.X(idx, ch);
    if invertLog1p && hint == "log1p"
        x = expm1(x);
    end
    if isBool
        x = double(x > 0.5);
    end
end

x = x(:);
valid = valid(:);
if numel(valid) ~= numel(x)
    valid = true(numel(x), 1);
end
if isnan(ChunkSet.sessions{sess}.fps)
    ChunkSet.sessions{sess}.fps = Seq.fps; %#ok<STRNU>
end
end

function dyad = local_load_dyad_from_session(repoRoot, sessionTable, rowIdx)
featurePath = resolve_repo_path(repoRoot, string(sessionTable.feature_output_file(rowIdx)));
S = load(featurePath, 'dyad', 'status');
assert(isfield(S, 'dyad'), 'run_multiscale_chunking:MissingDyad', ...
    'Missing dyad in feature file: %s', featurePath);
if isfield(S, 'status')
    assert(string(S.status) == "success", ...
        'run_multiscale_chunking:BadDyadStatus', ...
        'Feature file does not report success: %s', featurePath);
end
dyad = S.dyad;
end

function local_assert_canonical_dyad(dyad)
[featureNames, featureMeta] = default_dyad_feature_metadata();
assert(numel(dyad.featureNames) == numel(featureNames), ...
    'run_multiscale_chunking:DyadFeatureCountMismatch', ...
    'Dyad feature count differs from canonical feature metadata.');
assert(all(string(dyad.featureNames(:)) == string(featureNames(:))), ...
    'run_multiscale_chunking:DyadFeatureOrderMismatch', ...
    'Dyad feature names differ from canonical feature order.');
assert(all(dyad.featureMeta.ClusteringCandidate == featureMeta.ClusteringCandidate), ...
    'run_multiscale_chunking:DyadFeatureMetaMismatch', ...
    'Dyad feature metadata differs from canonical clustering-candidate flags.');
end

function local_assert_seq_fps(Seq, params)
assert(abs(double(Seq.fps) - double(params.fps)) <= 1e-10, ...
    'run_multiscale_chunking:FpsMismatch', ...
    'Loaded feature session fps %.15g does not match run_06 config fps %.15g.', ...
    double(Seq.fps), double(params.fps));
end

function n = local_anchor_cap(params)
if params.run_mode == "smoke"
    n = params.smoke_max_anchor_manifest_anchors;
else
    n = params.max_anchor_manifest_anchors;
end
end

function idx = local_even_pick(idx, nWant)
idx = idx(:);
nWant = min(nWant, numel(idx));
if nWant <= 0
    idx = zeros(0, 1);
elseif nWant < numel(idx)
    pos = unique(round(linspace(1, numel(idx), nWant)));
    while numel(pos) < nWant
        missing = setdiff(1:numel(idx), pos, 'stable');
        pos(end + 1) = missing(1); %#ok<AGROW>
    end
    idx = idx(sort(pos(:)));
end
end

function label = local_scale_band_label(sec, params)
if sec < params.scale_band_micro_max_sec
    label = "micro";
elseif sec < params.scale_band_motif_max_sec
    label = "motif";
else
    label = "context";
end
end

function profile = local_summary_profile_for_scale(scaleRow, params)
sec = double(scaleRow.chunk_sec(1));
band = local_scale_band_label(sec, params);
profile = struct();
profile.temporal_band = band;
if isfield(params, 'summary_band_adaptive') && logical(params.summary_band_adaptive)
    switch band
        case "micro"
            profile.temporal_bins = params.summary_micro_temporal_bins;
            profile.dct_coeffs = params.summary_micro_dct_coeffs;
        case "motif"
            profile.temporal_bins = params.summary_motif_temporal_bins;
            profile.dct_coeffs = params.summary_motif_dct_coeffs;
        otherwise
            profile.temporal_bins = params.summary_context_temporal_bins;
            profile.dct_coeffs = params.summary_context_dct_coeffs;
    end
    profile.profile = band + "_adaptive_bins" + string(profile.temporal_bins) + ...
        "_dct" + string(profile.dct_coeffs);
else
    profile.temporal_bins = params.summary_temporal_bins;
    profile.dct_coeffs = params.summary_dct_coeffs;
    profile.profile = "global_bins" + string(profile.temporal_bins) + ...
        "_dct" + string(profile.dct_coeffs);
end
end

function value = local_struct_string_field(S, fieldName, defaultValue)
if isfield(S, fieldName)
    value = string(S.(fieldName));
else
    value = string(defaultValue);
end
end

function value = local_table_string(T, varName, rowIdx, defaultValue)
if ismember(varName, T.Properties.VariableNames)
    value = string(T.(varName)(rowIdx));
else
    value = string(defaultValue);
end
end

function value = local_table_double(T, varName, rowIdx, defaultValue)
if ismember(varName, T.Properties.VariableNames)
    value = double(T.(varName)(rowIdx));
else
    value = defaultValue;
end
end

function local_ensure_dir(pathText)
if ~exist(pathText, 'dir')
    mkdir(pathText);
end
end
