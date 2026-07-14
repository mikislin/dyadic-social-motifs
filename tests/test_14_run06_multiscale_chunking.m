function tests = test_14_run06_multiscale_chunking
tests = functiontests(localfunctions);
end

function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));
end

function testRun06ConfigLoads(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
params = load_multiscale_chunking_config(fullfile(repoRoot, 'config', 'multiscale_chunking_config.csv'));

verifyEqual(testCase, params.fps, 80);
verifyTrue(testCase, params.min_chunk_valid_frac > 0);
verifyTrue(testCase, params.min_chunk_feature_finite_frac > 0);
verifyEqual(testCase, params.anchor_selection_rule, ...
    "dyad.frameMask + transformed feature availability + fixed stride + deterministic time-even materialization; no condition cohort arena drug genotype or outcome labels");
verifyTrue(testCase, any(params.run_mode == ["smoke", "full"]));
verifyEqual(testCase, params.output_dir, "derived/chunks_motif_discovery_smoke");
verifyFalse(testCase, params.allow_smoke_production_output);
verifyFalse(testCase, params.save_chunk_mat_artifact);
verifyTrue(testCase, params.use_scale_summary_shards);
verifyTrue(testCase, params.reuse_scale_summary_shards);
verifyEqual(testCase, params.scale_summary_shard_dir, "scale_summary_shards");
verifyEqual(testCase, params.scale_band_micro_max_sec, 0.8);
verifyEqual(testCase, params.scale_band_motif_max_sec, 2.5);
verifyGreaterThan(testCase, params.max_scale_sec, 20);
verifyEqual(testCase, params.score_n_pcs, 12);
verifyEqual(testCase, params.score_variance_target_pct, 90);
verifyEqual(testCase, params.score_representation_mode, "multiresolution");
verifyGreaterThanOrEqual(testCase, params.summary_temporal_bins, 4);
verifyGreaterThanOrEqual(testCase, params.summary_dct_coeffs, 2);
verifyGreaterThanOrEqual(testCase, params.stability_bootstraps, 1);
verifyGreaterThan(testCase, params.embedding_dim_context_max_pcs, params.embedding_dim_context_min_pcs);
verifyEqual(testCase, params.primary_scale_rule, "stable_and_dimension_supported");
verifyGreaterThanOrEqual(testCase, params.scale_specific_max_anchors_per_scale, 1);
verifyGreaterThanOrEqual(testCase, params.pca_loading_stability_splits, 1);
verifyGreaterThanOrEqual(testCase, params.pca_loading_stability_threshold, 0);
verifyGreaterThan(testCase, params.event_turn_threshold_deg, 0);
verifyGreaterThanOrEqual(testCase, params.arena_sensitivity_n_embedding_pcs, 1);
end

function testRun06BandAdaptiveSummaryDefaults(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
params = load_multiscale_chunking_config(fullfile(repoRoot, 'config', 'multiscale_chunking_config.csv'));

verifyTrue(testCase, logical(params.summary_band_adaptive));
verifyEqual(testCase, params.summary_micro_temporal_bins, 6);
verifyEqual(testCase, params.summary_micro_dct_coeffs, 4);
verifyEqual(testCase, params.summary_motif_temporal_bins, 12);
verifyEqual(testCase, params.summary_motif_dct_coeffs, 8);
verifyEqual(testCase, params.summary_context_temporal_bins, 6);
verifyEqual(testCase, params.summary_context_dct_coeffs, 4);
verifyGreaterThan(testCase, params.summary_motif_temporal_bins, params.summary_micro_temporal_bins);
verifyGreaterThan(testCase, params.summary_motif_dct_coeffs, params.summary_micro_dct_coeffs);
end

function testRun06SmokeComparisonArtifactsIfPresent(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
baseRoot = fullfile(repoRoot, 'derived', 'run06_smoke_experiments', 'baseline_bins6_dct4');
motifRoot = fullfile(repoRoot, 'derived', 'run06_smoke_experiments', 'rich_bins12_dct8');
basePath = fullfile(baseRoot, 'embedding_dimension_audit.csv');
motifPath = fullfile(motifRoot, 'embedding_dimension_audit.csv');

if ~(isfile(basePath) && isfile(motifPath))
    runnerText = string(fileread(fullfile(repoRoot, 'chunks', 'run_multiscale_chunking_and_scale_selection.m')));
    verifyTrue(testCase, contains(runnerText, 'summary_motif_temporal_bins'));
    verifyTrue(testCase, contains(runnerText, 'summary_motif_dct_coeffs'));
    return
end

Base = readtable(basePath, 'TextType', 'string');
Motif = readtable(motifPath, 'TextType', 'string');
Base = Base(string(Base.initial_band) == "motif", :);
Motif = Motif(string(Motif.initial_band) == "motif", :);
verifyFalse(testCase, isempty(Base));
verifyFalse(testCase, isempty(Motif));
verifyGreaterThan(testCase, median(Motif.n_summary_dims), median(Base.n_summary_dims));
verifyGreaterThan(testCase, median(Motif.effective_dim, 'omitnan'), median(Base.effective_dim, 'omitnan'));
end

function testRun06FullEnvRoutesToProductionRoot(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
oldMode = getenv('RUN06_CHUNK_RUN_MODE');
oldOutput = getenv('RUN06_CHUNK_OUTPUT_DIR');
cleanup = onCleanup(@() i_restore_env(oldMode, oldOutput)); %#ok<NASGU>

setenv('RUN06_CHUNK_RUN_MODE', 'full');
setenv('RUN06_CHUNK_OUTPUT_DIR', '');
params = load_multiscale_chunking_config(fullfile(repoRoot, 'config', 'multiscale_chunking_config.csv'));

verifyEqual(testCase, params.run_mode, "full");
verifyEqual(testCase, params.output_dir, "derived/chunks_motif_discovery");
verifyFalse(testCase, contains(params.output_dir, "smoke"));
end

function testRun06RejectsExplicitFullModeWithSmokeRoot(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
oldMode = getenv('RUN06_CHUNK_RUN_MODE');
oldOutput = getenv('RUN06_CHUNK_OUTPUT_DIR');
cleanup = onCleanup(@() i_restore_env(oldMode, oldOutput)); %#ok<NASGU>

setenv('RUN06_CHUNK_RUN_MODE', 'full');
setenv('RUN06_CHUNK_OUTPUT_DIR', 'derived/chunks_motif_discovery_smoke');

verifyError(testCase, ...
    @() load_multiscale_chunking_config(fullfile(repoRoot, 'config', 'multiscale_chunking_config.csv')), ...
    'load_multiscale_chunking_config:BadFullOutputDir');
end

function testAnchorFinderUsesMaskAndFiniteFeatures(testCase)
Seq = i_mock_seq(240, 20);
Seq.validMask(55:85) = false;
Seq.X(150:170, 2) = NaN;

Anchor = find_condition_blind_chunk_anchors(Seq, [0.5; 1.0], ...
    'strideSec', 0.25, ...
    'anchorMode', "center", ...
    'minValidFrac', 0.90, ...
    'minFeatureFiniteFrac', 0.95, ...
    'requireAnchorFrameValid', true);

verifyFalse(testCase, any(Anchor.anchor_frame >= 55 & Anchor.anchor_frame <= 85));
verifyTrue(testCase, all(Anchor.min_scale_valid_frac >= 0.90));
verifyTrue(testCase, all(Anchor.min_scale_feature_finite_frac >= 0.95));
verifyTrue(testCase, all(Anchor.anchor_frame_valid));
end

function testAnchorSubsetSelectionIgnoresProvenanceLabels(testCase)
C = table();
C.session_index = [ones(8,1); 2*ones(8,1)];
C.raw_index = [101*ones(8,1); 102*ones(8,1)];
C.anchor_time_s = [(1:8)'; (1:8)'];
C.condition_id = ["A"; "B"; "C"; "D"; "E"; "F"; "G"; "H"; ...
    "H"; "G"; "F"; "E"; "D"; "C"; "B"; "A"];
C.arena_label = repmat(["big"; "small"], 8, 1);

S1 = select_condition_blind_anchor_subset(C, 6, 'minAnchorsPerSession', 2);
C.condition_id = flipud(C.condition_id);
C.arena_label = flipud(C.arena_label);
S2 = select_condition_blind_anchor_subset(C, 6, 'minAnchorsPerSession', 2);

verifyEqual(testCase, S1.raw_index, S2.raw_index);
verifyEqual(testCase, S1.anchor_time_s, S2.anchor_time_s);
verifyEqual(testCase, S1.session_index, S2.session_index);
verifyTrue(testCase, all(S1.anchor_selection_rule == ...
    "deterministic_time_even_by_session_then_global_no_condition_labels"));
end

function testRun06PaperScriptHasNoLocalFunctions(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
scriptPath = fullfile(repoRoot, 'paper', 'run_06_build_multiscale_chunks_and_select_scales.m');
txt = string(fileread(scriptPath));

verifyFalse(testCase, contains(txt, newline + "function "));
verifyTrue(testCase, contains(txt, "run_multiscale_chunking_and_scale_selection"));
end

function testRun06DeclaresScaleSurveyOutputs(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
runnerPath = fullfile(repoRoot, 'chunks', 'run_multiscale_chunking_and_scale_selection.m');
figPath = fullfile(repoRoot, 'chunks', 'make_run06_chunk_qc_figures.m');
runnerText = string(fileread(runnerPath));
figureText = string(fileread(figPath));

verifyTrue(testCase, contains(runnerText, "scale_session_anchor_coverage_audit.csv"));
verifyTrue(testCase, contains(runnerText, "scale_anchor_coverage_audit.csv"));
verifyTrue(testCase, contains(runnerText, "scale_summary_shard_manifest.csv"));
verifyTrue(testCase, contains(runnerText, "build_scale_summary_shard_from_anchor_manifest"));
verifyTrue(testCase, contains(runnerText, "local_summary_profile_for_scale"));
verifyTrue(testCase, contains(runnerText, "summary_profile"));
verifyTrue(testCase, contains(runnerText, "local_drop_summary_matrices"));
verifyFalse(testCase, contains(runnerText, "[ChunkSet, anchorManifest] = build_chunkset_from_anchor_manifest"));
verifyTrue(testCase, contains(runnerText, "embedding_dimension_audit.csv"));
verifyTrue(testCase, contains(runnerText, "scale_selection_stability.csv"));
verifyTrue(testCase, contains(runnerText, "primary_operational_scales.csv"));
verifyTrue(testCase, contains(runnerText, "primary_scale_specific_anchor_manifest.csv"));
verifyTrue(testCase, contains(runnerText, "primary_chunk_event_summary_audit.csv"));
verifyTrue(testCase, contains(runnerText, "pca_loading_stability.csv"));
verifyTrue(testCase, contains(runnerText, "arena_sensitivity_audit.csv"));
verifyTrue(testCase, contains(runnerText, "n_common_all_scale_candidate_anchors"));
verifyTrue(testCase, contains(runnerText, "SmokeProductionOutputBlocked"));
verifyTrue(testCase, contains(runnerText, "save_chunk_mat_artifact"));
verifyTrue(testCase, contains(runnerText, "Skipped ignored MAT artifact by config"));
verifyTrue(testCase, contains(figureText, "scale_specific_anchor_coverage_audit"));
verifyTrue(testCase, contains(figureText, "embedding_dimension_audit"));
verifyTrue(testCase, contains(figureText, "scale_selection_stability"));
verifyTrue(testCase, contains(figureText, "primary_operational_scale_decision"));
verifyTrue(testCase, contains(figureText, "motif_band_dimension_stability_audit"));
verifyTrue(testCase, contains(figureText, "pca_loading_stability"));
verifyTrue(testCase, contains(figureText, "primary_chunk_event_summary_audit"));
verifyTrue(testCase, contains(figureText, "arena_embedding_sensitivity_audit_only"));
verifyTrue(testCase, contains(figureText, "n_pcs_90pct"));
verifyTrue(testCase, contains(figureText, "nRows = numel(scales);"));
end

function testLegacyConditionAwareTokensNotPorted(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
files = [
    string(fullfile(repoRoot, 'paper', 'run_06_build_multiscale_chunks_and_select_scales.m'))
    string(fullfile(repoRoot, 'chunks', 'run_multiscale_chunking_and_scale_selection.m'))
    string(fullfile(repoRoot, 'chunks', 'select_condition_blind_anchor_subset.m'))
    string(fullfile(repoRoot, 'chunks', 'find_condition_blind_chunk_anchors.m'))
    string(fullfile(repoRoot, 'chunks', 'build_scale_summary_shard_from_anchor_manifest.m'))
    string(fullfile(repoRoot, 'chunks', 'summarize_multiresolution_chunks.m'))
    string(fullfile(repoRoot, 'chunks', 'build_primary_scale_specific_anchor_manifest.m'))
    string(fullfile(repoRoot, 'chunks', 'summarize_primary_chunk_events.m'))
    string(fullfile(repoRoot, 'chunks', 'make_run06_chunk_qc_figures.m'))
    string(fullfile(repoRoot, 'scaleinfo', 'score_multiscale_chunk_bank.m'))
    string(fullfile(repoRoot, 'scaleinfo', 'estimate_scale_selection_stability.m'))
    string(fullfile(repoRoot, 'scaleinfo', 'estimate_pca_loading_stability.m'))
    string(fullfile(repoRoot, 'scaleinfo', 'audit_arena_embedding_sensitivity.m'))
    ];
txt = "";
for i = 1:numel(files)
    txt = txt + newline + string(fileread(files(i)));
end

forbidden = ["minAnchorsPerCondition", "conditionQuota", "rare_aware", ...
    "rare_weighted", "bootstrap_stratum", "analysis_weight_raw"];
for i = 1:numel(forbidden)
    verifyFalse(testCase, contains(txt, forbidden(i)), ...
        "Forbidden legacy token was ported: " + forbidden(i));
end
end

function testMultiresolutionSummaryCompressesLongChunks(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));

Seq = i_mock_seq(2400, 80);
Sc = struct();
Sc.chunkSec = 30;
Sc.nFrames = 2400;
Sc.X = single(randn(10, 2400, 3));
Sc.Xraw = Sc.X;
Sc.valid = true(10, 2400);
Sc.meta = table(repmat(1, 10, 1), 'VariableNames', {'scale_index'});

ChunkSet = struct();
ChunkSet.featureNames = {'centroid_dist'; 'heading_diff_deg'; 'in_contact'};
ChunkSet.channelMeta = table( ...
    ["centroid_dist"; "heading_diff_deg_cos"; "in_contact"], ...
    ["centroid_dist"; "heading_diff_deg"; "in_contact"], ...
    ["continuous"; "circular_cos"; "boolean"], ...
    'VariableNames', {'ObsName','BaseFeature','ChannelType'});

[Xsummary, dimMeta, audit] = summarize_multiresolution_chunks(Sc, ChunkSet, ...
    'featureNames', string(ChunkSet.featureNames), ...
    'nTemporalBins', 4, ...
    'nDctCoeffs', 2, ...
    'includeDct', true);

verifySize(testCase, Xsummary, [10, height(dimMeta)]);
verifyLessThan(testCase, size(Xsummary, 2), 2400 * 3);
verifyGreaterThan(testCase, audit.compression_ratio_raw_to_summary, 20);
verifyTrue(testCase, all(dimMeta.labels_used_for_summary == "none"));
end

function Seq = i_mock_seq(T, fps)
Seq = struct();
Seq.time = (0:T-1)' ./ fps;
Seq.fps = fps;
Seq.validMask = true(T, 1);
Seq.X = [sin((1:T)' ./ 10), cos((1:T)' ./ 11), double(mod((1:T)', 7) == 0)];
end

function i_restore_env(oldMode, oldOutput)
setenv('RUN06_CHUNK_RUN_MODE', oldMode);
setenv('RUN06_CHUNK_OUTPUT_DIR', oldOutput);
end
