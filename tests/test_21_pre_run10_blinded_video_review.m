function tests = test_21_pre_run10_blinded_video_review
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repoRoot = string(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(repoRoot));
testCase.TestData.repoRoot = repoRoot;
end

function testConfigAndUnregisteredPreparationContract(testCase)
envNames = ["PRE_RUN10_VIDEO_REVIEW_MODE", ...
    "PRE_RUN10_VIDEO_REVIEW_OUTPUT_DIR", "PRE_RUN10_RENDER_VIDEOS"];
oldValues = i_capture_env(envNames);
cleanup = onCleanup(@() i_restore_env(envNames, oldValues)); %#ok<NASGU>
setenv('PRE_RUN10_VIDEO_REVIEW_MODE', 'smoke');
setenv('PRE_RUN10_VIDEO_REVIEW_OUTPUT_DIR', '');
setenv('PRE_RUN10_RENDER_VIDEOS', '0');
params = load_motif_candidate_video_review_config();

verifyEqual(testCase, params.run_mode, "smoke");
verifyEqual(testCase, params.active_candidate_limit, 1);
verifyEqual(testCase, params.tiles_per_candidate, 9);
verifyEqual(testCase, [params.central_tiles, ...
    params.core_random_tiles, params.boundary_random_tiles], [3 3 3]);
verifyEqual(testCase, params.expected_candidate_freeze_id, ...
    "run09_candidates_2e58d214683b_8817d942_0e956a2c");
verifyEqual(testCase, params.expected_membership_sha256, ...
    "2e58d214683badb0b86e746bf32a7944da664b67243e778ea2d2638249268442");
verifyFalse(testCase, params.render_videos);

registry = load_paper_step_registry(testCase.TestData.repoRoot);
verifyFalse(testCase, any(contains(registry.script_path, "pre_run10")));
paperPath = fullfile(testCase.TestData.repoRoot, 'paper', ...
    'pre_run10_prepare_blinded_motif_candidate_video_review.m');
paperText = string(fileread(paperPath));
verifyFalse(testCase, contains(paperText, newline + "function "));
verifyTrue(testCase, contains(paperText, ...
    "prepare_blinded_motif_candidate_video_review"));
end

function testFullSelectionIsDeterministicUniqueAndBlinded(testCase)
envNames = ["PRE_RUN10_VIDEO_REVIEW_MODE", ...
    "PRE_RUN10_VIDEO_REVIEW_OUTPUT_DIR", "PRE_RUN10_RENDER_VIDEOS"];
oldValues = i_capture_env(envNames);
cleanup = onCleanup(@() i_restore_env(envNames, oldValues)); %#ok<NASGU>
setenv('PRE_RUN10_VIDEO_REVIEW_MODE', 'full');
setenv('PRE_RUN10_VIDEO_REVIEW_OUTPUT_DIR', '');
setenv('PRE_RUN10_RENDER_VIDEOS', '0');
params = load_motif_candidate_video_review_config();
[Selection1, Summary1] = ...
    select_blinded_motif_candidate_video_examples( ...
    testCase.TestData.repoRoot, params);
[Selection2, Summary2] = ...
    select_blinded_motif_candidate_video_examples( ...
    testCase.TestData.repoRoot, params);

verifyEqual(testCase, height(Summary1), 9);
verifyEqual(testCase, height(Selection1), 81);
verifyEqual(testCase, Selection1.graph_node_id, Selection2.graph_node_id);
verifyEqual(testCase, Selection1.display_tile_index, ...
    Selection2.display_tile_index);
verifyEqual(testCase, Summary1, Summary2);
verifyEqual(testCase, numel(unique(Selection1.graph_node_id)), 81);
verifyEqual(testCase, groupsummary(Selection1, ...
    'motif_candidate_id').GroupCount, 9 * ones(9, 1));
verifyEqual(testCase, groupsummary(Selection1, ...
    {'motif_candidate_id', 'selection_stratum'}).GroupCount, ...
    3 * ones(27, 1));
verifyEqual(testCase, groupsummary(Selection1, ...
    'motif_candidate_id', 'numunique', ...
    'session_index').numunique_session_index, 9 * ones(9, 1));
verifyEqual(testCase, sort(groupsummary(Selection1, ...
    'motif_candidate_id', 'numunique', ...
    'display_tile_index').numunique_display_tile_index), ...
    9 * ones(9, 1));
verifyTrue(testCase, all(Selection1.experimental_labels_used == "none"));
verifyTrue(testCase, all( ...
    Selection1.selected_after_candidate_freeze));
verifyFalse(testCase, any(Selection1.used_for_membership));
verifyFalse(testCase, ...
    any(Selection1.used_for_interpretation_eligibility));
verifyFalse(testCase, ...
    any(Selection1.visual_appearance_used_for_selection));
verifyEqual(testCase, unique(Selection1.frozen_membership_sha256), ...
    params.expected_membership_sha256);

forbiddenColumns = ["condition", "cohort", "group", "drug", ...
    "genotype", "perturbation", "treatment", "outcome", ...
    "response", "arena", "annotation", "umap", "rare_stratum", ...
    "event"];
actualColumns = lower(string(Selection1.Properties.VariableNames));
verifyFalse(testCase, any(contains(actualColumns, forbiddenColumns)));
end

function testSelectionAndRendererDoNotReadInterpretationMetadata(testCase)
selectionText = lower(string(fileread(fullfile( ...
    testCase.TestData.repoRoot, 'validation', ...
    'select_blinded_motif_candidate_video_examples.m'))));
renderText = lower(string(fileread(fullfile( ...
    testCase.TestData.repoRoot, 'validation', ...
    'render_blinded_motif_candidate_video_mosaics.m'))));

forbiddenInputs = ["motif_candidate_annotation", ...
    "motif_candidate_behavioral_profile", ...
    "motif_candidate_event_profile", "graph_plot_pc", "umap"];
verifyFalse(testCase, any(contains(selectionText, forbiddenInputs)));
verifyFalse(testCase, any(contains(renderText, forbiddenInputs)));
verifyFalse(testCase, contains(renderText, "examples.session_id"));
verifyFalse(testCase, contains(renderText, "examples.selection_stratum"));
verifyFalse(testCase, contains(renderText, "provisional_annotation"));
end

function testSmokeContactSheetRenders(testCase)
envNames = ["PRE_RUN10_VIDEO_REVIEW_MODE", ...
    "PRE_RUN10_VIDEO_REVIEW_OUTPUT_DIR", "PRE_RUN10_RENDER_VIDEOS"];
oldValues = i_capture_env(envNames);
cleanupEnv = onCleanup(@() i_restore_env( ...
    envNames, oldValues)); %#ok<NASGU>
setenv('PRE_RUN10_VIDEO_REVIEW_MODE', 'smoke');
setenv('PRE_RUN10_VIDEO_REVIEW_OUTPUT_DIR', '');
setenv('PRE_RUN10_RENDER_VIDEOS', '0');
params = load_motif_candidate_video_review_config();
[Selection, Summary] = select_blinded_motif_candidate_video_examples( ...
    testCase.TestData.repoRoot, params);
outRoot = string(tempname);
mkdir(outRoot);
cleanupDir = onCleanup(@() i_remove_dir(outRoot)); %#ok<NASGU>
Audit = render_blinded_motif_candidate_video_mosaics( ...
    testCase.TestData.repoRoot, Selection, Summary, ...
    outRoot, params);

verifyEqual(testCase, height(Audit), 1);
verifyEqual(testCase, Audit.render_status, "success");
verifyFalse(testCase, Audit.video_rendered);
verifyTrue(testCase, isfile(Audit.contact_sheet_file));
verifyEqual(testCase, compute_file_sha256(Audit.contact_sheet_file), ...
    Audit.contact_sheet_sha256);
verifyGreaterThan(testCase, Audit.contact_sheet_file_bytes, 0);
verifyEqual(testCase, Audit.experimental_labels_used, "none");
end

function oldValues = i_capture_env(names)
oldValues = strings(size(names));
for i = 1:numel(names)
    oldValues(i) = string(getenv(names(i)));
end
end

function i_restore_env(names, values)
for i = 1:numel(names)
    setenv(names(i), values(i));
end
end

function i_remove_dir(pathText)
if isfolder(pathText)
    rmdir(pathText, 's');
end
end
