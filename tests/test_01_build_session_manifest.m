function tests = test_01_build_session_manifest
tests = functiontests(localfunctions);
end

function testManifestClassifiesSessionShapes(testCase)
tmp = tempname;
mkdir(tmp);
cleanup = onCleanup(@() local_remove_dir(tmp));

local_write_session(tmp, 1, rand(10,13,2), ones(10,13), 'session_single');
local_write_session(tmp, 2, rand(20,13,2,2), ones(20,13,2), 'session_dyad');
local_write_session(tmp, 3, rand(15,13,2,3), ones(15,13,3), 'session_three');
tracksSpurious = rand(15,13,2,3);
tracksSpurious(:,:,:,3) = NaN;
scoresSpurious = ones(15,13,3);
local_write_session(tmp, 4, tracksSpurious, scoresSpurious, 'session_spurious_three');

outPath = fullfile(tmp, 'manifest.csv');
M = build_session_manifest(tmp, '', ...
    'loadDbaseMetadata', false, ...
    'writePath', outPath, ...
    'verbose', false);

verifyEqual(testCase, height(M), 4);
verifyEqual(testCase, M.raw_index', [1 2 3 4]);
verifyEqual(testCase, M.n_animals', [1 2 3 3]);
verifyEqual(testCase, M.effective_n_animals', [1 2 3 2]);
verifyFalse(testCase, M.include_block1_dyadic(1));
verifyTrue(testCase, M.include_block1_dyadic(2));
verifyFalse(testCase, M.include_block1_dyadic(3));
verifyTrue(testCase, M.include_block1_dyadic(4));
verifyFalse(testCase, M.include_motif_discovery(1));
verifyTrue(testCase, M.include_motif_discovery(2));
verifyFalse(testCase, M.include_motif_discovery(3));
verifyTrue(testCase, M.include_motif_discovery(4));
verifyEqual(testCase, M.motif_exclusion_reason(1), "not_two_usable_animals");
verifyEqual(testCase, M.motif_exclusion_reason(2), "");
verifyTrue(testCase, M.include_block2_egocentric(1));
verifyTrue(testCase, M.include_block2_egocentric(2));
verifyFalse(testCase, M.include_block2_egocentric(3));
verifyTrue(testCase, M.include_block2_egocentric(4));
verifyEqual(testCase, M.exclusion_reason(3), "more_than_two_usable_animals_not_supported_current_pipeline");
verifyEqual(testCase, M.animal_qc_status(4), "reduce_to_dyad");
verifyEqual(testCase, M.animal_keep_indices(4), "1;2");
verifyTrue(testCase, isfile(outPath));
end

function testSelectAnalysisCohortExcludesAnesthetizedMotifSessions(testCase)
M = table();
M.include_motif_discovery = [true; false; false];
M.include_block2_egocentric = [true; true; true];
M.paper_include = [true; true; true];
M.social_context = ["freely_moving_dyad"; "anesthetized_partner"; "single"];
M.mouse_type_1 = ["WT"; "MUT"; "WT"];
M.mouse_type_2 = ["WT"; "ANES"; ""];
M.mouse_type_3 = [""; ""; ""];

verifyEqual(testCase, select_analysis_cohort(M, "motif_discovery"), ...
    [true; false; false]);
verifyEqual(testCase, select_analysis_cohort(M, "block2_egocentric"), ...
    [true; true; true]);
verifyEqual(testCase, select_analysis_cohort(M, "anesthetized_context"), ...
    [false; true; false]);
end

function testPaperCohortLabelsArePublicationSafe(testCase)
M = table();
M.raw_index = [21; 32; 53; 130; 205];
M.session_file = "session_" + compose('%04d', M.raw_index) + ".mat";
M.condition = ["WT_WT"; "M_WT"; "M_A"; "WT_DCZPFC"; "M_DCZcDCN"];
M.arena = ["B"; "B"; "B"; "S"; "S"];
M.arena_label = ["big"; "big"; "big"; "small"; "small"];
M.mouse_type_1 = ["WT"; "MUT"; "MUT"; "WT"; "MUT"];
M.mouse_type_2 = ["WT"; "WT"; "ANES"; "WT_DREADD_PFC"; "WT"];

M = apply_paper_cohort_definitions(M);

verifyEqual(testCase, M.analysis_group_id(1), "cohort1|big|WT_WT");
verifyEqual(testCase, M.condition_label(1), "WT Pair");
verifyEqual(testCase, M.condition_label(2), "Mixed Pair");
verifyEqual(testCase, M.cohort_label(4), "Cohort 3 small arena WT DCZ PFC");
verifyEqual(testCase, M.defined_contrast_label(5), ...
    "Cohort 3 DCN chemogenetic extension");
verifyTrue(testCase, all(M.label_text_safe_for_matlab));

labelVars = string(M.Properties.VariableNames);
labelVars = labelVars(endsWith(labelVars, "_label") | labelVars == "plot_group_label");
for v = labelVars
    verifyFalse(testCase, any(contains(string(M.(v)), "_")), ...
        sprintf('Display label column %s contains underscores.', v));
end
end

function testPreprocessingConfigLoadsSleapBodypartMetadata(testCase)
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
cfg = load_preprocessing_pipeline_config( ...
    fullfile(repoRoot, 'config', 'preprocessing_pipeline_config.csv'));

expectedNames = ["nose","neck","Lear","Rear","LFpaw","RFpaw","LHpaw", ...
    "RHpaw","body","tailBase","tailMid","tailTip","bodyMid"];
expectedLabels = ["Nose","Neck","L ear","R ear","LF paw","RF paw", ...
    "LH paw","RH paw","Body","Tail base","Tail mid","Tail tip","Body mid"];
expectedFields = ["nose","neck","leftEar","rightEar","lfPaw","rfPaw", ...
    "lhPaw","rhPaw","body","tailBase","tailMid","tailTip","midBody"];

verifyEqual(testCase, cfg.bodypoints.bodypoint_names, expectedNames);
verifyEqual(testCase, cfg.bodypoints.bodypoint_labels, expectedLabels);
verifyEqual(testCase, cfg.bodypoints.bodypoint_matlab_fields, expectedFields);
verifyEqual(testCase, cfg.bodypoints.bodypoint_metadata.node_index', 1:13);
verifyTrue(testCase, isfile(cfg.bodypoints.bodypoint_metadata_path));
end

function local_write_session(rootDir, idx, tracks, scores, sessionId)
sessionRaw = struct();
sessionRaw.SLEAPtracks = tracks;
sessionRaw.SLEAPscores = scores;
sessionRaw.time = (0:size(tracks,1)-1)' ./ 80;
sessionRaw.excludedFrames = false(size(tracks,1),1);
sessionRaw.session_id = sessionId;
save(fullfile(rootDir, sprintf('session_%04d.mat', idx)), 'sessionRaw');
end

function local_remove_dir(pathName)
if isfolder(pathName)
    rmdir(pathName, 's');
end
end
