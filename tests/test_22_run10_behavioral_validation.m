function tests = test_22_run10_behavioral_validation
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repoRoot = string(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(repoRoot));
testCase.TestData.repoRoot = repoRoot;
envNames = ["RUN10_VALIDATION_RUN_MODE","RUN10_VALIDATION_PHASE", ...
    "RUN10_VALIDATION_OUTPUT_DIR","RUN10_RATINGS_INPUT_FILE"];
oldValues = strings(size(envNames));
for i = 1:numel(envNames)
    oldValues(i) = string(getenv(envNames(i)));
    setenv(envNames(i),'');
end
testCase.TestData.envNames = envNames;
testCase.TestData.oldValues = oldValues;
end

function teardownOnce(testCase)
for i = 1:numel(testCase.TestData.envNames)
    setenv(testCase.TestData.envNames(i), ...
        testCase.TestData.oldValues(i));
end
end

function testConfigRegistryAndThinPaperContract(testCase)
params = load_motif_candidate_behavioral_validation_config();
verifyEqual(testCase,params.run_mode,"smoke");
verifyEqual(testCase,params.phase,"full");
verifyEqual(testCase,height(params.feature_specs),16);
verifyEmpty(testCase,intersect(params.discovery_value_pose_nodes, ...
    params.validation_pose_nodes));
verifyEqual(testCase,params.expected_membership_sha256, ...
    "2e58d214683badb0b86e746bf32a7944da664b67243e778ea2d2638249268442");
verifyEqual(testCase,numel(params.expected_eligible_candidate_ids),9);

registry = load_paper_step_registry(testCase.TestData.repoRoot);
row = registry.step_id=="run_10";
verifyEqual(testCase,nnz(row),1);
verifyTrue(testCase,registry.smoke_safe(row));
paperPath = fullfile(testCase.TestData.repoRoot,'paper', ...
    'run_10_validate_motif_candidates_behaviorally.m');
paperText = string(fileread(paperPath));
verifyFalse(testCase,contains(paperText,newline+"function "));
verifyTrue(testCase,contains(paperText, ...
    "run_condition_blind_candidate_behavioral_validation"));
end

function testFreezeGateAndAdversarialPins(testCase)
params = load_motif_candidate_behavioral_validation_config();
Validation = validate_run10_candidate_freeze( ...
    testCase.TestData.repoRoot,params,'');
verifyTrue(testCase,all(Validation.gate.pass));
verifyEqual(testCase,height(Validation.membership),50548);
verifyEqual(testCase,numel(Validation.candidate_ids),42);
verifyEqual(testCase,numel(Validation.eligible_candidate_ids),9);

badHash = params;
badHash.expected_membership_sha256 = string(repmat('0',1,64));
verifyError(testCase,@() validate_run10_candidate_freeze( ...
    testCase.TestData.repoRoot,badHash,''), ...
    'validate_run10_candidate_freeze:FreezeMismatch');
badId = params;
badId.expected_candidate_freeze_id = "run09_candidates_wrong";
verifyError(testCase,@() validate_run10_candidate_freeze( ...
    testCase.TestData.repoRoot,badId,''), ...
    'validate_run10_candidate_freeze:FreezeMismatch');

tmp = string(tempname);
mkdir(tmp);
cleanup = onCleanup(@() i_remove(tmp)); %#ok<NASGU>
fields = ["membership_file","candidate_definition_file", ...
    "candidate_hierarchy_file","candidate_topology_file", ...
    "candidate_ambiguity_file","run09_manifest_file", ...
    "run09_freeze_validation_file"];
sourceRoot = resolve_repo_path(testCase.TestData.repoRoot, ...
    params.run09_output_dir);
for i = 1:numel(fields)
    copyfile(fullfile(sourceRoot,params.(fields(i))), ...
        fullfile(tmp,params.(fields(i))));
end
membershipPath = fullfile(tmp,params.membership_file);
[fid,message] = fopen(membershipPath,'a');
assert(fid>=0,message);
fprintf(fid,'\n');
fclose(fid);
modified = params;
modified.run09_output_dir = tmp;
verifyError(testCase,@() validate_run10_candidate_freeze( ...
    testCase.TestData.repoRoot,modified,''), ...
    'validate_run10_candidate_freeze:FreezeMismatch');
end

function testLineageFirewallAndValueBlindSampling(testCase)
params = load_motif_candidate_behavioral_validation_config();
Validation = validate_run10_candidate_freeze( ...
    testCase.TestData.repoRoot,params,'');
[Registry,Overlap] = build_run10_validation_feature_registry( ...
    testCase.TestData.repoRoot,params,Validation);
verifyTrue(testCase,all(Registry.included_in_primary_validation));
verifyFalse(testCase,any(Registry.raw_pose_node_overlap));
verifyFalse(testCase,any( ...
    Registry.exact_discovery_feature_name_overlap));
verifyFalse(testCase,any( ...
    Registry.renamed_or_deterministic_discovery_transform));
verifyGreaterThan(testCase,height(Overlap),10000);
verifyTrue(testCase,all( ...
    Overlap.prohibited_from_primary_validation));
verifyFalse(testCase,any( ...
    Overlap.used_for_run10_validation_status));
verifyTrue(testCase,all(Overlap.experimental_labels_used=="none"));

Sample1 = build_run10_validation_sample_manifest( ...
    testCase.TestData.repoRoot,params,Validation,Registry);
shuffled = Validation;
stream = RandStream('mt19937ar','Seed',123);
shuffled.membership = shuffled.membership( ...
    randperm(stream,height(shuffled.membership)),:);
Sample2 = build_run10_validation_sample_manifest( ...
    testCase.TestData.repoRoot,params,shuffled,Registry);
verifyEqual(testCase,sort(Sample1.graph_node_id( ...
    Sample1.sample_selected)),sort(Sample2.graph_node_id( ...
    Sample2.sample_selected)));
verifyTrue(testCase,all(Sample1.experimental_labels_used=="none"));
verifyFalse(testCase,any(Sample1.arena_used_for_sampling));
verifyFalse(testCase,any(Sample1.graph_stability_used_for_sampling));
verifyFalse(testCase,any(Sample1.annotation_used_for_sampling));
verifyFalse(testCase,any(Sample1.visual_appearance_used_for_sampling));
end

function testSyntheticFeatureExtractionIsDeterministic(testCase)
params = load_motif_candidate_behavioral_validation_config();
tmp = string(tempname);
mkdir(tmp);
cleanup = onCleanup(@() i_remove(tmp)); %#ok<NASGU>
posePath = fullfile(tmp,'pose.mat');
T = 200;
tracks = nan(T,13,2,2);
t = (1:T)';
for animal = 1:2
    for node = 1:13
        tracks(:,node,1,animal) = t*0.03+node*2+animal*20+ ...
            sin(t/11+node);
        tracks(:,node,2,animal) = node*3+animal*15+ ...
            cos(t/13+node/2);
    end
end
sessionPreproc = struct();
sessionPreproc.clean.tracks = tracks;
sessionPreproc.qc.badframes = false(T,2);
sessionPreproc.params.data.fps = 80;
sessionPreproc.params.data.pixel_size_mm = 0.50761421319797;
save(posePath,'sessionPreproc');

Registry = params.feature_specs;
Registry.included_in_primary_validation = true(height(Registry),1);
Validation = struct('membership_sha256', ...
    params.expected_membership_sha256);
Sample = table(1,"MC_M000750_C0001",1,true,1,2,0.2215,1, ...
    "session_synth",1,100,1.25,string(posePath),T,true,"none", ...
    'VariableNames', {'graph_node_id','motif_candidate_id', ...
    'candidate_local_index','eligible_for_behavioral_interpretation', ...
    'embedding_row_id','scale_index','chunk_sec','session_index', ...
    'session_id','raw_index','anchor_frame','anchor_time_s', ...
    'preprocess_output_file','source_frame_count','sample_selected', ...
    'experimental_labels_used'});
M1 = extract_run10_independent_behavioral_measures( ...
    testCase.TestData.repoRoot,params,Validation,Registry,Sample);
M2 = extract_run10_independent_behavioral_measures( ...
    testCase.TestData.repoRoot,params,Validation,Registry,Sample);
features = cellstr(params.feature_specs.feature_name);
verifyEqual(testCase,M1{:,features},M2{:,features},'AbsTol',0);
verifyEqual(testCase,M1.source_pose_sha256,M2.source_pose_sha256);
verifyTrue(testCase,all(M1.experimental_labels_used=="none"));
verifyTrue(testCase,all(M1.run05_to_run09_values_used=="none"));
end

function testGroupedFoldsAreDeterministicAndLeakFree(testCase)
params = load_motif_candidate_behavioral_validation_config();
nSessions = 12;
sessionIndex = repelem((1:nSessions)',5);
sessionId = "S"+string(sessionIndex);
eligible = true(numel(sessionIndex),1);
Measurement = table(sessionIndex,sessionId,eligible, ...
    'VariableNames', {'session_index','session_id', ...
    'eligible_for_automated_analysis'});
F1 = build_run10_blocked_validation_folds(Measurement,params);
F2 = build_run10_blocked_validation_folds(Measurement,params);
verifyEqual(testCase,F1,F2);
verifyEqual(testCase,numel(unique(F1.session_index)),height(F1));
verifyTrue(testCase,all(F1.session_in_exactly_one_test_fold));
verifyFalse(testCase,any(F1.candidate_labels_used_for_fold_assignment));
verifyTrue(testCase,all(F1.experimental_labels_used=="none"));
end

function testAbsentAndInvalidRatings(testCase)
params = load_motif_candidate_behavioral_validation_config();
tmp = string(tempname);
mkdir(tmp);
cleanup = onCleanup(@() i_remove(tmp)); %#ok<NASGU>
params.ratings_input_file = fullfile(tmp,'missing.csv');
[Ratings,Audit] = ingest_run10_blinded_ratings(tmp,params,struct());
verifyTrue(testCase,isempty(Ratings));
verifyEqual(testCase,Audit.ingestion_status,"awaiting_ratings");

bad = table("","RVW-X","candidate","",NaN,"","","","", ...
    NaN,NaN,3,"", ...
    'VariableNames', {'rater_id','review_id','rating_level','tile_label', ...
    'behavioral_coherence_score','interactive_impression', ...
    'contact_configuration','relative_orientation','motion_impression', ...
    'ambiguity_score','confidence_score', ...
    'candidate_overall_coherence_score','free_text_description'});
badPath = fullfile(tmp,'bad.csv');
writetable(bad,badPath);
params.ratings_input_file = badPath;
Packet = struct();
Packet.template = bad;
Packet.manifest = table("RVW-X",'VariableNames',{'review_id'});
verifyError(testCase,@() ingest_run10_blinded_ratings( ...
    tmp,params,Packet), ...
    'ingest_run10_blinded_ratings:MissingRater');
end

function testSmokeOutputsAreManifestedIfPresent(testCase)
root = fullfile(testCase.TestData.repoRoot,'derived', ...
    'motif_candidate_behavioral_validation_smoke');
manifestPath = fullfile(root,'run10_output_manifest.csv');
assumeTrue(testCase,isfile(manifestPath));
Manifest = readtable(manifestPath,'Delimiter',',','TextType','string');
verifyFalse(testCase,any(Manifest.record_type=="output_artifact" & ...
    contains(lower(Manifest.artifact_id),"membership")));
files = Manifest(ismember(Manifest.record_type,[ ...
    "frozen_or_audit_input","pose_measurement_input", ...
    "source_provenance","output_artifact","figure"]),:);
for i = 1:height(files)
    verifyTrue(testCase,isfile(files.artifact_path(i)));
    verifyEqual(testCase,compute_file_sha256(files.artifact_path(i)), ...
        files.sha256(i));
end
verifyTrue(testCase,all(Manifest.no_experimental_labels_used));
status = readtable(fullfile(root, ...
    'motif_candidate_behavioral_validation_status.csv'), ...
    'Delimiter',',','TextType','string');
verifyFalse(testCase,any(status.paper_claim_allowed));
verifyTrue(testCase,all(status.run10_validation_status( ...
    logical(status.run09_graph_interpretation_eligible))== ...
    "smoke_only_not_for_paper_claim"));
end

function i_remove(pathText)
if isfolder(pathText)
    rmdir(pathText,'s');
end
end
