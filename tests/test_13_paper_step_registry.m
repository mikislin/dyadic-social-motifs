function tests = test_13_paper_step_registry
tests = functiontests(localfunctions);
end

function testPaperStepRegistryLoadsAndValidates(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));

registry = load_paper_step_registry(repoRoot);

required = ["step_id", "script_path", "enabled", "default_audit_mode", ...
    "smoke_safe", "smoke_env", "full_env", "expected_outputs", ...
    "description", "script_abs_path"];
verifyTrue(testCase, all(ismember(required, string(registry.Properties.VariableNames))));
verifyEqual(testCase, numel(unique(registry.step_id)), height(registry));
verifyTrue(testCase, all(ismember("run_0" + string(1:9), registry.step_id)));
verifyTrue(testCase, all(isfile(registry.script_abs_path)));
verifyTrue(testCase, all(startsWith(registry.script_path, "paper/")));
verifyTrue(testCase, all(ismember(registry.default_audit_mode, ["smoke", "full"])));
end

function testSmokeRegistryProtectsSharedOutputFeatureStep(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));

registry = load_paper_step_registry(repoRoot);

verifyTrue(testCase, registry.smoke_safe(registry.step_id == "run_02"));
verifyFalse(testCase, registry.smoke_safe(registry.step_id == "run_05"));
verifyEqual(testCase, registry.smoke_env(registry.step_id == "run_05"), ...
    "RUN05_FEATURE_RUN_MODE=smoke;RUN05_SKIP_EXISTING=true");
verifyEqual(testCase, registry.full_env(registry.step_id == "run_05"), ...
    "RUN05_FEATURE_RUN_MODE=full");
verifyTrue(testCase, registry.smoke_safe(registry.step_id == "run_06"));
verifyEqual(testCase, registry.smoke_env(registry.step_id == "run_06"), ...
    "RUN06_ANCHOR_MANIFEST_MODE=primary;RUN06_CHUNK_RUN_MODE=smoke;RUN06_CHUNK_OUTPUT_DIR=derived/chunks_motif_discovery_smoke");
verifyEqual(testCase, registry.full_env(registry.step_id == "run_06"), ...
    "RUN06_ANCHOR_MANIFEST_MODE=primary;RUN06_CHUNK_RUN_MODE=full;RUN06_CHUNK_OUTPUT_DIR=derived/chunks_motif_discovery");
verifyTrue(testCase, registry.smoke_safe(registry.step_id == "run_07"));
verifyEqual(testCase, registry.smoke_env(registry.step_id == "run_07"), ...
    "RUN07_ANCHOR_MANIFEST_MODE=primary;RUN07_EMBEDDING_RUN_MODE=smoke;RUN07_EMBEDDING_OUTPUT_DIR=derived/embedding_motif_discovery_smoke");
verifyEqual(testCase, registry.full_env(registry.step_id == "run_07"), ...
    "RUN07_ANCHOR_MANIFEST_MODE=primary;RUN07_EMBEDDING_RUN_MODE=full;RUN07_EMBEDDING_OUTPUT_DIR=derived/embedding_motif_discovery");
verifyTrue(testCase, registry.smoke_safe(registry.step_id == "run_08"));
verifyEqual(testCase, registry.smoke_env(registry.step_id == "run_08"), ...
    "RUN08_ANCHOR_MANIFEST_MODE=primary;RUN08_GRAPH_RUN_MODE=smoke;RUN08_GRAPH_OUTPUT_DIR=derived/graph_motif_discovery_smoke;RUN08_EMBEDDING_INPUT_DIR=derived/embedding_motif_discovery_smoke");
verifyEqual(testCase, registry.full_env(registry.step_id == "run_08"), ...
    "RUN08_ANCHOR_MANIFEST_MODE=primary;RUN08_GRAPH_RUN_MODE=full;RUN08_GRAPH_OUTPUT_DIR=derived/graph_motif_discovery;RUN08_EMBEDDING_INPUT_DIR=derived/embedding_motif_discovery");
verifyTrue(testCase, registry.smoke_safe(registry.step_id == "run_09"));
verifyEqual(testCase, registry.smoke_env(registry.step_id == "run_09"), ...
    "RUN09_CANDIDATE_RUN_MODE=smoke;RUN09_CANDIDATE_OUTPUT_DIR=derived/motif_candidate_discovery_smoke");
verifyEqual(testCase, registry.full_env(registry.step_id == "run_09"), ...
    "RUN09_CANDIDATE_RUN_MODE=full;RUN09_CANDIDATE_OUTPUT_DIR=derived/motif_candidate_discovery");
end

function testStepCallAuditToolExists(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
toolPath = fullfile(repoRoot, 'tools', 'run_step_call_audit.m');

verifyTrue(testCase, isfile(toolPath));
toolText = string(fileread(toolPath));
verifyTrue(testCase, contains(toolText, 'paper_step_registry.csv'));
verifyTrue(testCase, contains(toolText, 'step_function_call_audit.csv'));
verifyTrue(testCase, contains(toolText, 'function_inventory_with_step_usage.csv'));
end
