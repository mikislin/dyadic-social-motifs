function tests = test_18_matlab_igraph_backend
tests = functiontests(localfunctions);
end

function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));
end

function testInstalledBackendAuditIsDeterministicAndRejectsWeights(testCase)
assumeNotEmpty(testCase, which('igraph.cluster'), ...
    'matlab-igraph is not installed on this MATLAB path.');

outRoot = string(tempname);
mkdir(outRoot);
cleanup = onCleanup(@() i_remove_dir(outRoot)); %#ok<NASGU>
Audit = validate_matlab_igraph_backend(outRoot);

weighted = Audit.fixtures(Audit.fixtures.check_id == ...
    "weighted_cpm_partition", :);
unitControl = Audit.fixtures(Audit.fixtures.check_id == ...
    "unit_weight_cpm_control", :);
seedCheck = Audit.fixtures(Audit.fixtures.check_id == ...
    "same_seed_exact_reproduction", :);

verifyFalse(testCase, weighted.pass);
verifyTrue(testCase, unitControl.pass);
verifyTrue(testCase, seedCheck.pass);
verifyLessThan(testCase, weighted.adapter_or_returned_quality, ...
    weighted.independent_reference_quality);
verifyEqual(testCase, weighted.membership, unitControl.membership);

weightedCapability = Audit.capabilities(Audit.capabilities.capability == ...
    "weighted_leiden_cpm", :);
verifyFalse(testCase, weightedCapability.validated_for_run09);
verifyFalse(testCase, Audit.decision.matlab_igraph_selected);

verifyEqual(testCase, height(Audit.manifest), 8);
verifyTrue(testCase, all(strlength(Audit.manifest.sha256) == 64));
verifyTrue(testCase, isfile(Audit.manifest_path));
end

function i_remove_dir(pathText)
if isfolder(pathText)
    rmdir(pathText, 's');
end
end
