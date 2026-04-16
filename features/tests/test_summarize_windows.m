function tests = test_summarize_windows
tests = functiontests(localfunctions);
end

function testWindowSummaries(testCase)
T = 100;
dyad = struct();
dyad.time_s = (0:T-1)'/80;
dyad.frameMask = true(T,1);
dyad.featureNames = {'centroid_dist','in_contact','heading_diff_deg'};
dyad.X = [linspace(10,20,T)' repmat([zeros(10,1); ones(10,1)],5,1) linspace(-90,90,T)'];
dyad.featureMeta = table( ...
    dyad.featureNames(:), {'distance';'contact';'orientation'}, [false;false;false], [false;false;true], [false;true;false], ["log1p";"binary";"circular"], ...
    'VariableNames', {'Name','Family','IsDirected','IsCircular','IsBoolean','TransformHint'});

W = summarize_windows(dyad, 80, [0.5 1.0]);
verifyFalse(testCase, isempty(W.table));
verifyTrue(testCase, any(contains(W.table.Properties.VariableNames, 'centroid_dist__mean')));
verifyTrue(testCase, any(contains(W.table.Properties.VariableNames, 'in_contact__occupancy')));
verifyTrue(testCase, any(contains(W.table.Properties.VariableNames, 'heading_diff_deg__circmean')));
end
