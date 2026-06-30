function G = local_sort_group_summary(G, styles)
conditionOrder = nan(height(G), 1);
for i = 1:height(G)
    conditionOrder(i) = find(styles.conditionIds == string(G.condition_id(i)), 1, 'first');
end
conditionOrder(~isfinite(conditionOrder)) = height(styles.table) + 1;
G.condition_order = conditionOrder;
sortVars = {'condition_order'};
if ismember('arena_condition_label', G.Properties.VariableNames)
    sortVars = [sortVars, {'arena_condition_label'}];
elseif ismember('plot_group_label', G.Properties.VariableNames)
    sortVars = [sortVars, {'plot_group_label'}];
end
G = sortrows(G, sortVars);
end