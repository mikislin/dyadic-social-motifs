function rgb = local_condition_rgb(styles, conditionId)
idx = find(styles.conditionIds == string(conditionId), 1, 'first');
if isempty(idx)
    rgb = [0.35 0.35 0.35];
else
    rgb = styles.colors(idx,:);
end
end