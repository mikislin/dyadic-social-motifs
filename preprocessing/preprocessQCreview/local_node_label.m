function label = local_node_label(cfg, node)
labels = cfg.bodypoints.bodypoint_labels;
if node <= numel(labels)
    label = string(labels(node));
else
    label = "Point " + string(sprintf('%02d', node));
end
end