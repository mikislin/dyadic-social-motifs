function release_run09_matlab_igraph_bridge(Bridge)
%RELEASE_RUN09_MATLAB_IGRAPH_BRIDGE Restore process path state.

clear mexIgraphDispatcher
if isstruct(Bridge) && isfield(Bridge, 'bridge_root') && ...
        contains(string(path), string(Bridge.bridge_root))
    rmpath(Bridge.bridge_root);
end
if isstruct(Bridge) && isfield(Bridge, 'old_system_path')
    setenv('PATH', Bridge.old_system_path);
end
rehash;
end
