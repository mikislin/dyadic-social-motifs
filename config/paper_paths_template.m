function paths = paper_paths_template(varargin)
%PAPER_PATHS_TEMPLATE Backward-compatible wrapper for paper_paths.
%
% 
% See config/paper_paths.m for environment-variable and local override options.

paths = paper_paths(varargin{:});
end
