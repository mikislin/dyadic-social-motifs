function P = validate_preprocessing_params(P)
%VALIDATE_PREPROCESSING_PARAMS Fill defaults and validate preprocessing params.
%
% This version merges user params into defaults first, which avoids the long
% chain of duplicated isfield checks and keeps validator behavior aligned
% with default_preprocessing_params.

if nargin < 1 || isempty(P)
    P = default_preprocessing_params();
    return
end

P = merge_structs(default_preprocessing_params(), P);

%% Meta / data
assert(isfinite(P.data.fps) && P.data.fps > 0, 'P.data.fps must be > 0');
assert(isfinite(P.data.pixel_size_mm) && P.data.pixel_size_mm > 0, ...
    'P.data.pixel_size_mm must be > 0');
assert(P.data.n_coords_expected == 2, 'P.data.n_coords_expected must be 2');

%% Confidence
validModes = {'fixed','fixed_by_node','adaptive_per_node'};
assert(any(strcmpi(char(P.confidence.mode), validModes)), ...
    'P.confidence.mode must be one of: fixed, fixed_by_node, adaptive_per_node');
assert(isscalar(P.confidence.threshold) && P.confidence.threshold >= 0 && P.confidence.threshold <= 1, ...
    'P.confidence.threshold must be in [0,1]');
assert(isscalar(P.confidence.prctile) && P.confidence.prctile >= 0 && P.confidence.prctile <= 100, ...
    'P.confidence.prctile must be in [0,100]');
assert(isscalar(P.confidence.min_threshold) && P.confidence.min_threshold >= 0 && P.confidence.min_threshold <= 1, ...
    'P.confidence.min_threshold must be in [0,1]');
assert(isscalar(P.confidence.max_mask_fraction) && P.confidence.max_mask_fraction >= 0 && P.confidence.max_mask_fraction <= 1, ...
    'P.confidence.max_mask_fraction must be in [0,1]');
assert(isscalar(P.confidence.relaxed_threshold) && P.confidence.relaxed_threshold >= 0 && P.confidence.relaxed_threshold <= 1, ...
    'P.confidence.relaxed_threshold must be in [0,1]');
if ~isempty(P.confidence.threshold_by_node)
    assert(isvector(P.confidence.threshold_by_node), 'P.confidence.threshold_by_node must be a vector or []');
    assert(all(isfinite(P.confidence.threshold_by_node)), 'P.confidence.threshold_by_node must be finite');
    assert(all(P.confidence.threshold_by_node >= 0 & P.confidence.threshold_by_node <= 1), ...
        'P.confidence.threshold_by_node values must be in [0,1]');
end

%% Jump
validInterpMethods = {'linear','nearest','pchip','makima','spline'};
assert(isscalar(P.jump.diff_quantile) && P.jump.diff_quantile > 0 && P.jump.diff_quantile < 1, ...
    'P.jump.diff_quantile must be in (0,1)');
assert(isscalar(P.jump.max_iter) && P.jump.max_iter >= 1 && mod(P.jump.max_iter,1) == 0, ...
    'P.jump.max_iter must be a positive integer');
assert(isscalar(P.jump.min_valid_points) && P.jump.min_valid_points >= 2, ...
    'P.jump.min_valid_points must be >= 2');
assert(isscalar(P.jump.max_disp_mm_per_frame) && isfinite(P.jump.max_disp_mm_per_frame) && P.jump.max_disp_mm_per_frame > 0, ...
    'P.jump.max_disp_mm_per_frame must be > 0');
assert(isscalar(P.jump.max_disp_px_per_frame) && isfinite(P.jump.max_disp_px_per_frame) && P.jump.max_disp_px_per_frame > 0, ...
    'P.jump.max_disp_px_per_frame must be > 0');
assert(isscalar(P.jump.dilate_frames) && P.jump.dilate_frames >= 0 && mod(P.jump.dilate_frames,1) == 0, ...
    'P.jump.dilate_frames must be a nonnegative integer');
assert(any(strcmpi(char(P.jump.method), validInterpMethods)), ...
    'P.jump.method must be one of: linear, nearest, pchip, makima, spline');
assert(any(strcmpi(char(P.jump.interpolate_method), validInterpMethods)), ...
    'P.jump.interpolate_method must be one of: linear, nearest, pchip, makima, spline');
if ~isempty(P.jump.max_disp_mm_per_frame_by_node)
    assert(isvector(P.jump.max_disp_mm_per_frame_by_node) && all(P.jump.max_disp_mm_per_frame_by_node > 0), ...
        'P.jump.max_disp_mm_per_frame_by_node must contain positive values');
end
if ~isempty(P.jump.max_disp_px_per_frame_by_node)
    assert(isvector(P.jump.max_disp_px_per_frame_by_node) && all(P.jump.max_disp_px_per_frame_by_node > 0), ...
        'P.jump.max_disp_px_per_frame_by_node must contain positive values');
end

%% Gaps
assert(isscalar(P.gaps.max_gap_frames) && P.gaps.max_gap_frames >= 0 && mod(P.gaps.max_gap_frames,1) == 0, ...
    'P.gaps.max_gap_frames must be a nonnegative integer');
assert(isscalar(P.gaps.min_points_for_interp) && P.gaps.min_points_for_interp >= 2 && mod(P.gaps.min_points_for_interp,1) == 0, ...
    'P.gaps.min_points_for_interp must be an integer >= 2');
if ~isempty(P.gaps.max_gap_frames_by_node)
    assert(isvector(P.gaps.max_gap_frames_by_node) && all(P.gaps.max_gap_frames_by_node >= 0) && ...
        all(mod(P.gaps.max_gap_frames_by_node,1) == 0), ...
        'P.gaps.max_gap_frames_by_node must contain nonnegative integers');
end
assert(any(strcmpi(char(P.gaps.method), validInterpMethods)), ...
    'P.gaps.method must be one of: linear, nearest, pchip, makima, spline');

%% Smoothing
validSmoothMethods = {'lowpass','none'};
assert(any(strcmpi(char(P.smooth.method), validSmoothMethods)), ...
    'P.smooth.method must be one of: lowpass, none');
assert(isscalar(P.smooth.low_pass_frq) && P.smooth.low_pass_frq > 0, ...
    'P.smooth.low_pass_frq must be > 0');
assert(isscalar(P.smooth.n_pad) && P.smooth.n_pad >= 0 && mod(P.smooth.n_pad,1) == 0, ...
    'P.smooth.n_pad must be a nonnegative integer');
assert(isscalar(P.smooth.min_segment_length) && P.smooth.min_segment_length >= 1 && mod(P.smooth.min_segment_length,1) == 0, ...
    'P.smooth.min_segment_length must be a positive integer');
assert(isscalar(P.smooth.repair_break_halfwidth) && P.smooth.repair_break_halfwidth >= 0 && mod(P.smooth.repair_break_halfwidth,1) == 0, ...
    'P.smooth.repair_break_halfwidth must be a nonnegative integer');
if ~isempty(P.smooth.low_pass_frq_by_node)
    assert(isvector(P.smooth.low_pass_frq_by_node) && all(P.smooth.low_pass_frq_by_node > 0), ...
        'P.smooth.low_pass_frq_by_node must contain positive values');
end

%% QC
fracNames = {'frame_bad_node_frac_thresh','max_interp_frac_per_frame','max_jump_frac_per_frame', ...
    'max_geom_frac_per_frame','max_arena_frac_per_frame','max_lowconf_frac_per_frame'};
for i = 1:numel(fracNames)
    v = P.qc.(fracNames{i});
    assert(isscalar(v) && v >= 0 && v <= 1, 'P.qc.%s must be in [0,1]', fracNames{i});
end
assert(isscalar(P.qc.expand_severe_badframes) && P.qc.expand_severe_badframes >= 0 && mod(P.qc.expand_severe_badframes,1) == 0, ...
    'P.qc.expand_severe_badframes must be a nonnegative integer');
assert(numel(P.qc.body_length_nodes) == 2, 'P.qc.body_length_nodes must have 2 entries');
assert(size(P.qc.geometry_node_pairs,2) == 2, 'P.qc.geometry_node_pairs must be N x 2');
assert(numel(P.qc.body_length_prctile_range) == 2 && P.qc.body_length_prctile_range(1) < P.qc.body_length_prctile_range(2), ...
    'P.qc.body_length_prctile_range must be [low high]');
assert(numel(P.qc.geometry_prctile_range) == 2 && P.qc.geometry_prctile_range(1) < P.qc.geometry_prctile_range(2), ...
    'P.qc.geometry_prctile_range must be [low high]');
assert(numel(P.qc.arena_percentile) == 2 && P.qc.arena_percentile(1) < P.qc.arena_percentile(2), ...
    'P.qc.arena_percentile must be [low high]');
assert(isscalar(P.qc.arena_anchor_node) && P.qc.arena_anchor_node >= 1 && mod(P.qc.arena_anchor_node,1) == 0, ...
    'P.qc.arena_anchor_node must be a positive integer');
assert(isscalar(P.qc.arena_margin_px) && P.qc.arena_margin_px >= 0, ...
    'P.qc.arena_margin_px must be >= 0');

%% Output / debug
assert(isscalar(P.output.audit_plot_dpi) && P.output.audit_plot_dpi > 0, ...
    'P.output.audit_plot_dpi must be > 0');

end

function out = merge_structs(defaults, user)
out = defaults;
if ~isstruct(user)
    error('validate_preprocessing_params:BadInput', 'Input P must be a struct');
end
fields = fieldnames(user);
for i = 1:numel(fields)
    f = fields{i};
    if isfield(defaults, f) && isstruct(defaults.(f)) && isstruct(user.(f))
        out.(f) = merge_structs(defaults.(f), user.(f));
    else
        out.(f) = user.(f);
    end
end
end
