function P = default_preprocessing_params()
%DEFAULT_PREPROCESSING_PARAMS Default parameters for dyadic SLEAP preprocessing.
%
% Returns a nested struct used by preprocess_session and helper functions.
% This version removes duplicated assignments and keeps field names aligned
% with downstream consumers.

P.meta.version = 'v0.4.2';
P.meta.created_by = 'MK';
P.meta.description = 'Preprocessing4dyadicInteraction';

%% Data
P.data.fps = 80;
P.data.pixel_size_mm = 1/1.97;
P.data.noi = [1,2,5,6,7,8,10,12];
P.data.n_animals_expected = 2;
P.data.n_coords_expected = 2;
P.data.n_nodes_expected = 13;

%% Confidence masking
P.confidence.enabled = true;
P.confidence.require_scores = false;
P.confidence.mode = 'adaptive_per_node';   % 'fixed', 'fixed_by_node', 'adaptive_per_node'
P.confidence.prctile_method = 'prctile';
P.confidence.threshold = 0.2;
P.confidence.threshold_by_node = [];
P.confidence.adaptive_per_node = true;
P.confidence.prctile = 2;
P.confidence.min_threshold = 0.03;
P.confidence.max_mask_fraction = 0.50;
P.confidence.relaxed_threshold = 0.01;

%% Jump filtering
P.jump.enabled = true;
P.jump.diff_quantile = 0.995;
P.jump.max_iter = 20;
P.jump.min_valid_points = 10;
P.jump.use_excludedFrames_mask = true;
P.jump.method = 'linear';
P.jump.interpolate_method = 'linear';
P.jump.dilate_frames = 1;
P.jump.max_disp_mm_per_frame = 25;
P.jump.max_disp_px_per_frame = P.jump.max_disp_mm_per_frame / P.data.pixel_size_mm;
P.jump.max_disp_mm_per_frame_by_node = [];
P.jump.max_disp_px_per_frame_by_node = [];

%% Gap interpolation
P.gaps.enabled = true;
P.gaps.method = 'linear';
P.gaps.max_gap_frames = 8;
P.gaps.fill_ends = false;
P.gaps.min_points_for_interp = 2;
P.gaps.max_gap_frames_by_node = [];
if P.data.n_nodes_expected >= 13
    P.gaps.max_gap_frames_by_node = 8 * ones(1, P.data.n_nodes_expected);
    P.gaps.max_gap_frames_by_node([9 10 13]) = 15;
    P.gaps.max_gap_frames_by_node([1 2]) = 8;
    P.gaps.max_gap_frames_by_node([5 6 7 8]) = 5;
end

%% Smoothing
P.smooth.enabled = true;
P.smooth.method = 'lowpass';
P.smooth.low_pass_frq = 8;
P.smooth.n_pad = 20;
P.smooth.min_segment_length = 2 * P.smooth.n_pad + 5;
P.smooth.low_pass_frq_by_node = [];
if P.data.n_nodes_expected >= 13
    P.smooth.low_pass_frq_by_node = 8 * ones(1, P.data.n_nodes_expected);
    P.smooth.low_pass_frq_by_node([9 10 13]) = 6;
    P.smooth.low_pass_frq_by_node([1 2]) = 8;
    P.smooth.low_pass_frq_by_node([5 6 7 8]) = 10;
end
P.smooth.break_at_repaired_points = true;
P.smooth.repair_break_halfwidth = 1;

%% QC
P.qc.enabled = true;
P.qc.frame_bad_node_frac_thresh = 0.25;
P.qc.max_interp_frac_per_frame = 0.35;
P.qc.max_jump_frac_per_frame = 0.35;
P.qc.max_geom_frac_per_frame = 0.35;
P.qc.max_arena_frac_per_frame = 0.35;
P.qc.max_lowconf_frac_per_frame = 0.50;
P.qc.expand_severe_badframes = 0;
P.qc.body_length_nodes = [1 10];
P.qc.body_length_prctile_range = [1 99];
P.qc.require_finite_body_length = true;
P.qc.geometry_node_pairs = [1 10; 2 10; 5 10; 6 10; 7 10; 8 10];
P.qc.geometry_prctile_range = [1 99];
P.qc.geometry_use_confidence_for_blame = true;
P.qc.arena_anchor_node = 10;
P.qc.arena_percentile = [1 99];
P.qc.arena_margin_px = 40;

%% Output control
P.output.return_raw = true;
P.output.keep_original_nan_mask = true;
P.output.store_full_masks = false;
P.output.make_plots = false;
P.output.audit_plot_dpi = 220;

%% Debug
P.debug.enabled = true;
P.debug.store_intermediate = false;
P.debug.verbose = false;
end
