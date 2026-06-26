function manifest = build_session_manifest(sessionDir, dbasePath, varargin)
%BUILD_SESSION_MANIFEST Build the session manifest.
%
% The manifest is the stable bridge from raw/session files to analyses.
% It keeps Block 1 dyadic-motif inclusion separate from Block 2 egocentric
% context inclusion.

p = inputParser;
p.addParameter('writePath', '', @(x)ischar(x) || isstring(x));
p.addParameter('loadDbaseMetadata', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('verbose', true, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});
P = p.Results;

if nargin < 2
    dbasePath = '';
end

files = discover_session_files(sessionDir);
n = height(files);

meta = local_empty_metadata(n);
if P.loadDbaseMetadata && strlength(string(dbasePath)) > 0 && isfile(dbasePath)
    meta = local_read_dbase_metadata(dbasePath, n);
end

session_id = strings(n,1);
n_frames = nan(n,1);
n_nodes = nan(n,1);
n_coords = nan(n,1);
n_animals = nan(n,1);
effective_n_animals = nan(n,1);
animal_keep_indices = strings(n,1);
animal_finite_fraction = strings(n,1);
animal_qc_status = strings(n,1);
animal_qc_notes = strings(n,1);
has_scores = false(n,1);
score_dims = strings(n,1);
has_excluded_frames = false(n,1);
excluded_frame_count = nan(n,1);

for i = 1:n
    sessionPath = fullfile(sessionDir, files.session_file(i));
    S = load(sessionPath, 'sessionRaw');
    assert(isfield(S, 'sessionRaw'), 'build_session_manifest:MissingSessionRaw', ...
        'File does not contain sessionRaw: %s', files.session_path(i));
    raw = S.sessionRaw;
    assert(isfield(raw, 'SLEAPtracks'), 'build_session_manifest:MissingTracks', ...
        'sessionRaw.SLEAPtracks is required: %s', files.session_path(i));

    if isfield(raw, 'session_id') && ~isempty(raw.session_id)
        session_id(i) = string(raw.session_id);
    else
        session_id(i) = "session_" + compose('%04d', files.raw_index(i));
    end

    animalQC = qc_session_animals(raw);

    sz = size(raw.SLEAPtracks);
    n_frames(i) = sz(1);
    n_nodes(i) = sz(2);
    n_coords(i) = sz(3);
    if numel(sz) >= 4
        n_animals(i) = sz(4);
    else
        n_animals(i) = 1;
    end
    effective_n_animals(i) = animalQC.n_animals_effective;
    animal_keep_indices(i) = local_numeric_list(animalQC.keep_indices);
    animal_finite_fraction(i) = local_numeric_list(round(animalQC.finite_fraction_by_animal, 4));
    animal_qc_status(i) = string(animalQC.status);
    animal_qc_notes(i) = string(animalQC.notes);

    if isfield(raw, 'SLEAPscores') && ~isempty(raw.SLEAPscores)
        has_scores(i) = true;
        score_dims(i) = local_dims_string(size(raw.SLEAPscores));
    end

    if isfield(raw, 'excludedFrames') && ~isempty(raw.excludedFrames)
        has_excluded_frames(i) = true;
        excluded_frame_count(i) = nnz(logical(raw.excludedFrames(:)));
    else
        excluded_frame_count(i) = 0;
    end
end

mouse_codes = meta.mouse_codes;
mouse_type_1 = strings(n,1);
mouse_type_2 = strings(n,1);
mouse_type_3 = strings(n,1);
for i = 1:n
    codes = mouse_codes{i};
    if numel(codes) >= 1, mouse_type_1(i) = local_mouse_type(codes(1)); end
    if numel(codes) >= 2, mouse_type_2(i) = local_mouse_type(codes(2)); end
    if numel(codes) >= 3, mouse_type_3(i) = local_mouse_type(codes(3)); end
end

arena_label = strings(n,1);
for i = 1:n
    arena_label(i) = local_arena_label(meta.arena(i));
end

conditionInfo = local_parse_conditions(meta.condition);

isDyadic = effective_n_animals == 2;
isAnesthetizedPartner = conditionInfo.social_context == "anesthetized_partner" | ...
    mouse_type_1 == "ANES" | mouse_type_2 == "ANES" | mouse_type_3 == "ANES";
include_motif_discovery = isDyadic & ~isAnesthetizedPartner;
include_block1_dyadic = include_motif_discovery;
include_block2_egocentric = effective_n_animals == 1 | effective_n_animals == 2;
paper_include = include_motif_discovery | include_block2_egocentric;
exclusion_reason = strings(n,1);
exclusion_reason(effective_n_animals > 2) = "more_than_two_usable_animals_not_supported_current_pipeline";
exclusion_reason(effective_n_animals < 1) = "no_usable_animals";
exclusion_reason(~paper_include & exclusion_reason == "") = "unsupported_session_shape";
motif_exclusion_reason = strings(n,1);
motif_exclusion_reason(effective_n_animals ~= 2) = "not_two_usable_animals";
motif_exclusion_reason(isDyadic & isAnesthetizedPartner) = ...
    "anesthetized_center_special_egocentric_context";
motif_exclusion_reason(include_motif_discovery) = "";

manifest = table();
manifest.raw_index = files.raw_index;
manifest.session_file = files.session_file;
manifest.session_path = files.session_path;
manifest.session_rel_path = files.session_rel_path;
manifest.session_source_dir = files.session_source_dir;
manifest.session_id = session_id;
manifest.dbase_fileID = meta.fileID;
manifest.condition = meta.condition;
manifest.condition_primary = conditionInfo.primary;
manifest.social_context = conditionInfo.social_context;
manifest.drug = conditionInfo.drug;
manifest.dreadd_target = conditionInfo.dreadd_target;
manifest.dreadd_control_cohort = conditionInfo.dreadd_control_cohort;
manifest.arena = meta.arena;
manifest.arena_label = arena_label;
manifest.mouse_codes = strings(n,1);
for i = 1:n
    manifest.mouse_codes(i) = local_numeric_list(mouse_codes{i});
end
manifest.mouse_type_1 = mouse_type_1;
manifest.mouse_type_2 = mouse_type_2;
manifest.mouse_type_3 = mouse_type_3;
manifest.n_animals = n_animals;
manifest.effective_n_animals = effective_n_animals;
manifest.animal_keep_indices = animal_keep_indices;
manifest.animal_finite_fraction = animal_finite_fraction;
manifest.animal_qc_status = animal_qc_status;
manifest.animal_qc_notes = animal_qc_notes;
manifest.n_frames = n_frames;
manifest.n_nodes = n_nodes;
manifest.n_coords = n_coords;
manifest.has_scores = has_scores;
manifest.score_dims = score_dims;
manifest.has_excluded_frames = has_excluded_frames;
manifest.excluded_frame_count = excluded_frame_count;
manifest.file_bytes = files.file_bytes;
manifest.include_block1_dyadic = include_block1_dyadic;
manifest.include_motif_discovery = include_motif_discovery;
manifest.motif_exclusion_reason = motif_exclusion_reason;
manifest.include_block2_egocentric = include_block2_egocentric;
manifest.paper_include = paper_include;
manifest.exclusion_reason = exclusion_reason;
manifest.file_modified = files.file_modified;

manifest = apply_paper_cohort_definitions(manifest);

writePath = string(P.writePath);
if strlength(writePath) > 0
    outDir = fileparts(char(writePath));
    if ~isempty(outDir) && ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    writetable(manifest, writePath);
end

if P.verbose
    local_print_summary(manifest);
end
end

function meta = local_empty_metadata(n)
meta = struct();
meta.fileID = strings(n,1);
meta.arena = strings(n,1);
meta.condition = strings(n,1);
meta.mouse_codes = repmat({[]}, n, 1);
end

function meta = local_read_dbase_metadata(dbasePath, n)
meta = local_empty_metadata(n);
fields = {'fileID','arena','condition','mouse'};
raw = struct();
for i = 1:numel(fields)
    raw.(fields{i}) = h5read(dbasePath, ['/dbase/' fields{i}]);
end

nMeta = numel(raw.fileID);
nRead = min(n, nMeta);
for i = 1:nRead
    meta.fileID(i) = local_to_string(raw.fileID{i});
    meta.arena(i) = local_to_string(raw.arena{i});
    meta.condition(i) = local_to_string(raw.condition{i});
    meta.mouse_codes{i} = double(raw.mouse{i}(:)');
end
end

function s = local_to_string(x)
if isempty(x)
    s = "";
elseif isstring(x)
    s = x(1);
elseif ischar(x)
    s = string(x);
elseif isnumeric(x)
    vals = x(:)';
    if all(vals >= 0 & vals <= 65535 & vals == round(vals))
        s = string(char(vals));
    else
        s = string(strtrim(sprintf('%g ', vals)));
    end
else
    s = string(x);
end
s = strtrim(s);
end

function s = local_mouse_type(code)
if isempty(code) || ~isfinite(code)
    s = "";
elseif code == 1
    s = "WT";
elseif code == 2
    s = "MUT";
elseif code == 3
    s = "WT_DREADD_PFC";
elseif code == 4
    s = "WT_DREADD_DCN";
elseif code == 0
    s = "ANES";
else
    s = "code_" + string(code);
end
end

function s = local_arena_label(code)
code = upper(strtrim(string(code)));
if code == "B"
    s = "big";
elseif code == "S"
    s = "small";
else
    s = "";
end
end

function info = local_parse_conditions(condition)
n = numel(condition);
info = struct();
info.primary = strings(n,1);
info.social_context = strings(n,1);
info.drug = strings(n,1);
info.dreadd_target = strings(n,1);
info.dreadd_control_cohort = false(n,1);

for i = 1:n
    c = upper(strtrim(string(condition(i))));
    if startsWith(c, "WT")
        info.primary(i) = "WT";
    elseif startsWith(c, "M")
        info.primary(i) = "MUT";
    end

    if c == "WT"
        info.social_context(i) = "single";
    elseif endsWith(c, "_A")
        info.social_context(i) = "anesthetized_partner";
    elseif contains(c, "DCZ") || contains(c, "SAL")
        info.social_context(i) = "chemogenetic_context";
    else
        info.social_context(i) = "freely_moving_dyad";
    end

    if contains(c, "DCZ")
        info.drug(i) = "DCZ";
    elseif contains(c, "SAL")
        info.drug(i) = "SAL";
    else
        info.drug(i) = "none";
    end

    if contains(c, "PFC")
        info.dreadd_target(i) = "PFC";
    elseif contains(c, "DCN")
        info.dreadd_target(i) = "DCN";
    else
        info.dreadd_target(i) = "none";
    end

    info.dreadd_control_cohort(i) = contains(c, "DCZC") || contains(c, "SALC");
end
end

function s = local_numeric_list(x)
if isempty(x)
    s = "";
else
    s = string(strjoin(compose('%g', x), ';'));
end
end

function s = local_dims_string(sz)
s = string(strjoin(compose('%d', sz), 'x'));
end

function local_print_summary(manifest)
fprintf('Session manifest rows: %d\n', height(manifest));
fprintf('Block 1 dyadic sessions: %d\n', nnz(manifest.include_block1_dyadic));
fprintf('Block 2 egocentric-context sessions: %d\n', nnz(manifest.include_block2_egocentric));
fprintf('Excluded/special sessions: %d\n', nnz(~manifest.paper_include));
disp(groupsummary(manifest, 'n_animals'));
disp(groupsummary(manifest, 'effective_n_animals'));
disp(groupsummary(manifest, 'animal_qc_status'));
if all(strlength(manifest.condition) > 0)
    disp(groupsummary(manifest, {'arena','condition'}));
end
end
