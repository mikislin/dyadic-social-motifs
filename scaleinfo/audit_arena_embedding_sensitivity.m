function T = audit_arena_embedding_sensitivity(ScaleScore, featureMeta, varargin)
%AUDIT_ARENA_EMBEDDING_SENSITIVITY Post hoc arena/feature-family audit.
%
% Arena labels are used only after condition-blind fitting to quantify
% possible domain shifts. They are not used for scale scoring or selection.

p = inputParser;
p.addParameter('nEmbeddingPCs', 4, @(x)isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('topNFeatures', 8, @(x)isnumeric(x) && isscalar(x) && x >= 1);
p.parse(varargin{:});
P = p.Results;

T = table();
if ~isstruct(ScaleScore) || ~isfield(ScaleScore, 'embeddingByScale')
    return
end

featureNames = string(featureMeta.Name(:));
featureFamily = string(featureMeta.Family(:));

for s = 1:numel(ScaleScore.embeddingByScale)
    E = ScaleScore.embeddingByScale{s};
    dimMeta = E.dimMeta;
    if isempty(dimMeta) || ~ismember('base_feature', dimMeta.Properties.VariableNames)
        continue
    end
    base = string(dimMeta.base_feature);
    family = repmat("unknown", numel(base), 1);
    [tf, loc] = ismember(base, featureNames);
    family(tf) = featureFamily(loc(tf));

    coeff = double(E.coeff);
    if isempty(coeff)
        continue
    end
    energy = sum(coeff.^2, 2);
    totalEnergy = sum(energy, 'omitnan');
    if ~(isfinite(totalEnergy) && totalEnergy > 0)
        totalEnergy = 1;
    end

    distanceFraction = sum(energy(family == "distance"), 'omitnan') ./ totalEnergy;
    contactFraction = sum(energy(family == "contact"), 'omitnan') ./ totalEnergy;
    kinematicsFraction = sum(energy(family == "kinematics"), 'omitnan') ./ totalEnergy;
    orientationFraction = sum(energy(family == "orientation"), 'omitnan') ./ totalEnergy;

    featureLevels = unique(base, 'stable');
    featureEnergy = nan(numel(featureLevels), 1);
    for f = 1:numel(featureLevels)
        featureEnergy(f) = sum(energy(base == featureLevels(f)), 'omitnan') ./ totalEnergy;
    end
    [featureEnergySorted, ord] = sort(featureEnergy, 'descend');
    topFeature = "";
    topFeatureFraction = NaN;
    topFeatureList = "";
    if ~isempty(ord)
        topFeature = featureLevels(ord(1));
        topFeatureFraction = featureEnergySorted(1);
        nTop = min(P.topNFeatures, numel(ord));
        parts = strings(nTop, 1);
        for i = 1:nTop
            parts(i) = featureLevels(ord(i)) + ":" + compose('%.4f', featureEnergySorted(i));
        end
        topFeatureList = strjoin(parts, ";");
    end

    arenaShift = NaN;
    arenaLevelsText = "";
    nArenaA = NaN;
    nArenaB = NaN;
    if istable(E.anchorTable) && ismember('arena_label', E.anchorTable.Properties.VariableNames)
        arenaAll = string(E.anchorTable.arena_label);
        labelOk = strlength(arenaAll) > 0 & ~ismissing(arenaAll);
        labelOk = labelOk(:) & (1:numel(arenaAll))' <= size(E.score, 1);
        arena = arenaAll(labelOk);
        levels = unique(arena, 'stable');
        arenaLevelsText = strjoin(levels, ";");
        if numel(levels) >= 2
            Xall = double(E.score(:, 1:min(P.nEmbeddingPCs, size(E.score, 2))));
            X = Xall(labelOk, :);
            a = arena == levels(1);
            b = arena == levels(2);
            nArenaA = nnz(a);
            nArenaB = nnz(b);
            d = nan(1, size(X, 2));
            for pc = 1:size(X, 2)
                xa = X(a, pc);
                xb = X(b, pc);
                pooled = iqr([xa(isfinite(xa)); xb(isfinite(xb))]);
                if ~(isfinite(pooled) && pooled > 0)
                    pooled = 1;
                end
                d(pc) = (median(xb, 'omitnan') - median(xa, 'omitnan')) ./ pooled;
            end
            arenaShift = sqrt(sum(d.^2, 'omitnan'));
        end
    end

    one = table();
    one.scale_index = local_scale_index(E, s);
    one.chunk_sec = E.chunk_sec;
    one.initial_band = local_scale_band(ScaleScore, s);
    one.n_embedding_pcs_audited = min(P.nEmbeddingPCs, size(E.score, 2));
    one.distance_family_loading_fraction = distanceFraction;
    one.contact_family_loading_fraction = contactFraction;
    one.kinematics_family_loading_fraction = kinematicsFraction;
    one.orientation_family_loading_fraction = orientationFraction;
    one.top_base_feature_by_loading = topFeature;
    one.top_base_feature_loading_fraction = topFeatureFraction;
    one.top_feature_loading_fractions = topFeatureList;
    one.arena_levels = arenaLevelsText;
    one.n_arena_level_a = nArenaA;
    one.n_arena_level_b = nArenaB;
    one.embedding_arena_shift_iqr_units = arenaShift;
    one.audit_only_not_selection = true;
    one.labels_used_for_embedding_fit = "none";
    one.arena_used_for_embedding_fit = false;
    one.arena_used_for_posthoc_audit = true;
    one.condition_used_for_embedding_fit = false;
    T = [T; one]; %#ok<AGROW>
end
end

function idx = local_scale_index(E, fallback)
idx = fallback;
if isfield(E, 'dimensionAudit') && istable(E.dimensionAudit) && ...
        ismember('scale_index', E.dimensionAudit.Properties.VariableNames)
    idx = double(E.dimensionAudit.scale_index(1));
end
end

function band = local_scale_band(ScaleScore, s)
band = "";
if isfield(ScaleScore, 'scaleTable') && istable(ScaleScore.scaleTable) && ...
        ismember('initial_band', ScaleScore.scaleTable.Properties.VariableNames)
    band = string(ScaleScore.scaleTable.initial_band(s));
end
end
