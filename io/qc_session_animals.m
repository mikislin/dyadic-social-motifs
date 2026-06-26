function qc = qc_session_animals(sessionRaw, varargin)
%QC_SESSION_ANIMALS Decide which tracked animal identities are usable.
%
% Some exported sessions contain an extra SLEAP tracks
% This function records raw/effective animal counts and the
% indices that downstream preprocessing should keep.

p = inputParser;
p.addParameter('minFiniteFrac', 0.05, @(x)isscalar(x) && x >= 0 && x <= 1);
p.addParameter('targetMaxAnimals', 2, @(x)isscalar(x) && x >= 1);
p.parse(varargin{:});
P = p.Results;

tracks = double(sessionRaw.SLEAPtracks);
if ndims(tracks) == 3
    tracks = reshape(tracks, size(tracks,1), size(tracks,2), size(tracks,3), 1);
end

nAnimalsRaw = size(tracks,4);
finiteFrac = nan(1,nAnimalsRaw);
for a = 1:nAnimalsRaw
    ta = tracks(:,:,:,a);
    finiteFrac(a) = mean(isfinite(ta(:)));
end

usable = finiteFrac >= P.minFiniteFrac;
keep = find(usable);
status = "ok";
notes = "";

if isempty(keep)
    status = "no_usable_animals";
    notes = "no_tracked_identity_passed_finite_fraction_threshold";
elseif numel(keep) > P.targetMaxAnimals
    status = "more_than_two_usable_animals";
    notes = "true_or_unresolved_multi_animal_session_not_reduced";
elseif nAnimalsRaw > numel(keep) && numel(keep) == P.targetMaxAnimals
    status = "reduce_to_dyad";
    notes = "extra_track_identity_removed_low_finite_fraction";
elseif nAnimalsRaw > numel(keep)
    status = "remove_low_finite_identity";
    notes = "low_finite_track_identity_removed";
end

qc = struct();
qc.n_animals_raw = nAnimalsRaw;
qc.n_animals_effective = numel(keep);
qc.keep_indices = keep(:)';
qc.finite_fraction_by_animal = finiteFrac;
qc.status = status;
qc.notes = notes;
qc.min_finite_fraction = P.minFiniteFrac;
end
