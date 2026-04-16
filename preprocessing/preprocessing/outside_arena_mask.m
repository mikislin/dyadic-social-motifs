function outside = outside_arena_mask(xy, arena)
%OUTSIDE_ARENA_MASK Flag samples outside estimated arena.

outside = false(size(xy,1),1);
if ~isfield(arena, 'isValid') || ~arena.isValid
    return
end
finiteMask = all(isfinite(xy),2);
outside(finiteMask) = xy(finiteMask,1) < arena.xlim(1) | xy(finiteMask,1) > arena.xlim(2) | ...
                    xy(finiteMask,2) < arena.ylim(1) | xy(finiteMask,2) > arena.ylim(2);
end
