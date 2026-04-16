function runs = find_runs(mask)
%FIND_RUNS Return Nx2 [start stop] runs for a logical mask.
mask = logical(mask(:));
d = diff([false; mask; false]);
starts = find(d == 1);
stops = find(d == -1) - 1;
runs = [starts stops];
end
