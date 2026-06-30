function badframe_fraction = local_badframe_fraction(sessionPreproc, stats)
if isstruct(stats) && isfield(stats, 'badframeFraction')
    badframe_fraction = stats.badframeFraction;
else
    badframe_fraction = mean(sessionPreproc.qc.badframes(:));
end
end