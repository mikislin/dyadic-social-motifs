function f = local_fraction(num, den)
if den > 0
    f = num ./ den;
else
    f = NaN;
end
end