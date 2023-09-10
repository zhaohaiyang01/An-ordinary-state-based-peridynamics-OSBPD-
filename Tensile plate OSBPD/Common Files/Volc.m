function fac = Volc(idist,Geome)
delta = Geome.delta;
radij = Geome.radij;
if (idist <= delta - radij)
    fac = 1.0;
elseif (idist <= delta + radij)
    fac = (delta + radij - idist) / (2.0 * radij);
else
    fac = 0.0;
end