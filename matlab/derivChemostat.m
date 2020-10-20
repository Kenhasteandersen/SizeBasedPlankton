function dudt = derivChemostat(t, y,p)
dudt = 0*y;

[tmp, tmp, tmp, tmp, tmp, dudt] = calllib('model','derivativeChemostat', ...
    double(p.L), double(p.T), double(p.d), double(p.N0), y, dudt);

end