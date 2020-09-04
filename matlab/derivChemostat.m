function dudt = derivChemostat(t, y,p)
%y = [p.N0, p.DOC0, p.B0];
dudt = 0*y;

[tmp, tmp, tmp, tmp, tmp, dudt] = calllib('model','derivativeChemostat', ...
    p.L, p.T, p.d, p.N0, y, dudt);

end