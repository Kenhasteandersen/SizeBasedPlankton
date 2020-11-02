p = parameters();

u = [p.N0, p.DOC0, p.B0]';
[dudt, rates] = calcrates(0,30,u,p,true);

plotrates(p,rates);