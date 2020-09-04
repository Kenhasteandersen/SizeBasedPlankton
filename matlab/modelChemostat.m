p = parameters;
p.mortHTL= 0.05;
x = 10;


[t, u] = ode45(@deriv, [0 100], [p.N0, p.DOC0, p.B0],[],p);
%%
N = u(:,1);
DOC = u(:,2);
B = u(:,3:end);
B(B<0) = 0;

figure(1)
clf
subplot(2,1,1)
surface(t,p.m,log10(B'));
shading flat
set(gca,'yscale','log')
%zlim(log10([1 500]))
caxis(log10([1 500]))
colorbar

subplot(2,1,2)
[dudt, rates] = calcrates(t(end), x, [N(end) DOC(end) B(end,:)]', p, true);
plotrates(p, rates)

function dudt = deriv(t,u,p)
    dudt = calcrates(t, 10,u, p);
    dudt(1) = dudt(1) + p.diff*(p.N0-u(1));
    dudt(2) = dudt(2) - p.diff*u(2);
end


