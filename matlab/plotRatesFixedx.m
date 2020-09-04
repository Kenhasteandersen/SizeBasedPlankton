function plotRatesFixedx(sol, p, iX, iB)

t = sol.t;

for j = 1:length(t)
    [dudt(j,:) r(j)] = calcrates(sol.t(j), sol.x(iX), [sol.N(j,iX), sol.DOC(j,iX), squeeze(sol.B(j,iX,:))']', p, true);
end

m = p.m(iB);
JF = [r.JF];
JL = [r.JL];
JLreal = [r.JLreal];
JN = [r.JN]*p.rhoCN;
JDOC = [r.JDOC];
Jtot = [r.Jtot];
JCloss = [r.JCloss];
JNloss = [r.JNloss]*p.rhoCN;
JR = [r.JR];

plot(t, Jtot(iB,:)/m, 'k-')
hold on
plot(t, JF(iB,:)/m,'r-')
plot(t, JL(iB,:)/m,'g--')
plot(t, JLreal(iB,:)/m,'g-')
plot(t, JN(iB,:)/m,'b-')
plot(t, JDOC(iB,:)/m,'m-')

plot(t, -JCloss(iB,:)/m,'m-')
plot(t, -JNloss(iB,:)/m,'b-')
plot(t, -JR(iB,:)/m,'r-')


plot(t,0*t,'k:')
hold off
