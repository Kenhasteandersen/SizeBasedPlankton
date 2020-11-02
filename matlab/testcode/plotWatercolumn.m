function plotWatercolumn(p,sol,time)

x = p.x;
t = p.t;

if (nargin==3)
    iTime = find(t>=time,1);
else
    iTime = length(t);
end
    


figure(1)
clf
subplot(1,2,1)
plot(sol.N(iTime,:)/10, -x,'b-',...
    sol.DOC(iTime,:), -x, 'm-', 'linewidth',3);
hold on
plot(sol.Bpico(iTime,:), -x, 'k-')
plot(sol.Bnano(iTime,:), -x, 'k-','linewidth',2)
plot(sol.Bmeso(iTime,:), -x, 'k-','linewidth',3)
ylabel('Depth (m)')
xlim([0 100])
%set(gca,'xscale','log')

subplot(1,2,2)
for i = 1:length(x)
    [dudt, rates(i)] = calcrates(t(iTime), x(i), squeeze(sol.y(iTime,i,:)), p, true);
    col(i,:,:) = 3*[rates(i).JF./p.m, rates(i).JLreal./p.m, rates(i).JDOC./p.m];
    %col(i,:,:) = [rates(i).JF, rates(i).JLreal, rates(i).JDOC];
    %col(i,:,:) = col(i,:,:)./max(squeeze(col(i,:,:))');
    %strategy(i,:) = calcstrategy(p,rates);
end

surface(log10(p.m),-x,0*(log10(squeeze(sol.B(iTime,:,:)))), col)
shading interp
hold on
contour(log10(p.m), -x, (log10(squeeze(sol.B(iTime,:,:)))), [-3:.5:2],'w'),
xlabel('log_{10} mass (\mu g_C)')



figure(2)
clf
subplot(5,1,1)
surface(t,-x,sol.N');
shading flat
hold on
plot3(t(iTime)*[1,1], ylim(),[100 100],'w:', 'linewidth',3)
axis("tight")
colorbar

subplot(5,1,2)
surface(t,-x,sol.DOC');
caxis([0 1])
shading flat
hold on
plot3(t(iTime)*[1,1], ylim(),[100 100],'w:', 'linewidth',3)
axis("tight")
colorbar

subplot(5,1,3)
surface(t,-x,sol.Bpico');
caxis([0 100])
shading flat
hold on
plot3(t(iTime)*[1,1], ylim(),[100 100],'w:', 'linewidth',3)
axis("tight")
colorbar

subplot(5,1,4)
surface(t,-x,sol.Bnano');
caxis([0 100])
shading flat
hold on
plot3(t(iTime)*[1,1], ylim(),[100 100],'w:', 'linewidth',3)
axis("tight")
colorbar

subplot(5,1,5)
surface(t,-x,sol.Bmeso');
caxis([0 100])
shading flat
hold on
plot3(t(iTime)*[1,1], ylim(),[100 100],'w:', 'linewidth',3)
axis("tight")
colorbar