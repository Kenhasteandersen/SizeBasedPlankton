% runWatercolumnAnnual;

%
% Implicit with various dt:
%
dt = [1, 0.25, 0.05, 0.01];
col = {'y','g','b','k'};
figure(3)
clf
for i = 1:length(dt)
    p.dt = dt(i);
    
    [sol, Bpico, Bnano, Bmicro] = analyse(@modelWatercolumnImplicit,p);
    
    figure(3)
    plot(sol.t, Bpico,col{i}, 'linewidth',1)
    hold on
    plot(sol.t, Bnano,col{i}, 'linewidth',2)
    plot(sol.t, Bmicro,col{i}, 'linewidth',3)
    drawnow
end
%%
% Compare implicit, implicitPC and full:
%
p.tEnd = 365*2;

p.dt = 0.1;
[sol, Bpico, Bnano, Bmicro] = analyse(@modelWatercolumnImplicit,p);
figure(3)
clf
plot(sol.t, Bpico,'g', 'linewidth',1)
hold on
plot(sol.t, Bnano,'g', 'linewidth',2)
plot(sol.t, Bmicro,'g', 'linewidth',3)
drawnow

[sol, Bpico, Bnano, Bmicro] = analyse(@modelWatercolumnImplicitPC,p);
figure(3)
plot(sol.t, Bpico,'b', 'linewidth',1)
hold on
plot(sol.t, Bnano,'b', 'linewidth',2)
plot(sol.t, Bmicro,'b', 'linewidth',3)
drawnow

[sol, Bpico, Bnano, Bmicro] = analyse(@modelWatercolumn,p);
figure(3)
plot(sol.t, Bpico,'k', 'linewidth',1)
hold on
plot(sol.t, Bnano,'k', 'linewidth',2)
plot(sol.t, Bmicro,'k', 'linewidth',3)
drawnow

xlim(365*[1 2])
%%
% Implicit with various dx and fixed dt:
%
nGrid = [10, 30, 100];
col = {'y','g','b','k'};
p.dt = 0.1;

figure(3)
clf
for i = 1:length(nGrid)
    p.nGrid = nGrid(i);
    
    [sol, Bpico, Bnano, Bmicro] = analyse(@modelWatercolumnImplicit,p);
    
    figure(3)
    plot(sol.t, Bpico,col{i}, 'linewidth',1)
    hold on
    plot(sol.t, Bnano,col{i}, 'linewidth',2)
    plot(sol.t, Bmicro,col{i}, 'linewidth',3)
    drawnow
    xlim(365*[1 2])
end
%%
% Implicit with various dx and adjusted dt
%
nGrid = [10, 30, 100];
col = {'y','g','b','k'};
p.dt = 0.1;

figure(3)
clf
for i = 1:length(nGrid)
    p.nGrid = nGrid(i);
    dx = max(p.xObs)/p.nGrid;
    p.dt = min(0.1, 3.75*dx^2/500);
    p.dt
    
    
    [sol, Bpico, Bnano, Bmicro] = analyse(@modelWatercolumnImplicit,p);
    
    figure(3)
    plot(sol.t, Bpico,col{i}, 'linewidth',1)
    hold on
    plot(sol.t, Bnano,col{i}, 'linewidth',2)
    plot(sol.t, Bmicro,col{i}, 'linewidth',3)
    drawnow
    xlim(365*[1 2])
end


%%
%  Full with different grids:
%
nGrid = [20, 30, 100];
col = {'y','g','b','k'};

figure(3)
clf
for i = 1:length(nGrid)
    p.nGrid = nGrid(i);
    [sol, Bpico, Bnano, Bmicro] = analyse(@modelWatercolumn,p);
    figure(3)
    plot(sol.t, Bpico,'k', 'linewidth',1)
    hold on
    plot(sol.t, Bnano,'k', 'linewidth',2)
    plot(sol.t, Bmicro,'k', 'linewidth',3)
    drawnow
    
    xlim(365*[1 2])
end

%%
%  Full with different tolerences:
%
tol = [1e-6 1e-3 1e-2 0.1];
col = {'y','g','b','k'};
p.nGrid = 20;

figure(3)
clf
for i = 1:length(tol)
    [sol, Bpico, Bnano, Bmicro] = analyse(@modelWatercolumn,p,tol(i));
    figure(3)
    plot(sol.t, Bpico,'k', 'linewidth',1)
    hold on
    plot(sol.t, Bnano,'k', 'linewidth',2)
    plot(sol.t, Bmicro,'k', 'linewidth',3)
    drawnow
    
    xlim(365*[1 2])
end


function [sol, Bpico, Bnano, Bmicro] = analyse(func,p,tol)
tic

if (nargin==2)
    [p,sol] = func(p);
else
    [p,sol] = func(p,tol);
end
toc

Bpico = sum(sol.Bpico'/p.nGrid);
Bnano = sum(sol.Bnano'/p.nGrid);
Bmicro = sum(sol.Bmeso'/p.nGrid);
end