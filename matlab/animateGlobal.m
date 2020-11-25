function animateGlobal(sim, sFilename)
if nargin==1
    sFilename = "Global";
end

figure(1)
clf
plotGlobalSimple(sim,1)
ax = gca;
ax.NextPlot = 'replaceChildren';
%%
n = length(sim.t);
F(n) = struct('cdata',[],'colormap',[]);

for i = 1:n
    plotGlobalSimple(sim,i)
    figure(1)
    drawnow
    F(i) = getframe(figure(1));
    disp(i);
end
%%
v = VideoWriter(sFilename, 'MPEG-4');
open(v);
for i = 1:n
    writeVideo(v,F(i));
end
close(v)
