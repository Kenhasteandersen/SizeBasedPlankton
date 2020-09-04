plotWatercolumn(p,sol,1)
figure(1)
ax = gca;
ax.NextPlot = 'replaceChildren';

clear F
F(length(p.t/5)) = struct('cdata',[],'colormap',[]);

t = p.t(p.t>(p.t(end)-365));
t = t(1:5:end);
for j = 1:length(t)
    plotWatercolumn(p,sol,t(j))
    figure(1)
    drawnow
    F(j) = getframe(figure(1));
end

v = VideoWriter('annual','MPEG-4');
open(v);
for i = 1:length(F)
    writeVideo(v,F(i));
end
close(v)
