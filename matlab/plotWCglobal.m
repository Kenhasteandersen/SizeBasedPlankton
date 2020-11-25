function plotWCglobal(lat,lon,sim,iTime);
% lat: South is negative
% lon: West is negative
load(p.pathBoxes,'Xbox','Ybox','Zbox');


if (nargin()==3)
    iTime = length(sim.t);
end


Xbox(Xbox>180) = Xbox(Xbox>180)-360;

[val,idx1] = min(abs(sqrt((Xbox-lon).^2+(Ybox-lat).^2)));

XYbox = [Xbox,Ybox];

tmp = ismember(XYbox,XYbox(idx1,:),'rows');
idx = find(tmp==1);


figure 
plot(sim.N(idx,iTime),Zbox(idx))
axis ij
title(['Nitrogen, lat ', num2str(lat),', lon ', num2str(lon)])
ylabel('Depth (m)')
xlabel('Concentration (\mugN l^{-1})')

figure 
plot(sim.DOC(idx,iTime),Zbox(idx))
axis ij
title(['DOC, lat ', num2str(lat),', lon ', num2str(lon)])
ylabel('Depth (m)')
xlabel('Concentration (\mugC l^{-1})')

figure 
plot(sim.B(idx,:,iTime),Zbox(idx))
axis ij
title(['Nitrogen, lat ', num2str(lat),', lon ', num2str(lon)])
ylabel('Depth (m)')
xlabel('Concentration (\mugN l^{-1})')

end