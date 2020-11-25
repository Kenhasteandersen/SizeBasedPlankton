function [idx, z] = calcGlobalWatercolumn(lat, lon, sim)

% lat: South is negative
% lon: West is negative
load(sim.p.pathBoxes,'Xbox','Ybox','Zbox');



Xbox(Xbox>180) = Xbox(Xbox>180)-360;
[val,idx1] = min(abs(sqrt((Xbox-lon).^2+(Ybox-lat).^2)));
XYbox = [Xbox,Ybox];

idx = find(ismember(XYbox,XYbox(idx1,:),'rows')==1);
z = Zbox(idx);