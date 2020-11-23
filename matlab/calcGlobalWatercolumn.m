function [idx, z] = calcGlobalWatercolumn(lat, lon)

% lat: South is negative
% lon: West is negative
load('../TMs/MITgcm/Matrix5/Data/boxes.mat','Xbox','Ybox','Zbox');
load('../TMs/MITgcm/grid.mat','z');
% x(x>180) = x(x>180)-360;
% [tmp,idxX] = min(abs(lon-x)); 
% [tmp,idxY] = min(abs(lat-y)); 


Xbox(Xbox>180) = Xbox(Xbox>180)-360;
[val,idx1] = min(abs(sqrt((Xbox-lon).^2+(Ybox-lat).^2)));
XYbox = [Xbox,Ybox];

idx = find(ismember(XYbox,XYbox(idx1,:),'rows')==1);
z = Zbox(idx);