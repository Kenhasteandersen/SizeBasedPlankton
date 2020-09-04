%
% Run water column with gotm annual simulation
%

%
% Read gotm file
%
sFile = '../gotm/annual/annual.nc';

x = ncread(sFile, 'z');
x = -x(end:-1:1);
t = ncread(sFile, 'time')/(60*60*24);  % In days
diff = squeeze(ncread(sFile, 'nuh'))*(60*60*24);
T = squeeze(ncread(sFile, 'temp'));
%
% Use the two last years:
%
p = parameters;
p.obs = true;
ixTime = (t>(max(t)-2*365));
p.diffObs = double(diff(:,ixTime));
p.diffObs(end,:) = p.diffObs(end-1,:);
p.TObs = double(T(:, ixTime));
p.tObs = double(t(ixTime));
p.tObs = p.tObs-p.tObs(1);
p.xObs = double(x);
%
% Limit the diffusivity
%
p.diffObs(p.diffObs > 500) = 500;
%%
tic
p.nGrid = 30;
p.tEnd = 2*365;
%[p,sol] = modelWatercolumn(p,0.01);

%p.epsilonL = 1;
[p,sol] = modelWatercolumnImplicitPC(p);
%[p,sol] = modelWatercolumnImplicitPC(p);
toc