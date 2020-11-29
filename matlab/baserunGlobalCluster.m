c = parcluster('dcc R2019a');
c.AdditionalProperties.EmailAddress = 'kha@aqua.dtu.dk';
c.AdditionalProperties.MemUsage = '15GB';
c.AdditionalProperties.ProcsPerNode = 0;
c.AdditionalProperties.WallTime = '4:00';
c.saveProfile
%
clust = parcluster('dcc R2019a');
numW=10;    % Exactly the number of nodes times the number of processors per cores requested
parpool(clust, numW);

p = parametersGlobal(10,2);
load(p.pathInit);
p.tEnd = 5*365;
sim = simulateGlobal(p,sim);
saveGlobal(sim);
save('tmp','sim','-v7.3');
%baserunGlobal

%load ../TMs/globalInitMITgcm_ECCO.mat
%p = parametersGlobal(10,2);
%p.tEnd = 3*365;

%
% 80 cores: 11 minutes
% 20 cores: 4:01
% 10 cores: 4:18 min
% 5 cores:  5:24
