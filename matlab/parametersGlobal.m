function p = parametersGlobal(n, pathTM)
%
% pathTM:
% 1 = MITgcm_2.8
% 2 = MITgcm_ECCO
%

if (nargin==0)
    n = 10; % Minimal number for correct resolution of the size spectrum
end
% Basic model parameters:
p = parameters(n);
%
% Set load paths
%
if (nargin==1 || nargin==0 || pathTM == 1)
    p.pathMatrix0   = '../TMs/MITgcm/Matrix5/TMs/matrix_nocorrection_0';
    p.pathMatrix1   =  '../TMs/MITgcm/Matrix5/TMs/matrix_nocorrection_';
    p.pathBoxes     = '../TMs/MITgcm/Matrix5/Data/boxes.mat';
    p.pathGrid      = '../TMs/MITgcm/grid.mat';
    p.pathConfigData = '../TMs/MITgcm/config_data.mat';
    p.pathTemp      = '../TMs/MITgcm/BiogeochemData/Theta_bc.mat'; 
    p.pathN0        = '../TMs/MITgcm_N0';
    p.pathInit      = '../TMs/globalInitMITgcm';

elseif pathTM == 2
    p.pathMatrix0 = '../TMs/MITgcm_ECCO/Matrix1/TMs/matrix_nocorrection_0';
    p.pathMatrix1 =  '../TMs/MITgcm_ECCO/Matrix1/TMs/matrix_nocorrection_';
    p.pathBoxes = '../TMs/MITgcm_ECCO/Matrix1/Data/boxes.mat';
    p.pathGrid = '../TMs/MITgcm_ECCO/grid.mat';
    p.pathConfigData = '../TMs/MITgcm_ECCO/config_data.mat';
    p.pathTemp = '../TMs/MITgcm_ECCO/BiogeochemData/Theta_bc.mat'; 
    p.pathN0    = '../TMs/MITgcm_ECCO_N0';
    p.pathInit      = '../TMs/globalInitMITgcm_ECCO';

end


% Numerical parameters:
p.tEnd = 365; % In days
p.tSave = 365/12; % How often to save results (monthly)
p.dt = 0.1; % For Euler time stepping
p.bParallel = true;
p.bTransport = true;

% Light environment (??):
p.EinConv = 4.57; % conversion factor from W m^-2 to \mu mol s^-1 m^-2 (Thimijan & Heins 1983)
p.PARfrac = 0.4; % ??
p.kw = 0.1; % 0.4 Camila ??
% 
% load('../TMs/MITgcm/Matrix5/Data/boxes.mat','Ybox','Zbox');
% p.L = zeros(length(Ybox),730);
% for i = 1:730
%     p.L(:,i) = EinConv*PARfrac*daily_insolation(0,Ybox,i/2,1).*exp(-kw*Zbox);
% end
% p.L(p.L<1) = 0;


