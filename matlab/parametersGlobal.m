function p = parametersGlobal(n)
if (nargin==0)
    n = 10; % Minimal number for correct resolution of the size spectrum
end
% Basic model parameters:
p = parameters(n);

% Numerical parameters:
p.tEnd = 365; % In days
p.tSave = 365/12; % How often to save results (monthly)
p.dt = 0.1; % For Euler time stepping

% Light environment (??):
EinConv = 4.57; % conversion factor from W m^-2 to \mu mol s^-1 m^-2 (Thimijan & Heins 1983)
PARfrac = 0.4; % 
kw = 0.1; % 0.4 Camila ??

 load('../TMs/MITgcm/Matrix5/Data/boxes.mat','Ybox','Zbox');
 p.L = zeros(length(Ybox),730);
 for i = 1:730
     p.L(:,i) = EinConv*PARfrac*daily_insolation(0,Ybox,i/2,1).*exp(-kw*Zbox);
 end
 p.L(p.L<1) = 0;
 

