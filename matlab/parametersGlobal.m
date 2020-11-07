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
p.EinConv = 4.57; % ??
p.PARfrac = 0.4; % ??
p.kw = 0.1; % 0.4 Camila ??