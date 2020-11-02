function [p,sol] = modelWatercolumn(p, tol)

if (nargin==0) 
    p = parameters;
end

if (nargin==1)
    tol = 1e-3;
end

p.mortHTL = 0.005;

x = linspace(0,max(p.xObs),p.nGrid);
dx = x(2)-x(1);
p.x = x;
p.xmid = p.x(1:end-1) + diff(p.x)/2;

%
% Interpolate observations
%
if isfield(p,'obs')
    p.t = 0:p.tEnd; %(floor(min(p.tObs)+1)):(floor(max(p.tObs)));
    p.diff = interp2(p.tObs, p.xObs, p.diffObs, p.t, p.xmid');
    p.TObsInterp = interp2(p.tObs, p.xObs, p.TObs, p.t, p.xmid');
else
    t = 0:p.tEnd;
    p.t = t;
end
%
% Integrate:
%
options = odeset('abstol',tol, 'reltol',tol);
%options = odeset();

y = pdepe(0, @pdePDE, @pdeIC, @pdeBC, x, p.t, options);
%
% Extract results
% 
sol.y = y;
sol.N = squeeze(y(:,:,1));
sol.DOC = squeeze(y(:,:,2));
sol.B = y(:,:,3:end);
sol.B(sol.B<0)=0;

d = 10000 * 1.5 * (p.m*1e-6).^(1/3);
sol.Bpico = sum(sol.B(:,:,d<2),3);
sol.Bnano = sum(sol.B(:,:,d>=2 & d<20),3);
sol.Bmeso = sum(sol.B(:,:,d>=20),3);

sol.t = p.t;

plotWatercolumn(p,sol);




% --------------------------------------------------------------
    function [c,flux,source] = pdePDE(x,t,u,dudx)
        c = ones(2+p.n,1);
        %
        % Calculate diffusion and temperature from observations:
        %
        if isfield(p,'obs')
            ixTime = max(1,round(t));
            ixX = round((x-p.xmid(1))/dx+1);
            diff = p.diff(ixX, ixTime);
            p.T = p.TObsInterp(ixX, ixTime);
        else
            diff = p.diff;
        end
        
        flux = diff * dudx; % Only diffusion
        source = calcrates(t,x,u,p);
    end


% --------------------------------------------------------------
    function u0 = pdeIC(x)
        u0 = [p.N0; 0; ones(p.n,1)];
    end

% --------------------------------------------------------------
    function [pl,ql,pr,qr] = pdeBC(xl,ul,xr,ur,t)
        pl = zeros(2+p.n,1);
        ql = ones(2+p.n,1);
        pr = [ur(1)-p.N0; zeros(1+p.n,1)];%[0; 0];%
        qr = [0; ones(1+p.n,1)];%[1;1]; %
    end

end