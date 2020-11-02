function [p,sol] = modelWatercolumnImplicit(p)

if (nargin==0)
    p = parameters;
end
p.mortHTL = 0.005;

x = linspace(0, p.Depth, p.nGrid);
dx = x(2)-x(1);
p.x = x;
p.xmid = p.x(1:end-1) + dx/2;
dt = p.dt; 
nGrid = p.nGrid;

u00 = [p.N0, p.DOC0, p.B0];
for i = 1:length(x)
    u0(i,:) = u00;
end

p.v = 0;
p.vD = 1;
%
% Interpolate observations
%
if isfield(p,'obs')
    p.t = dt:dt:p.tEnd;
    xmidDiff = [p.x(1)-dx/2 p.xmid p.xmid(end)+dx/2];
    p.diff = interp2(p.tObs, p.xObs, p.diffObs, p.t, xmidDiff');
    p.TObsInterp = interp2(p.tObs, p.xObs, p.TObs, p.t, p.xmid');
else
    p.obs = false;
    t = 0:p.tEnd;
    p.t = t;
end
% -------------------
% Main loop:
% -------------------
u= zeros(length(p.t), nGrid, 12);
u(1,:,:)=u0;
%u(1,end/2:end,:)=0;

tic
for iTime = 2:p.tEnd/dt
    t = iTime*dt;
    %
    % Determine the diffusivity and surface light as a function of z at the time t:
    %
    if (p.obs)
        Dif = dt/(dx^2)*(p.diff(:,iTime));
    else
        Dif(1:nGrid+1) = dt/(dx^2)*p.diff;
    end
    
    nGrid = length(p.x);
    % Nutrient:
    a(1:nGrid) = -Dif(1:nGrid);
    b(1:nGrid) = 1 + Dif(1:nGrid)+Dif(2:(nGrid+1));
    c(1:nGrid) = -Dif(2:(nGrid+1));
    Snutrient(1:nGrid) = 0;
    % Nutrient BCs:
    b(1) = 1 + Dif(2); % Closed
    Snutrient(end) = -c(end)*p.N0; % Fixed concentration at the bottom
    % Assemble nutrient matrix:
    ANutrient = matrix(a,b,c); %
    
    % Phytoplankton:
    adv = p.v*dt/dx;
    a(1:nGrid) = -Dif(1:nGrid) - adv;
    b(1:nGrid) = 1 + adv + Dif(1:nGrid)+Dif(2:nGrid+1);
    c(1:nGrid) = -Dif(2:nGrid+1);
    s(1:nGrid) = 0;
    % Phytoplankton BCs:
    b(1) = 1 + adv + Dif(2); % Closed
    % Assemble phytoplankton matrix:
    APhyto = matrix(a,b,c);
    
    % Detritus:
    adv = p.vD*dt/dx;
    a(1:nGrid) = -Dif(1:nGrid) - adv;
    b(1:nGrid) = 1 + adv + Dif(1:nGrid)+Dif(2:nGrid+1);
    c(1:nGrid) = -Dif(2:nGrid+1);
    s(1:nGrid) = 0;
    % Detritus BCs:
    b(1) = 1 + adv + Dif(2); % Closed
    
    % Assemble matrix:
    ADetritus = matrix(a,b,c);
    % --------------------------------------
    % Calculate derivatives for the reactions:
    % --------------------------------------
    %
    % Growth rates:
    %
    for i = 1:length(p.x)
        dudt(i,:) = calcrates(t,x(i), squeeze(u(iTime-1,i,:)), p); %max(zeros(12,1), squeeze(u(iTime-1,i,:))), p);
    end
    %
    % Integrate advection and diffusion implicitly:
    % (predictor)
    u(iTime,:,1) = (ANutrient \ (dudt(:,1)*dt + u(iTime-1,:,1)' + Snutrient'));
    u(iTime,:,2) = (ANutrient \ (dudt(:,2)*dt + u(iTime-1,:,2)'));
    for i = 1:10
        u(iTime,:,2+i) = (APhyto \ (dudt(:,2+i)*dt + u(iTime-1,:,2+i)'));
    end
end
toc
%
% Extract results
%
sol.t = p.t;
sol.x = x;
sol.y = u;
sol.N = squeeze(u(:,:,1));
sol.DOC = squeeze(u(:,:,2));
sol.B = u(:,:,3:end);
sol.B(sol.B<0)=0;
%
d = 10000 * 1.5 * (p.m*1e-6).^(1/3);
sol.Bpico = sum(sol.B(:,:,d<2),3);
sol.Bnano = sum(sol.B(:,:,d>=2 & d<20),3);
sol.Bmeso = sum(sol.B(:,:,d>=20),3);

plotWatercolumn(p,sol);


% --------------------------------------
%  Function to calculate light:
% --------------------------------------
    function I = calcLight(P,I0)
        I = I0 * ...
            exp( -cumsum(param.kP*P)*param.dz ...  % Self-shading
            - param.k*z);                          % Extinction by water
    end

% --------------------------------------
%  Assemble the matrix from the three diagonals
% --------------------------------------
    function M = matrix(A,B,C)
        M = zeros(nGrid);
        M(1:nGrid+1:nGrid*nGrid) =B;
        M(2:nGrid+1:nGrid*nGrid) = A(2:end);
        M(nGrid+1:nGrid+1:nGrid*nGrid) = C(1:end-1);
        M = sparse(M);
    end



% --------------------------------------------------------------
    function [c,flux,source] = pdePDE(x,t,u,dudx)
        c = ones(2+p.n,1);
        %
        % Calculate diffusion and temperature from observations:
        %
        if isfield(p,'obs')
            ixTime = max(1,round(t));
            ixX = round((x-p.xmid(1))/dx+1);
            Dif = p.diff(ixX, ixTime);
            p.T = p.TObsInterp(ixX, ixTime);
        else
            Dif = p.diff;
        end
        
        flux = Dif * dudx; % Only diffusion
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