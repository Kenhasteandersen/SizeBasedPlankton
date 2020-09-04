function [dudt, rates] = calcrates(t,x,u, p, bReturnRates)

if (nargin==4)
    bReturnRates = false;
end

%persistent N DOC B L ANmT JmaxT JR JL JN JF JNtot JLreal JCtot Jtot 
%persistent JCloss_feeding JCloss_photouptake JNlossLiebig JClossLiebig
%persistent JNloss JCloss mortpred dBdt dNdt dDOCdt  

N = u(1);%max(0,u(1));
DOC = u(2);%max(0,u(2));
B = u(3:end);
Bpositive = max(1e-8,B);

L = p.L*exp(-0.1*x)*(1-p.Lamplitude*cos(2*pi*t/365));  % Light
%
% Temperature corrections:
%
ANmT = p.ANm*1.5^(p.T/10-1);
JmaxT = p.Jmax*2^(p.T/10-1);
JR = p.Jresp*2^(p.T/10-1);
%
% Uptakes
%
if (N<0)
    JN = zeros(10,1);
else
    JN =   JmaxT/p.rhoCN .* ANmT*N ./ (JmaxT/p.rhoCN + ANmT*N); % Diffusive nutrient uptake
end
% in units of N/time
if (DOC<0)
    JDOC=zeros(10,1);
else
    JDOC = JmaxT .* ANmT*DOC ./ (JmaxT + ANmT*DOC); % Diffusive DOC uptake, units of C/time
    %JDOC = JmaxT .* (ANmT*DOC).^4 ./ (JmaxT.^4 + (ANmT*DOC).^4); % Diffusive DOC uptake, units of C/time
end

JL =   p.epsilonL * JmaxT .* p.ALm*L ./ (JmaxT + p.ALm*L);  % Photoharvesting

F = p.theta * Bpositive;
JF = p.epsilonF * JmaxT .* p.AFm.*F ./ (JmaxT + p.AFm.*F);        % Feeding

% Total nitrogen uptake:
JNtot = JN+JF/p.rhoCN; % In units of N

% Down-regulation of light uptake:
JLreal = min(JL, max(0, JNtot*p.rhoCN - (JDOC-JR)));

JCtot = JLreal+JF+JDOC-JR; % Total carbon untake

Jtot = min( JCtot, JNtot*p.rhoCN );  % Liebigs law; units of C
%
% Losses:
%
JCloss_feeding = (1-p.epsilonF)/p.epsilonF*JF; % Incomplete feeding (units of carbon per time)
JCloss_photouptake = (1-p.epsilonL)/p.epsilonL*JLreal;
JNlossLiebig = max(0, JNtot*p.rhoCN-JCtot)/p.rhoCN;  % N losses from Liebig
JClossLiebig = max(0, JCtot-JNtot*p.rhoCN); % C losses from Liebig, not counting losses from photoharvesting
%JClossLiebig = pmax(0, Jtot - JNtot*rhoCN) % C losses from Liebig, not counting losses from photoharvesting
%JClossLiebig = pmin(JClossLiebig, JDOC) % However, light surplus is not leaked but is downregulated

JNloss = JCloss_feeding/p.rhoCN + JNlossLiebig;
JCloss = JCloss_feeding + JCloss_photouptake + JClossLiebig;

%if (sum(c(JNloss,JCloss,B)<0))
%  browser()
%
% Mortality:
%
mortpred =  p.theta' * (JF/p.epsilonF.*Bpositive./p.m./F);
%
% System:
%
dBdt = (Jtot./p.m  - (mortpred + p.mort2*B + p.mortHTL.*(p.m>=p.mHTL))).*B;
%dBdt((B<1e-2) & (dBdt<0)) = 0; % Impose a minimum concentration even if it means loss of mass balance
mortloss = sum(B.*(p.mort2.*B + p.mortHTL.*(p.m>=p.mHTL)));
dNdt   =  - sum(JN.*B./p.m)   + sum(JNloss.*B./p.m) + p.remin*mortloss/p.rhoCN;
dDOCdt =  - sum(JDOC.*B./p.m) + sum(JCloss.*B./p.m) + p.remin*mortloss;


%if (floor(x)==30)
%    fprintf('%f %f %f %f: %f %f\n', DOC, - sum(JDOC.*B./p.m), sum(JCloss.*B./p.m), p.remin*mortloss, ...
%        DOC/-dDOCdt, dDOCdt);
%end


dudt = [dNdt; dDOCdt; dBdt];

if bReturnRates
   rates.JF = JF;
   rates.JN = JN;
   rates.JDOC = JDOC;
   rates.JL = JL;
   rates.JLreal = JLreal;
   rates.Jtot = Jtot;
   rates.JR = JR;
   rates.mortpred = mortpred;
   rates.JCloss = JCloss;
   rates.JNloss = JNloss;
end