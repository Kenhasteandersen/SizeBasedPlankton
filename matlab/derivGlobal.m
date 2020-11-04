function dudt = derivGlobal(t,y,p)
%y = [p.N0, p.DOC0, p.B0];

dNdt = zeros(1,p.nb);
dDOCdt = zeros(1,p.nb);
dBdt = zeros(1,p.nb*p.n);
%Btmp = zeros(p.n,1);
nb = p.nb; % number of grid cells
n = p.n;    % number of plankton size classes


for i = 1:nb
    
    Ntmp    = y(i);
    DOCtmp  = y(nb+i);
    Btmp = y(i+2*nb:nb:end);
               
   
    %ICtmp = [Ntmp;DOCtmp;Btmp];
    dudtTmp = zeros(n+2,1);
    
    Ltmp = p.L(i); %current light level
    Ttmp = p.T(i);  %current temperature
    if any(isnan([Ntmp;DOCtmp;Btmp]))
        warning('NaNs before running current grid box');        
        keyboard
    end

    [tmp, tmp, tmp, tmp, tmp, dudtTmp] = calllib('model','derivativeChemostat', ...
        Ltmp, Ttmp,0,0, [Ntmp;DOCtmp;Btmp], dudtTmp);
    if any(isnan(dudtTmp))
        warning('NaNs after running current grid box');
        keyboard
    end
      
  
    dNdt(i) = dudtTmp(1);
    dDOCdt(i) = dudtTmp(2);
    dBdt(i:nb:end) = dudtTmp(3:end);
    
    
end
dudt = [dNdt,dDOCdt,dBdt]';

if ~isempty(find(isnan(dudt),1))
    keyboard
end
end