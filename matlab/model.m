function model

p = parameters;

x = linspace(0,100,30);
t = 0:1000;

options = odeset();
tic
sol = pdepe(0, @pdePDE, @pdeIC, @pdeBC, x, t, options, p);
toc


% A surface plot is often a good way to study a solution.
figure(1)
clf
subplot(2,1,1)
surface(t,-x,u1')
ylabel('Depth')
shading flat

subplot(2,1,2)
surface(t,-x,log10(u2'))
xlabel('Time')
ylabel('Depth')
shading flat







% --------------------------------------------------------------
    function [c,flux,source] = pdePDE(x,t,u,dudx)
        c = [1; 1];
        flux = dudx; % Only diffusion
        source = calcrates(x,t,u);
    end


% --------------------------------------------------------------
    function u0 = pdeIC(x)
        u0 = [Nbottom+0*x; 0; ones(1,10)];
    end

% --------------------------------------------------------------
    function [pl,ql,pr,qr] = pdeBC(xl,ul,xr,ur,t)
        pl = [0; 0];
        ql = [1; 1];
        pr = [ur(1)-Nbottom; 0];%[0; 0];%
        qr = [0; 1];%[1;1]; %
    end

end