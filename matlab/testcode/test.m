tspan = colon(0,0.01,200);
b = 0.083; % [offspring/ year/ individual] = [1/year] 10/70 ->this is only for females ->divide by 2 to include the whole population
d = 0.017; % [1/year] 1/70
r = b-d; % [1/year]
y0 = 100;
K = 1842; % [spermwhales]

k = @(t) (b*seasonal_repr(t)-d)/.000006;

opts = odeset('Maxstep', 1e-2);
%f = @(t,n) (b*seasonal_repr(t)-d)*n*(1-n/K);  % with constant K
 f = @(t,n) (b*seasonal_repr(t)-d)*n*(1-n/k(t)); % K as a function of t


[t,n] = ode45(f, tspan, y0, opts);


% Analytical solution
T = linspace(tspan(1), tspan(end), 100);
N = y0.*K.*exp(r.*T) ./ (K + y0.*(exp(r.*T)-1));

clf
plot(t, n, 'x')
hold on
plot(T,N)
xlabel('Time [years]')
ylabel('Number of spermwhales')
% legend('Seasonal Variation', 'Analytical Solution', 'Location', 'southeast')
grid on
hold off

function z=seasonal_repr(t)
period = 0.25; % if it is 1 whales reproduce all year long (constant birth rate)
% if period 0.25, whales reproduce 1/4 of the year -> so z=4 to keep birth
% rate at the value of 10/70/2

if(t-floor(t)< period)
    z = 1/period;
else
    z = 0;
end
end
