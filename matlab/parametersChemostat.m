
function p = parametersChemostat(n)
if (nargin==0)
    n = 25;
end
p = parameters(n);
p.latitude=0;

p.d = 0.05;
