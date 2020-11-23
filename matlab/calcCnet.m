function Cnet = calcCnet(p,N,DOC,B,L,T)

setParameters(p);

u = [N, DOC, B];
Cnet = 0;
Cnet=calllib('model','calcCnet',Cnet,u,L,T);