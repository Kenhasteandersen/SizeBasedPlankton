function saveGlobal(s)

sim.N = s.N(:,end);
sim.DOC = s.DOC(:,end);
sim.B = s.B(:,:,end);
sim.t = 0;
sim.p = s.p;

save('../TMs/globalInit.mat','sim');
