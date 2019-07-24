function v0 = v0update(bus, branch)
%%%
%%%			v0 = u0update(bus, branch)
%%% 
%%% Outputs the v0 vector for the multiphase distflow formulation
%%% for bus and branch see distflow_multi

ephasing = {branch.phase}.';
a = exp(-1i*2*pi/3);
avect = [1; a; a^2];
vref = bus(1).vref*avect;
vref = vec(vref*vref');

pnum = sqrt(length(vref));
idx  = reshape(1:length(vref), pnum, pnum);
cidx = cell2mat(cellfun(@(y) vec(idx(y,y)), ephasing, 'UniformOutput', false));
I0   = sparse(1:length(cidx), cidx, 1, length(cidx), numel(idx));
v0   = I0*vref;
