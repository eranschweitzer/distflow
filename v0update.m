function v0 = v0update(bus, branch)
%%%
%%%			v0 = u0update
%%% phasing is a n-1 x 1 cell array where each entry contains a vector with 
%%% the phasing of the given branch. For example, if branch 7 (that one whose
%%% `to` node is 8) has phases A and B then: `phasing{7} = [1,2]`

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
