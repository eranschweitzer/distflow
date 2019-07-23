function K = updateKmat(bus, branch, conn)
%%%
%%%  K = updateKmat(bus, branch, conn)
%%%  
%%%  Returns the K matrix used in the multiphase distflow.
%%%  This function is mainly needed if opt.mats_gen=1 is used
%%%  with distflow_multi, and the constant impedance loads have changed.
%%%
%%%  bus and branch are arrays, see distflow_multi
%%%  conn is a structure of matrices, returned from distflow_multi


ephasing = {branch.phase}.';
nphasing = {bus.phase}.';
Zconj = cellfun(@conj , {branch.Z},'UniformOutput', false);

if isfield(bus, 'yd') || isfield(bus, 'Ysh')
    [yl, ylc] = yload(bus);
    K = kdiag(Zconj, 'eye', ephasing)*conn.B*conn.TE*kdiag(yl, 'eye', nphasing(2:end)) + ...
        kdiag('eye', {branch.Z}, ephasing)*conn.B*conn.TE*kdiag('eye', ylc, nphasing(2:end));
else
    K = 0;
end
