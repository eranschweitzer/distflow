function A = kdiag(x,y, phases)
%%%
%%%   A = kdiag(x,y, phases)
%%%
%%%   Operator that forms a diagonal matrix of kroneger products.
%%%   x and y should be cells with matrix entries, 'eye', or 'gamma', 'gamma_conj'
%%%   phases is a cell array of the phases present each entry.
%%%   For use see the distflow_multi code.


n  = length(phases);
ridx = cell(n,1);
cidx = cell(n,1);
vidx = cell(n,1);
ptr  = 0;
for k = 1:n
    try
        xtmp = kdiag_tmpmat(x{k}, phases{k});
    catch ME
        if strcmp(ME.identifier, 'MATLAB:cellRefFromNonCell')
            xtmp = kdiag_tmpmat(x, phases{k});
        else
            rethrow(ME)
        end
    end
    try
        ytmp = kdiag_tmpmat(y{k}, phases{k});
    catch ME
        if strcmp(ME.identifier, 'MATLAB:cellRefFromNonCell')
            ytmp = kdiag_tmpmat(y, phases{k});
        else
            rethrow(ME)
        end
    end
    [rtmp,ctmp,vtmp] = find(kron(xtmp, ytmp));
	ridx{k} = ptr + rtmp;
	cidx{k} = ptr + ctmp;
	vidx{k} = vtmp;

	ptr = ptr + length(phases{k})^2;
end
A = sparse(cell2mat(ridx), cell2mat(cidx), cell2mat(vidx), ptr, ptr);

function xtmp = kdiag_tmpmat(x, phi)
if strcmp(x, 'eye')
    xtmp = eye(length(phi));
elseif strcmp(x, 'gamma')
    xtmp = gamma_phi(phi);
elseif strcmp(x, 'gammac')
    xtmp = conj(gamma_phi(phi));
else
    xtmp = x;
end
