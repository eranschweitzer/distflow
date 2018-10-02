function bus = distflow_multi(bus, branch, alpha)
%%%    bus is an nx1 structure array with fields
%%%         - vref scalar reference voltage magnitude (p.u.) (source is
%%%         assumed to be 3 phase balanced)
%%%         - phase which a vector containing phases on bus.
%%%         - sy phi x 1 vector of per unit Y connected constant power load
%%%         (complex).
%%%         - yd 3x1 vector of delta connected constant admittance load
%%%         order is [ab, bc, ca]
%%%
%%%    branch is a n-1x1 structre array with fields
%%%         - R phase x phase matrix of branch resistance
%%%         - X phase x phase matrix of branch reactance
%%%         - phase a vector containing the phases of the branch.
%%%

a = exp(-1i*2*pi/3);
avect = [1; a; a^2];
idx = reshape(1:9,3,3);
if nargin < 3
    alpha = 0.5;
end
if nargin == 0
    %%% Demo
    %%% a1 ---- b1 ----- c1
    %%% a2 ---- b2 ----- c2 ----- e2
    %%% a3 ---- b3
    %%%    \--------d3
    
    zself = 0.01 +1i*0.1;
    zmut  = zself/2;
    zsamp = zmut*ones(3,3);
    zsamp(diag(idx)) = zself;
    nphasing = {[1,2,3], [1,2,3], [1,2], 3, 2}.';
    ephasing = nphasing(2:end);
    zarray   = cellfun(@(x) zsamp(x,x), ephasing,'UniformOutput', false);

    branch   = struct('f', {1, 2, 1, 3}.', 't', {2, 3, 4, 5}.', ...
    'Z', zarray, 'phase', ephasing);
    
    bus      = struct('phase', nphasing, ...
                      'sy', {[0 0 0].', [1 0.9, 1.1].', [0.25 0.25].', 0.1, 0.3}.',...
                      'yd', {[0 0 0].', 0 , [1/(1+1i*0.2), 0, 0].', 0, 0}.',...
                      'vref', 1.02);
    
else
    ephasing = {branch.phase}.';
    nphasing = {bus.phase}.';
end

%% form matrices
vref = bus(1).vref*avect;
conn = connmats(bus,branch);
sigma = getsigma(bus);
Zconj = cellfun(@conj , {branch.Z},'UniformOutput', false);
if ~isfield(branch, 'Y')
    Y    = cellfun(@(x) inv(x), {branch.Z}, 'UniformOutput', false);
else
    Y = {branch.Y};
end
Yconj= cellfun(@conj, Y, 'UniformOutput', false);
% [branch.Zc] = tmp{:};
v0   = v0vec(vec(vref*vref'), ephasing);
zeta = kdiag(Zconj, 'eye', ephasing)*conn.B*conn.TE*kdiag('eye','gamma',nphasing(2:end));
eta  = kdiag('eye', {branch.Z}, ephasing) *conn.B*conn.TE*kdiag('gammac','eye',nphasing(2:end));
if isfield(bus, 'yd')
    [yd, ydc] = ydelta(bus);
    K = kdiag(Zconj, 'eye', ephasing)*conn.B*conn.TE*kdiag(yd, 'eye', nphasing(2:end)) + ...
        kdiag('eye', {branch.Z}, ephasing)*conn.B*conn.TE*kdiag('eye', ydc, nphasing(2:end));
else
    K = 0;
end
Gamma = kdiag(Zconj, 'eye', ephasing)*(conn.B - conn.I)*kdiag(Yconj, 'eye', ephasing) + ...
        kdiag('eye', {branch.Z}, ephasing)*(conn.B - conn.I)*kdiag('eye', Y, ephasing);
Beta  = 2*alpha*conn.I - (1-2*alpha)*Gamma;

%% solve
nu = (Beta*conn.M - K)\(Beta*conn.M*v0 + zeta*sigma + eta*conj(sigma));
v2 = conn.U*nu;
if max(abs(imag(v2))) > 1e-8
    warning('distflow_multi: imaginary entries in v^2 with magnitude larger than 1e-8 found.\n\t These are discarded in the result')
end
v = sqrt(real(v2));
bus = updatebus(bus,v);


%% Utility functions
function S = connmats(bus, branch)
%%% f and t are vectors with bus indices, nphasing is a cell array with the
%%% phase vectors for **all the nodes** 

idx  = reshape(1:9,3,3);
ephasing = {branch.phase}.';
nphasing = {bus.phase}.';
f = [branch.f].';
t = [branch.t].';

% node mapping returning global index
ridx = cell2mat(cellfun(@(x,y) y*ones(length(x)^2,1), nphasing, num2cell(1:length(nphasing)).','UniformOutput',false));
cidx = cell2mat(cellfun(@(x,y) vec(idx(x,x)), nphasing, 'UniformOutput',false));
S.nidx = sparse(ridx,cidx,1:length(cidx),length(nphasing), numel(idx));

% edge mapping returning global index
ridx = cell2mat(cellfun(@(x,y) y*ones(length(x)^2,1), ephasing, num2cell(1:length(ephasing)).','UniformOutput',false));
cidx = cell2mat(cellfun(@(x,y) vec(idx(x,x)), ephasing, 'UniformOutput',false));
S.eidx = sparse(ridx, cidx, 1:length(cidx), length(ephasing), numel(idx));

S.E = max(S.eidx(end,:));
S.N = max(S.nidx(end,:));

cidx = cell2mat(cellfun(@(x,y) full(S.nidx(y,vec(idx(x,x)))).', nphasing(t), num2cell(f), 'UniformOutput', false));
S.F  = sparse(1:S.E, cidx,1, S.E, S.N);

cidx = cell2mat(cellfun(@(x,y) full(S.nidx(y,vec(idx(x,x)))).', nphasing(t), num2cell(t), 'UniformOutput', false));
S.T  = sparse(1:S.E, cidx, 1, S.E, S.N);

tmp  = min(S.nidx(2,:));
S.M  = S.F(:,tmp:end) - S.T(:,tmp:end);
S.TE = S.T(:,tmp:end);

S.I = sparse(1:S.E,1:S.E,1);

S.B = inv(S.I - S.T*S.F');

S.U = unvecd(nphasing(2:end));

function [yd, ydc] = ydelta(bus)
D = [1 -1 0; 0 1 -1; -1 0 1];
yd = cell(length(bus)-1, 1);
for k = 2:length(bus)
    if length(bus(k).phase) < 2
        if any(bus(k).yd ~= 0)
            warning('distflow_multi: Ignoring delta load on single phase bus.')
        end
        yd{k-1} = 0;
        continue
    end
    yd{k-1} = D(:,bus(k).phase)'*diag(conj(bus(k).yd))*D(:,bus(k).phase);
end
ydc = cellfun(@conj, yd, 'UniformOutput', false);

function sigma = getsigma(bus)

idx  = {1, [1,4].', [1,5,9].'};
ridx = cell(length(bus)-1,1);
vidx = cell(length(bus)-1,1);
ptr = 0;
for k = 2:length(bus)
    if ~all(bus(k).sy == 0)
        vidx{k} = bus(k).sy;
        ridx{k} = ptr + idx{length(bus(k).phase)};
    end
    ptr = ptr + length(bus(k).phase)^2;
end
sigma = sparse(cell2mat(ridx), 1, cell2mat(vidx), ptr, 1);

function g = gamma_phi(phases)
a = exp(-1i*2*pi/3);
gamma = [1  , a^2, a  ;
         a  , 1  , a^2;
         a^2, a  , 1];
g = gamma(phases,phases);

function [v0, I0] = v0vec(vref, ephasing)
%%% phasing is a n-1 x 1 cell array where each entry contains a vector with 
%%% the phasing of the given branch. For example, if branch 7 (that one whose
%%% `to` node is 8) has phases A and B then: `phasing{7} = [1,2]`

pnum = sqrt(length(vref));
idx  = reshape(1:length(vref), pnum, pnum);
cidx = cell2mat(cellfun(@(y) vec(idx(y,y)), ephasing, 'UniformOutput', false));
I0   = sparse(1:length(cidx), cidx, 1, length(cidx), numel(idx));
v0   = I0*vref;
% Iphi0 = speye(length(vref));
% x = cell2mat(cellfun(@(y) Iphi0(vec(idx(y,y)),:), phasing))*vref;

function A = kdiag(x,y, phases)
%%% x and y should be cells with matrix entries, 'eye', or 'gamma', 'gamma_conj'

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

function U = unvecd(phasing)
%%% for each block select the diagonal entries of the reshaped square
%%% matrix

n    = sum(cellfun(@length, phasing)); % number of phase variables
ridx = (1:n).';
cidx = cell(length(phasing),1);
idx  = {1, [1,4].', [1,5,9].'};
ptr = 0;
for k = 1:length(phasing)
    cidx{k} = ptr + idx{length(phasing{k})};
    ptr = ptr + length(phasing{k})^2;
end
U = sparse(ridx, cell2mat(cidx),1, n, ptr);

function bus = updatebus(bus,v)
ptr = 0;
for k = 2:length(bus)
    bus(k).vm = v(ptr + (1:length(bus(k).phase)));
    ptr = ptr + length(bus(k).phase);
end
