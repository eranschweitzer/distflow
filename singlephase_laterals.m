function laterals = singlephase_laterals(mpc)

%% input check
define_constants;
nb = size(mpc.bus,1);
nl = size(mpc.branch,1);
if ~all(mpc.bus(:,BUS_I) == (1:nb)')
    error('distflow_lossy: bus number must be consequtive starting at 1.')
end
if (nb - (nl + 1)) ~= 0
    error('distflow_lossy: case must be radial but nb=%d and nl=%d', nb, nl)
end
if mpc.bus(1,BUS_TYPE) ~= REF
    error('distflow_lossy: source bus must be first bus.')
end

%% connection matrics
F = sparse(1:nl, mpc.branch(:,F_BUS), 1, nl, nb);
T = sparse(1:nl, mpc.branch(:,T_BUS), 1, nl, nb);
M = F - T; %incidence matrix
L = M'*M;  %laplacian
d = diag(L); % node degree
A = sparse(1:nb, 1:nb, d, nb, nb) - L; %adjacency matrix
%% from laterals
% laterals are defined by the number of leaf nodes.
% essentially a lateral is a path to a leaf.


leaves   = find(d == 1);
leaves   = leaves(leaves ~= 1); % remove source node
laterals = struct('nid', {cell(length(leaves), 1)}, 'bid', {cell(length(leaves), 1)});
for l = 1:length(leaves)
    p = dijkstra(A, 1, leaves(l)); % path to source
    pedge = find(ismember(mpc.branch(:,T_BUS), p(2:end)));
    laterals.nid{l} = ensure_col_vect(p);
    laterals.bid{l} = ensure_col_vect(pedge);
end
