function [nmap, rmap, fnew, tnew] = NodeRelabling(f, t, root, varargin)
%%% Relables radial nodes so that the `to` node is always greater than the
%%% `from` node.
%%%
%%%     [nmap, rmap, fnew, tnew] = NodeRelabling(f, t, root, varargin)
%%%
%%% INPUTS:
%%%     f: list of `from` nodes
%%%     t: list of `to` nodes
%%%     root: root nodes (source node)
%%%     optional inputs:
%%%         'plots': generates plots of the feeder with original and mapped
%%%         labels
%%%
%%% OUTPUTS:
%%%     nmap: mapping vector nmap(i) -> new node id
%%%     rmap: reverse mapping vector: rmap(new node id) -> old id
%%%     fnew: new `from` list (this is simply `nmap(f)`)
%%%     tnew: new `to` list (this is simply `nmap(t)`)

showplots = ismember('plots', varargin);
if nargin == 0
    % demo inputs
    f = [6, 1, 1, 3, 4].';
    t = [1, 3, 4, 2, 5].';
    root = 6;
    showplots = true;
elseif nargin < 3
    error('NodeRelabling: At least 3 arguments are required (unless running demo mode with 0 arguments)')
end
f = ensure_col_vect(f);
t = ensure_col_vect(t);
nl = length(f);
nb = nl + 1;
if ~all(sort(unique([f;t])) == (1:nb).')
    error('NodeRelabling: nodes must be labled consecutively.')
end

F = sparse(1:nl, f, 1, nl, nb);
T = sparse(1:nl, t, 1, nl, nb);
M = F - T;

L = M'*M;
A = diag(diag(L)) - L;
%% 
if showplots
    G = graph(A);
    figure;
    subplot(1,3,1)
    plot(G,'Layout', 'layered', 'Sources', root)
    title('Original Lables')
    set(gca,'xtick',[], 'xticklabel',[], 'ytick',[], 'yticklabel',[])
end
%% initialize map
% nmap(i) -> new node id
nmap = zeros(nb,1);
nmap(root) = 1;
ptr = 2;
%% BFS
x0 = zeros(nb,1);
x0(root) = 1;
while ptr <= length(nmap)
    x = (A'*x0) | x0;
    new_nodes = find(x - x0);
    for k = new_nodes.'    
        nmap(k) = ptr;
        ptr = ptr + 1;
    end
    x0 = x;
end
%% reverse map
rmap = full(sparse(nmap(1:nb),1,1:nb));
%% test rmap
if ~all(rmap(nmap(f)) == f)
    error('NodeRelabling: rmap test failed.')
end
%% mapped outputs
if nargout > 2
    fnew = nmap(f);
    tnew = nmap(t);
end
%%
if showplots

    subplot(1,3,2)
    plot(G,'Layout', 'layered', 'Sources', root, 'NodeLabel', nmap(1:nb))
    title('New Lables')
    set(gca,'xtick',[])
    set(gca,'xtick',[], 'xticklabel',[], 'ytick',[], 'yticklabel',[])
    
    labels = cell(nb,1);
    for k = 1:nb
        labels{k} = sprintf('%d<->%d',nmap(k), rmap(nmap(k)));
    end
    subplot(1,3,3)
    plot(G,'Layout', 'layered', 'Sources', root, 'NodeLabel', labels)
    title('Reverse Mapped Lables')
    set(gca,'xtick',[])
    set(gca,'xtick',[], 'xticklabel',[], 'ytick',[], 'yticklabel',[])
end