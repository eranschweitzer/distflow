function [p,w] = dijkstra(A,s,t)
%%% Dijkstra's shortest path algorithm based on Cormen, "Introduction to
%%% Algorithms" 3rd edition, pp. 651.
%%%
%%% INPUTS: A: adjacency matrix
%%%         s: source node
%%%         t: target node (nan for all targets)
%%% OUTPUTS: p: path from s to t or list of predecesors
%%%          w: weight of path from s to t or list of weights

    d    = inf(size(A,1),1);
    pi   = nan(size(A,1),1);
    d(s) = 0;
    S = [];
    Q = (1:size(A,1)).';
    while ~isempty(Q)
        [u,Q] = extract_min(Q,d);
        S = [S; u]; %#ok<AGROW>
        for v = find(A(u,:))
            if d(v) > d(u) + A(u,v)
               d(v)  = d(u) + A(u,v);
               pi(v) = u;
            end
        end
        if any(S == t) 
           break 
        end
    end
    if ~isnan(t)
        p = unravel_path(s,t,pi);
        w = d(t);
    else
        p = pi;
        w = d;
    end
    
    function [u,Q] = extract_min(Q,d)
        [~,I] = min(d(Q));
        u = Q(I);
        Q(I) = [];
    end
    
end