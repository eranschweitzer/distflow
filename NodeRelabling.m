f = [6, 1, 1, 3, 4];
t = [1, 3, 4, 2, 5];

F = sparse(1:5, f, 1, 5, 6);
T = sparse(1:5, t, 1, 5, 6);
M = F - T;

L = M'*M;
A = diag(diag(L)) - L;
G = A + eye(6);
%% initialize map
% nmap(i) -> new node id
nmap = zeros(6,1);
nmap(6) = 1;
ptr = 2;
%% BFS
x0 = [0, 0 , 0 , 0 , 0, 1].';
x = x0;
while ptr <= length(nmap)
    x = A'*x;
    new_nodes = find(x);
    for k = new_nodes.'
        if nmap(k) == 0
            nmap(k) = ptr;
            ptr = ptr + 1;
        end
    end
end
fnew = nmap(f);
tnew = nmap(t);

% x1 = logical(A'*x0);
% new_nodes = find(x1);
% for k = new_nodes.'
%     if nmap(k) == 0
%         nmap(k) = ptr;
%         ptr = ptr + 1;
%     end
% end
% 
% x2 = logical(A'*x1);
% new_nodes = find(x2);
% for k = new_nodes.'
%     if nmap(k) == 0
%         nmap(k) = ptr;
%         ptr = ptr + 1;
%     end
% end
% 
% x3 = logical(A'*x2);
% new_nodes = find(x3);
% for k = new_nodes.'
%     if nmap(k) == 0
%         nmap(k) = ptr;
%         ptr = ptr + 1;
%     end
% end