% downstream matrix testing

[n,e] = single_feeder_gen(); % generate random feeder

N = max(n.id);

%% option 1
F = sparse(1:N-1, e.f, 1, N-1, N);
T = sparse(1:N-1, e.t, 1, N-1, N);

B1 = (speye(N-1) - T*F.')\speye(N-1);

%% option 2
B2 = eye(N-1);
for k = N:-1:2
    tmp = n.pred(k);
    while tmp > 1
        B2(tmp - 1, k-1) = 1;
        tmp = n.pred(tmp);
    end
end

%%
test = B1 - B2;
fprintf('Difference between option1 and option2: %d\n', sum(test(:)));