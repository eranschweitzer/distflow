%%%% fit parameters to mean alpha data so that
%%%% alpha = x1 - x2*exp(x3*n)
clear variables
close all
%% load data
load('alphaest_opt.mat')
%%
y     = stats.mean;
x     = feeder_sizes.';
beta0 = [0.497, 0.002, -0.01]';
f     = @(x,b) b(1) - b(2)*exp(b(3)*x);
Jf    = @(x,b) [ones(length(x),1),...
                -exp(b(3)*x),...
                -b(2)*x.*exp(b(3)*x)];
R2    = @(y,f) 1 - sum((y-f).^2)/sum((y-mean(y)).^2);
Sold = norm(y - f(x,beta0),2);
b_old = beta0;
while true
    J = Jf(x, b_old);
    beta = b_old + (J'*J)\J'*(y-f(x,b_old));
    Snew = norm(y - f(x,beta),2);
    if abs(Sold - Snew) < 1e-8
        break
    end
    b_old = beta;
    Sold  = Snew;
end

%%
figure;
plot(feeder_sizes, stats.mean,'o', 'linewidth', 2)
hold on;
plot(feeder_sizes, f(feeder_sizes, beta), 'linewidth', 1.5)
xlabel('# of nodes')
ylabel('\alpha')
