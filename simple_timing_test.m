% Time test
clear variables;
% load('IEEE_123.mat')
load('IEEE_123_QSTS.mat')

opt = struct('alpha_method', 12, 'alpha', [0.498, 0.5], 'mats_gen', 1);
[Beta, K, zeta, eta, v0, conn] = distflow_multi(Bus, Branch,opt);
sigma0 = getsigma(Bus);
sigma  = sparse(size(sigma0,1), 8760);
for t = 1:8760
    sigma(:,t) = getsigma(Bus,t);
end

A = Beta*conn.M - K;
c = Beta*conn.M*v0;

% vcheck = sqrt(real(conn.U*(A\(c + zeta*sigma + eta*conj(sigma)))));

b0 = conn.U*(A\c);
b1 = conn.U*(A\zeta);
b2 = conn.U*(A\eta);

% b1 = conn.U*(btmp1 + btmp2);
% b2 = conn.U*(btmp1 - btmp2);

ttot = 0;
for t = 1:8760
  ttmp = tic;
  v = sqrt(real(b0 + b1*sigma(:,t) + b2*conj(sigma(:,t))));
  % v = sqrt(real(b0 + b1*real(sigma) + b2*imag(sigma)));
  ttot = ttot + toc(ttmp);
end
ttot/t