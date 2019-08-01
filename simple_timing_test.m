% Time test
clear variables;
load('IEEE_123.mat')

opt = struct('alpha_method', 12, 'alpha', [0.498, 0.5], 'mats_gen', 1);
[Beta, K, zeta, eta, v0, conn] = distflow_multi(Bus, Branch,opt);
sigma = getsigma(Bus);

A = Beta*conn.M - K;
c = Beta*conn.M*v0;

vcheck = sqrt(real(conn.U*(A\(c + zeta*sigma + eta*conj(sigma)))));

b0 = conn.U*(A\c);
b1 = conn.U*(A\zeta);
b2 = conn.U*(A\eta);

% b1 = conn.U*(btmp1 + btmp2);
% b2 = conn.U*(btmp1 - btmp2);

ttot = 0;
for t = 1:10000
  ttmp = tic;
  v = sqrt(real(b0 + b1*sigma + b2*conj(sigma)));
  % v = sqrt(real(b0 + b1*real(sigma) + b2*imag(sigma)));
  ttot = ttot + toc(ttmp);
end
ttot/t