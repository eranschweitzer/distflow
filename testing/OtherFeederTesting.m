clear varaiables; close all;
%%
define_constants;
% optimal alpha/method setup
alpha = [0.4830, 0.4990];
% alpha = [0.499, 0.499];
% alpha = 0.499;
alpha_method = 7;
mpopt  = mpoption('out.all',0, 'verbose', 0);
%%
casenames = {'case18', 'case22', 'case33bw', 'case69', 'case85', 'case141'};
norm2lossy    = zeros(length(casenames), 3);
norm2lossless = zeros(length(casenames), 3);
maxlossy      = zeros(length(casenames), 3);
maxlossless   = zeros(length(casenames), 3);
for c = 1:length(casenames)
  %% load matpower
  mpc = loadcase(casenames{c});
  %% some modifcations
  mpc.branch = mpc.branch(mpc.branch(:,BR_STATUS) > 0, :);
  mpc.bus(:, [GS,BS]) = 0; % no shunt loads

  fz  = size(mpc.bus,1);
  nl  = size(mpc.branch, 1);
  if fz - nl - 1 ~= 0
    error('Not a radial system')
  end
  %% solve matpower
  r = runpf(mpc, mpopt);
  v  = r.bus(:, VM);
  pf = r.branch(:, PF);
  qf = r.branch(:, QF);
  %% renumber stuff
  f = mpc.branch(:,F_BUS);
  t = mpc.branch(:,T_BUS);
  root = mpc.bus(mpc.bus(:,BUS_TYPE)==REF, BUS_I);

  [nmap, rmap, fnew, tnew] = NodeRelabling(f, t, root);


  mpci = mpc;
  mpci.bus(:,BUS_I)    = nmap(mpc.bus(:,BUS_I));
  mpci.branch(:,F_BUS) = nmap(mpc.branch(:,F_BUS));
  mpci.branch(:,T_BUS) = nmap(mpc.branch(:,T_BUS));
  mpci.gen(:,GEN_BUS)   = nmap(mpc.gen(:,GEN_BUS));

  % sort mpci so bus number are consecutive
  [~,busmap] = sort(mpci.bus(:, BUS_I));
  rbusmap    = full(sparse(busmap ,1, 1:fz));
  mpci.bus = mpci.bus(busmap,:);
  %% lossless distflow
  [vlossless, pflossless, qflossless] = distflow_lossy(mpci);

  %% lossy distflow
  distflowopt = struct('alpha', alpha, 'alpha_method', alpha_method);
  [vlossy, pflossy, qflossy] = distflow_lossy(mpci, distflowopt);
  %% Order bus results like mpc
  vlossless = vlossless(rbusmap);
  vlossy    = vlossy(rbusmap);
  %% errors
  tmp = struct();
  tmp.lossy.v = v - vlossy;
  tmp.lossy.p = (pf - pflossy)/r.baseMVA;
  tmp.lossy.q = (qf - qflossy)/r.baseMVA;
  tmp.lossless.v = v - vlossless;
  tmp.lossless.p = (pf - pflossless)/r.baseMVA;
  tmp.lossless.q = (qf - qflossless)/r.baseMVA;

  % 2 norm 
  norm2lossy(c,:)    = [norm(tmp.lossy.v,2), norm(tmp.lossy.p,2), norm(tmp.lossy.q,2)]/sqrt(fz - 1);
  norm2lossless(c,:) = [norm(tmp.lossless.v,2), norm(tmp.lossless.p,2), norm(tmp.lossless.q,2)]/sqrt(fz - 1);

  % max error
  maxlossy(c,:)      = [max(abs(tmp.lossy.v)), max(abs(tmp.lossy.p)), max(abs(tmp.lossy.q))];
  maxlossless(c,:)   = [max(abs(tmp.lossless.v)), max(abs(tmp.lossless.p)), max(abs(tmp.lossless.q))];
end
%%
% figure;
% subplot(3,1,1)
% plot(1:fz, v, 1:fz, vlossless, 1:fz, vlossy)
% legend('nr', 'df', 'ldf')
% ylabel('voltage [p.u.]')
% xlabel('bus #')
% 
% subplot(3,1,2)
% plot(1:nl, pf, 1:nl, pflossless, 1:nl, pflossy)
% legend('nr', 'df', 'ldf')
% ylabel('Pf [MW]')
% 
% subplot(3,1,3)
% plot(1:nl, qf, 1:nl, qflossless, 1:nl, qflossy)
% legend('nr', 'df', 'ldf')
% ylabel('Qf [MVAr]')
% 
fprintf(['                 v       Pf        Qf   \n',...
         'mean rms df : %5.4f   %5.4f    %5.4f\n',...
         'mean rms ldf: %5.4f   %5.4f    %5.4f\n',...
         'max  int df : %5.4f   %5.4f    %5.4f\n',...
         'max  inf ldf: %5.4f   %5.4f    %5.4f\n'], ...
         mean(norm2lossless,1), mean(norm2lossy,1), max(maxlossless,[],1), max(maxlossy,[],1));
       
fprintf(['                 v       Pf        Qf   \n',...
         'mean rms df : %5.4f   %5.4f    %5.4f\n',...
         'mean rms ldf: %5.4f   %5.4f    %5.4f\n',...
         'mean int df : %5.4f   %5.4f    %5.4f\n',...
         'mean inf ldf: %5.4f   %5.4f    %5.4f\n'], ...
         mean(norm2lossless,1), mean(norm2lossy,1), mean(maxlossless,1), mean(maxlossy,1));