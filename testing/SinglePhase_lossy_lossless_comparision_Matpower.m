% example on a Matlab feeder of comparison between lossy and lossless
% distflow
define_constants;

% optimal alpha/method setup
alpha = [0.4830, 0.4990];
alpha_method = 7;
mpopt  = mpoption('out.all',0, 'verbose', 0);
%%
casenames = {'case85'};
%% load matpower
mpc = loadcase(casenames{c});

% some modifcations
mpc.branch = mpc.branch(mpc.branch(:,BR_STATUS) > 0, :);
mpc.bus(:, [GS,BS]) = 0; % no shunt loads

%% solve matpower
r = runpf(mpc, mpopt);

 %% lossless distflow
rlossless = distflow_lossy(mpc);

%% lossy distflow
distflowopt = struct('alpha', alpha, 'alpha_method', alpha_method);
rlossy = distflow_lossy(mpc, distflowopt);

%% plot
figure;
singlephase_lateral_plot({r, rlossless, rlossy}, {'nr', 'df', 'ldf'});