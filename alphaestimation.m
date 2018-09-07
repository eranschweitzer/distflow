clear variables
close all
%%
define_constants;
feeder_sizes = 10:10:500; % number of nodes for samples
nsamples = 100;           % number of samples per feeder size

mpopt = mpoption('out.all',0, 'verbose', 0);
%% proccess samples
alpha = cell(length(feeder_sizes), 1);
for k = 1:length(feeder_sizes)
    fz = feeder_sizes(k);
    fprintf('Running samples of size %d (%d of %d)...\n',fz,k,length(feeder_sizes))
    tmp = cell(nsamples,1);
    for iter = 1:nsamples
        [n,e] = single_feeder_gen(fz);
        mpc = matpower_fmt(n,e,60);
        mpc = parallel_branch_join(mpc);
        
        r = runpf(mpc, mpopt);
        if ~r.success
            continue
        end
        vi = r.bus(r.branch(:,F_BUS), VM);
        vj = r.bus(r.branch(:,T_BUS), VM);
        
        tmp{iter} = (vi.*vj - vj.^2)./(vi.^2 - vj.^2);
    end
    alpha{k} = vertcat(tmp{:});
end
%% alpha statistics
stats.mean = cellfun(@mean, alpha);
stats.std  = cellfun(@std,  alpha);

%% save
save('alphaest.mat', 'alpha', 'feeder_sizes', 'nsamples', 'stats', '-v7.3')

%% plot histograms
idx = [1, round(length(feeder_sizes)/2), length(feeder_sizes)];
figure;
for k = 1:3
    subplot(3,1,k)
    histogram(alpha{idx(k)},'Normalization', 'pdf')
    xlabel('\alpha')
    ylabel('Density')
    titlestr = sprintf('Feeders with %d Nodes', feeder_sizes(idx(k)));
    title(titlestr)
end