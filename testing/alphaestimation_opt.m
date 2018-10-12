clear variables
close all
%%
% methods:
%   2: linear interpolation based on branch impedance
%   3: quadradic interpolation based on branch impedance
%   4: linear interpolation based on downstream power
%   5: quadratic interpolation based on downstream power
%   6: linear interpolation based on branch impedance TIMES downstream
%      power
%   7: quadratic interpolation based on branch impedance TIMES downstream
%      power
define_constants;
feeder_sizes = 10:10:600; % number of nodes for samples
alpha_range  = 0.4:0.0001:0.6;
[amin, amax] = meshgrid(alpha_range);
mtds = 2:7;
nsamples = 100;           % number of samples per feeder size

mpopt  = mpoption('out.all',0, 'verbose', 0);
mpopt2 = mpopt;
mpopt2.pf.alg = 'ISUM';
mpopt2.pf.radial.max_it = 500;
%% proccess samples
errs = cell(length(feeder_sizes), 1);
if isempty(gcp('nocreate'))
    parpool(min(length(feeder_sizes),60));
end
parfor k = 1:length(feeder_sizes)
    fz = feeder_sizes(k);
    fprintf('Running samples of size %d (%d of %d)...\n',fz,k,length(feeder_sizes))
    tmp = cell(nsamples,1);
    for iter = 1:nsamples
        [n,e] = single_feeder_gen(fz);
        mpc = matpower_fmt(n,e,60);
        mpc = parallel_branch_join(mpc);
        % remove charging shunts which are not modeled in distflow
        mpc.branch(:,BR_B) = 0; 
        
        % set transformer tap to 1 since these are also not modled in
        % distflow
        mpc.branch(1,TAP) = 1;
        % to avoid extremely large voltage drops due to the ill-setup
        % transformer, reduce its impedance artificailly.
        mpc.branch(1,[BR_R,BR_X]) = mpc.branch(1,[BR_R,BR_X])/4;
        
        % solve matpower case
        r = runpf(mpc, mpopt);
        if ~r.success
            r = runpf(mpc, mpopt2);
            if ~r.success
                fprintf('MATPOWER convergence failed: Feeder size %d, iter %d\n', fz, iter)
                continue
            end
        end
%         err = zeros(length(alpha_range), 1);
        err = zeros(numel(amin), length(mtds));
        for kk = 1:length(err)
%         for kk = 1:length(alpha_range)
            for midx = 1:length(mtds)
                mtd = mtds(midx);
                opt = struct('alpha', [amin(kk), amax(kk)], 'alpha_method', mtd);
%                 a = alpha_range(kk);
                [v, pf, qf] = distflow_lossy(r, opt);
                err(kk, midx) = (norm(r.bus(2:end,VM) - v(2:end), 2) +...
                      norm( (r.branch(:,PF) - pf)/r.baseMVA, 2) + ...
                      norm( (r.branch(:,QF) - qf)/r.baseMVA, 2)) / (fz-1);
            end
        end
%         [ev, idx]  = min(err);
%         tmp{iter} = alpha_range(idx);
        tmp{iter} = err;
    end
    errs{k} = tmp;
end
delete(gcp('nocreate'))
%% alpha statistics
% stats.mean = cellfun(@mean, alpha);
% stats.std  = cellfun(@std,  alpha);
% stats.median = cellfun(@median, alpha);
%% save
save('alphaest_opt_with_mtds.mat', 'errs', 'feeder_sizes', 'nsamples', 'alpha_range', 'mtds', '-v7.3')