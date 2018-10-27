function errs = alphaestimation_opt(feeder_sizes, fname, mtds)
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
if nargin < 3
	mtds = 2:7;
	if (nargin < 2) || isempty(fname)
		fname = sprintf('alphaest_opt_with_mtds_%dto%d.mat', min(feeder_sizes), max(feeder_sizes));
	end
end
define_constants;
%feeder_sizes = 10:10:600; % number of nodes for samples
if ismember(8, mtds) && ~isscalar(mtds)
  error('method 8 works differently than the rest, it cannot be run along with other methods.')
end
if ismember(8,mtds)
  alpha_range = 0:0.001:40;
  nelem = length(alpha_range);
else
  alpha_range  = 0.4:0.001:0.6;
  [amin, amax] = meshgrid(alpha_range);
  nelem = numel(amin);
end
nsamples = 100;           % number of samples per feeder size

mpopt  = mpoption('out.all',0, 'verbose', 0);
mpopt2 = mpopt;
mpopt2.pf.alg = 'ISUM';
mpopt2.pf.radial.max_it = 500;
%% proccess samples
errs = cell(length(feeder_sizes), 1);
if isempty(gcp('nocreate'))
    parpool(min(nelem,60));
end
for k = 1:length(feeder_sizes)
    fz = feeder_sizes(k);
    fprintf('Running samples of size %d (%d of %d)...\n',fz,k,length(feeder_sizes))
    tmp = cell(nsamples,1);
		status = 0.025; status_step = 0.025;
    for iter = 1:nsamples
				if (iter/nsamples) >= status
					fprintf('\t%0.1f%% complete.\n', status*100)
					status = status + status_step;
				end
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
				vtrue = r.bus(2:end,VM);
				ptrue = r.branch(:,PF);
				qtrue = r.branch(:,QF);
%         err = zeros(length(alpha_range), 1);
        err = zeros(nelem, length(mtds));
				if ismember(8, mtds)
        	parfor kk = 1:nelem
        		opt = struct('alpha', alpha_range(kk), 'alpha_method', 8);
        	  [v, pf, qf] = distflow_lossy(r, opt);
        	  err(kk) = (norm(vtrue - v(2:end), 2) +...
        	        norm( (ptrue - pf)/r.baseMVA, 2) + ...
        	        norm( (qtrue - qf)/r.baseMVA, 2)) / sqrt(fz-1);
					end
				else
        	parfor kk = 1:nelem
%       	  for kk = 1:length(alpha_range)
							tmperr = zeros(1, length(mtds));
        	    for midx = 1:length(mtds)
        	        mtd = mtds(midx);
       	          opt = struct('alpha', [amin(kk), amax(kk)], 'alpha_method', mtd);
%       	          a = alpha_range(kk);
        	        [v, pf, qf] = distflow_lossy(r, opt);
        	        tmperr(midx) = (norm(vtrue - v(2:end), 2) +...
        	              norm( (ptrue - pf)/r.baseMVA, 2) + ...
        	              norm( (qtrue - qf)/r.baseMVA, 2)) / sqrt(fz-1);
        	    end
							err(kk,:) = tmperr;
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
save(fname, 'errs', 'feeder_sizes', 'nsamples', 'alpha_range', 'mtds', '-v7.3')
