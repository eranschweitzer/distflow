if false
clear variables;
close all;
%%
savenamebase = 'alpha_validation_mtd8and7';
alpha =  {[0.4830, 0.4990], 6.524};
alpha_method = {7, 8};
%%
define_constants;
nsamples = 5000;
mpopt  = mpoption('out.all',0, 'verbose', 0);
mpopt2 = mpopt;
mpopt2.pf.alg = 'ISUM';
mpopt2.pf.radial.max_it = 500;


%%
zvect = zeros(nsamples,3);
nmtds = length(alpha_method);
% err = struct('fz', zeros(nsamples,1),...
%       'lossy', struct('norm2', zvect, 'max', zvect, 'avg', zvect, 'std', zvect),...
%       'lossless', struct('norm2', zvect, 'max', zvect, 'avg', zvect, 'std', zvect));
fzvect     = zeros(nsamples,1);
norm2lossy = cellfun(@(x) zvect, cell(1,nmtds),'UniformOutput', false); 
maxlossy   = cellfun(@(x) zvect, cell(1,nmtds), 'UniformOutput', false); 
avglossy   = cellfun(@(x) zvect, cell(1,nmtds), 'UniformOutput', false); 
stdlossy   = cellfun(@(x) zvect, cell(1,nmtds), 'UniformOutput', false);
norm2lossless = zvect; maxlossless = zvect; avglossless = zvect; stdlossless = zvect;

%%
%parforstatus(nsamples, 0.1, 1);
%if isempty(gcp('nocreate'))
%    parpool(60);
%end
for k = 1:nsamples
  [n,e] = single_feeder_gen();
  fz    = length(n.id);
	if mod(k,50) == 0
		fprintf('Running sample %d (size %d).\n', k, fz)
	end
  %% matpower case
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
          fprintf('MATPOWER convergence failed: Feeder size %d, iter %d\n', fz, k)
          continue
      end
  end
  v  = r.bus(:, VM);
  pf = r.branch(:, PF);
  qf = r.branch(:, QF);
  
  fzvect(k) = fz;
	%% lossless distflow
  [vlossless, pflossless, qflossless] = distflow_lossy(r);
	%------ errros
  tmp = struct();
  tmp.lossless.v = v - vlossless;
  tmp.lossless.p = (pf - pflossless)/r.baseMVA;
  tmp.lossless.q = (qf - qflossless)/r.baseMVA;

  % 2 norm 
  norm2lossless(k,:) = [norm(tmp.lossless.v,2), norm(tmp.lossless.p,2), norm(tmp.lossless.q,2)]/sqrt(fz - 1);
  
  % max error
  maxlossless(k,:)   = [max(abs(tmp.lossless.v)), max(abs(tmp.lossless.p)), max(abs(tmp.lossless.q))];
  
  % avg error
  avglossless(k,:)   = [mean(abs(tmp.lossless.v)), mean(abs(tmp.lossless.p)), mean(abs(tmp.lossless.q))];
  
  % std error
  stdlossless(k,:)   = [std(abs(tmp.lossless.v)), std(abs(tmp.lossless.p)), std(abs(tmp.lossless.q))];

  %% lossy distflow
	for mtd = 1:nmtds
		distflowopt = struct('alpha', alpha{mtd}, 'alpha_method', alpha_method{mtd});
  	[vlossy, pflossy, qflossy] = distflow_lossy(r, distflowopt);
  
  	%% errors 
  	tmp = struct();
  	tmp.lossy.v = v - vlossy;
  	tmp.lossy.p = (pf - pflossy)/r.baseMVA;
  	tmp.lossy.q = (qf - qflossy)/r.baseMVA;

  	% 2 norm 
  	norm2lossy{mtd}(k,:)    = [norm(tmp.lossy.v,2), norm(tmp.lossy.p,2), norm(tmp.lossy.q,2)]/sqrt(fz - 1);
  	
  	% max error
  	maxlossy{mtd}(k,:)      = [max(abs(tmp.lossy.v)), max(abs(tmp.lossy.p)), max(abs(tmp.lossy.q))];
  	
  	% avg error
  	avglossy{mtd}(k,:)      = [mean(abs(tmp.lossy.v)), mean(abs(tmp.lossy.p)), mean(abs(tmp.lossy.q))];
  	
  	% std error
  	stdlossy{mtd}(k,:)      = [std(abs(tmp.lossy.v)), std(abs(tmp.lossy.p)), std(abs(tmp.lossy.q))];
	end
end
%delete(gcp('nocreate'))

%% place results in a structure
err = struct('fz', fzvect,...
      'lossless', struct('norm2', norm2lossless, 'max', maxlossless, 'avg', avglossless, 'std', stdlossless));
for mtd = 1:nmtds
	err.(['lossy_mtd' num2str(alpha_method{mtd})]) =  struct('norm2', norm2lossy{mtd}, 'max', maxlossy{mtd}, 'avg', avglossy{mtd}, 'std', stdlossy{mtd});
end

%% filter result and save
mask = err.fz ~= 0;
err.fz = err.fz(mask);
for field = {'norm2', 'max', 'avg', 'std'}
  err.lossless.(field{:}) = err.lossless.(field{:})(mask,:);
	for mtd = 1:nmtds
		mtdname = ['lossy_mtd' num2str(alpha_method{mtd})];
  	err.(mtdname).(field{:}) = err.(mtdname).(field{:})(mask,:);
	end
end
end
%%
fprintf('\n\nSaving error results.\n')
distflowopt = struct('alpha', alpha, 'alpha_method', alpha_method);
save([savenamebase '.mat'], 'err', 'distflowopt')

%% error histograms
fprintf('Calculating error histograms.\n')
histout = struct();
f = {'v', 'P', 'Q'};
for t = {'max', 'norm2'}
	maxv = max(err.lossless.(t{:}), [], 1);
	for mtd = 1:nmtds
		mtdname = ['lossy_mtd' num2str(alpha_method{mtd})];
		maxv = max(maxv, max(err.(mtdname).(t{:}), [], 1));
	end
	for k = 1:3
		dfstr = ['lossless_', t{:}, f{k}];
		if strcmp(t{:}, 'norm2') || strcmp(f{k}, 'v')
			dx = 0.01;
		else
			dx = 0.05;
		end
		edges = 0:dx:maxv(k)+dx;
		histout.(dfstr) = ensure_col_vect(histcounts(err.lossless.(t{:})(:,k), edges));
		for mtd = 1:nmtds
			mtdname = ['lossy_mtd' num2str(alpha_method{mtd})];
			ldfstr = [mtdname '_' t{:}, f{k}];
			histout.(ldfstr) = ensure_col_vect(histcounts(err.(mtdname).(t{:})(:,k), edges));
		end
	end
end
fprintf('Saving error histograms.\n')
save([savenamebase '_histograms.mat'], 'histout');

%% summary statistics
fprintf('Calculating summary statistics.\n')
summary_stats = struct('lossless', ...
	[mean(err.lossless.norm2, 1); mean(err.lossless.max, 1); max(err.lossless.max,[],1)]);
for mtd = 1:nmtds
	mtdname = ['lossy_mtd' num2str(alpha_method{mtd})];
	summary_stats.(mtdname) = [mean(err.(mtdname).norm2, 1); mean(err.(mtdname).max, 1); max(err.(mtdname).max,[],1)];
end
fprintf('Saving summary statistics.\n')
save([savenamebase '_summarystats.mat'], 'summary_stats')
