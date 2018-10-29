%%% parse data and select the best slope m in method 8
clear variables;
close all;
%%
data = load('alphaest_opt_mtd8.mat');

%% setup structures
epsilon = struct('total', 0);
results.n           = zeros(length(data.feeder_sizes),1);
results.m           = zeros(length(data.feeder_sizes),1);
results.err         = zeros(length(data.feeder_sizes),1);
results.inds        = zeros(length(data.feeder_sizes),1);
inds_tot = 0;
%% loop over feeders sizes
for k = 1:length(data.feeder_sizes)
	fz = data.feeder_sizes(k);
	results.inds(k) = sum(~cellfun(@isempty, data.errs{k}));
	inds_tot = inds_tot + results.inds(k);
	tmp = sum(cat(2, data.errs{k}{:}), 2);
	epsilon.(['fz' num2str(fz)]) = tmp;

	[v, I] = min(tmp);
	results.n(k)   = fz;
	results.m(k)   = data.alpha_range(I);
	results.err(k) = v;

	%% total
	epsilon.total = epsilon.total + tmp;
end

%% cummulative resutls
[vtot, Itot] = min(epsilon.total);
cumulative_results = struct('err', vtot, 'm', data.alpha_range(Itot), 'inds', inds_tot);

%% total to per-size comparison
results.cmperr = zeros(length(data.feeder_sizes), 1);
for k = 1:length(data.feeder_sizes)
	fz = data.feeder_sizes(k);
	results.cmperr(k) = epsilon.(['fz' num2str(fz)])(Itot);
end

%% save
save('alpha_mtd8_parametrization.mat', 'results', 'cumulative_results')

writetable(struct2table(results), 'alpha_tuning_mtd8_per_size_comparison.csv')
writetable(struct2table(cumulative_results), 'alpha_tuning_mtd8_total_error.csv')

%% print results
fprintf('Total Results\n')
fprintf('-----------------------------------\n')
fprintf('m value: %0.4f\n', cumulative_results.m)
fprintf('Error   Inds  Error/Ind\n')
fprintf('------  ----  ----------\n')
fprintf('%6.3f  %4.1d  %10.7g\n', cumulative_results.err, cumulative_results.inds, cumulative_results.err/cumulative_results.inds)
fprintf('\n\n')
fprintf('Per-Size Results\n')
fprintf('-----------------------------------\n')
fprintf('fz     m     Error   Inds  Error/Ind   cmperr    cmperr/Ind  cmp-min/cmp*100\n')
fprintf('---  -----   ------  ----  ----------  --------  ----------  --------------\n')
for k = 1:length(data.feeder_sizes)
	fprintf('%3.1d  %6.3f  %6.4f  %4.1d  %10.4f  %6.4f  %10.4f  %10.3f%%\n',...
 	results.n(k), results.m(k), results.err(k), results.inds(k), results.err(k)/results.inds(k),...
	results.cmperr(k), results.cmperr(k)/results.inds(k), 100*(results.cmperr(k)/results.err(k) - 1));
end
