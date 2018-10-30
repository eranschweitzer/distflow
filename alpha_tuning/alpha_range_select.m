%%% parse alpha range sweep data and select best range and methods
clear variables;
close all;
%%
% for these files the error was divided by (fz-1) instead of 
% sqrt(fz-1), therefore the error needs to be renormalized, by 
% multiplying by sqrt(fz-1).
chng_error = {'alphaest_opt_with_mtds_10to200.mat',...
              'alphaest_opt_with_mtds_210to400.mat'};

fnames = dir('alphaest_opt_with_mtds*.mat');

%% collect data
epsilon = struct('total', 0);
inds_tot = 0;
inds = struct();
for f = {fnames.name}
	fprintf('Loading %s.\n', f{:})
	data = load(f{:});
	if ~exist('alpha_range','var')
		alpha_range = data.alpha_range;
	else
		if ~all(alpha_range == data.alpha_range)
			error('alpha range for file %s different from already loaded range.', f{:})
		end
	end
	if ~exist('tau', 'var')
		tau = data.mtds;
	else
		if ~all(tau == data.mtds)
			error('methods in file %s different from already loaded set.', f{:})
		end
	end

	fprintf('\tCalculating per size and total errors\n')
	for k = 1:length(data.feeder_sizes)
		%%%% per size
		fz = data.feeder_sizes(k);
		if strcmp(f{:}, chng_error{1}) && (fz == 30)
			continue
		end
		inds.(['fz' num2str(fz)]) = sum(~cellfun(@isempty, data.errs{k}));
		tmp = sum(cat(3, data.errs{k}{:}), 3);
		if ismember(f{:}, chng_error)
			fprintf('\tCorrecting error normalization by multiplying by sqrt(%d-1)\n', fz)
			tmp = tmp*sqrt(fz-1);
		end
		epsilon.(['fz' num2str(fz)]) = tmp;

		%%%% total
		epsilon.total = epsilon.total + tmp;
		inds_tot = inds_tot + inds.(['fz' num2str(fz)]);
	end
	
end

%% calculate minimums
fprintf('Calculating minimum errors and corresponding ranges')
results.tau = [0 tau]; %add method 0 (lossless)
results.alpha_range = alpha_range;
results.n    = zeros(length(fieldnames(epsilon)) - 1, 1);
results.inds = zeros(length(fieldnames(epsilon)) - 1, 1);
results.amin = zeros(length(results.n), length(tau)+1);
results.amax = zeros(length(results.n), length(tau)+1);
results.err  = zeros(length(results.n), length(tau)+1);
[amin,amax] = meshgrid(alpha_range);
Idf    = find( (amin==0.5) & (amax == 0.5) ); %lossless index

ptr = 1;
for f = fieldnames(epsilon).'
	% find row vector of minimums, each entry correspnds to one of the methods
	% and the index is to the particular alpha_min alpha_max pair.
	[v, I] = min(epsilon.(f{:}), [], 1);
	alpha_min = amin(I);
	alpha_max = amax(I);
	if abs(sum(epsilon.(f{:})(Idf,:))/length(tau) - epsilon.(f{:})(Idf,1)) > 1e-8
		error('lossless errors do not match for the different methods (field=%s)', f{:})
	end
	vdf    = epsilon.(f{:})(Idf, 1);
	if ~strcmp(f{:}, 'total')
		n = str2double(f{:}(3:end));
		results.n(ptr) = n;
		results.inds(ptr) = inds.(f{:});
		results.err(ptr,:) = [vdf v];
		results.amin(ptr,:) = [0.5 alpha_min];
		results.amax(ptr,:) = [0.5 alpha_max];
		ptr = ptr + 1;
	else
		results.total = struct('err', [vdf; v.'], 'amin', [0.5; alpha_min.'], 'amax', [0.5; alpha_max.'], 'inds', inds_tot);
	end	
end

%% sort by ascending size
[~,I]        = sort(results.n);
results.n    = results.n(I);
results.inds = results.inds(I);
results.err  = results.err(I,:);
results.amin = results.amin(I,:);
results.amax = results.amax(I,:);

%% best solution
[v, tauind] = min(results.total.err);
Ibest  = find( (amin==results.total.amin(tauind)) & (amax==results.total.amax(tauind)) );

results.total.best = results.tau(tauind);

fprintf('\nCummulative Results\n')
fprintf('Best Error: method=%d, error=%0.3f, inds=%d, error/inds=%0.4f\n', results.tau(tauind), v, results.total.inds,  v/results.total.inds)
fprintf(' tau    Error   Inds   Error/Ind \n')
fprintf('----  --------  ----  ----------\n')
for k = 1:length(results.tau)
	fprintf('%4.1d  %8.3f  %4.1d  %10.7g\n', results.tau(k), results.total.err(k), results.total.inds, results.total.err(k)/results.total.inds)
end
fprintf('\n')
%% comparative results
% collect results per size for the best solution and lossless
if tauind == 1
	fprintf('Not collecting cmperr since best result is lossless found in results.err(:,1).\n')
	results.cmperr = results.err(:,1);
else
	results.cmperr = zeros(length(results.n), 1);
	ptr = 1;
	for n = results.n.'
		% the error we are seeking is at index Ibest, tauind-1 
		results.cmperr(ptr, :) = epsilon.(['fz' num2str(n)])(Ibest, tauind - 1);
		ptr = ptr + 1;
	end
end
fprintf('\nPer-Size Results\n')
fprintf(' fz   Inds   errmin (tau=%d)  errcmp    err df   d_errcmp  d_err df\n', results.tau(tauind))
fprintf('----  ----  ---------------  --------  --------  --------  --------\n')
for k = 1:length(results.n)
fprintf('%4.1d  %4.1d  %15.5f  %8.5f  %8.5f  %8.3f  %8.3f\n', results.n(k), results.inds(k), ...
	results.err(k,tauind), results.cmperr(k), results.err(k,1), 100*(results.cmperr(k)/results.err(k,tauind) - 1), 100*(results.err(k,1)/results.err(k,tauind) - 1))
end

%% save
save('alpha_parametrization.mat', 'results');

%% remap tau
% the methods in the distflow function that are used here are 2--7. In the
% paper they are 1--6, the mapping is:
%   saved     paper
%     7   -->  1
%     5   -->  2
%     3   -->  3
%     6   -->  4
%     4   -->  5
%     2   -->  6
taumap = [0 6 3 5 2 4 1];

papertau = [0 taumap(results.tau(2:end))].';
%% tables
Ttot = array2table([papertau results.total.amin results.total.amax results.total.err, results.total.err/results.total.inds], ...
  'VariableNames', {'tau', 'amin', 'amax', 'err', 'errperind'});

% cmp err is the error using the range from the total results, min err is the minimum error for the 
% specific feeder size. df is the lossless error.
Tcmp = array2table([results.n, results.inds, results.err(:,tauind), results.cmperr, results.err(:,1)], 'VariableNames', {'fz', 'inds', 'minerr', 'cmperr', 'df'});

writetable(Ttot, 'alpha_tuning_total_error_analysis.csv')
writetable(Tcmp, 'alpha_tuning_per_size_comparison.csv')
