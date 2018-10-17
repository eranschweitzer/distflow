%%% parse alpha range sweep data and select best range and methods
clear variables;
close all;
%%
% for these files the error was divided by (fz-1) instead of 
% sqrt(fz-1), therefore the error needs to be renormalized, by 
% multiplying by sqrt(fz-1).
chng_error = {'alphaest_opt_with_mtds_10to200.mat',...
              'alphaest_opt_with_mtds_210to400.mat'};

fnames = {'alphaest_opt_with_mtds_10to200.mat',...
          'alphaest_opt_with_mtds_210to400.mat',...
          'alphaest_opt_with_mtds_410to450.mat',...
          'alphaest_opt_with_mtds_460to500.mat',...
          'alphaest_opt_with_mtds_510to550.mat',...
          'alphaest_opt_with_mtds_560to600.mat'};

%% collect data
epsilon = struct('total', 0);

for f = fnames
	fprintf('loading %s.\n', f{:})
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
	%%%% per size
	for k = 1:length(data.feeder_sizes)
		fz = data.feeder_sizes(k);
		tmp = sum(cat(3, data{k}{:}), 3);
		if ismember(f{:}, chng_error)
			tmp = tmp*sqrt(fz-1);
		end
		epsilon.(['fz' num2str(fz)]) = tmp;
	end
	
	%%%% total
	epsilon.total = epsilon.total + tmp;
end

%% calculate minimums
results.tau = tau;
results.alpha_range = alpha_range;
results.n = zeros(length(fieldnames(epsilon)), 1);
results.amin = zeros(length(results.n), length(tau));
results.amax = zeros(length(results.n), length(tau));
results.err    = zeros(length(results.n), length(tau));
[amin,amax] = meshgrid(alpha_range);

ptr = 1;
for f = fieldnames(epsilon).'
	[v, I] = min(epsilon.(f{:}), [], 1);
	alpha_min = amin(I);
	alpha_max = amax(I);
	if ~strcmp(f{:}, 'total')
		n = str2double(f{:}(3:end));

		results.n(ptr) = n;
		results.err(ptr,:) = v.';
		results.amin(ptr,:) = alpha_min.';
		results.amax(ptr,:) = alpha_max.';
		ptr = ptr + 1;
	else
		results.total = struct('err', v.', 'amin', alpha_min.', 'amax', alpha_max.');
	end	
end

%% sort by ascending size
[~,I]        = sort(results.n);
results.n    = results.n(I);
results.err  = results.err(I,:);
results.amin = results.amin(I,:);
results.amax = results.amax(I,:);

%% save
save('alpha_parametrization.mat', 'results');
