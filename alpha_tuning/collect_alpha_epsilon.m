function epsilon = collect_alpha_epsilon()
%%% wrapper function for collecting all the errors from data
%%% this is coppied from the alpha_range_select script.

% for these files the error was divided by (fz-1) instead of 
% sqrt(fz-1), therefore the error needs to be renormalized, by 
% multiplying by sqrt(fz-1).
chng_error = {'alphaest_opt_with_mtds_10to200.mat',...
              'alphaest_opt_with_mtds_210to400.mat'};

fnames = dir('alphaest_opt_with_mtds*.mat');

%% collect data
epsilon = struct('total', 0);

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
		tmp = sum(cat(3, data.errs{k}{:}), 3);
		if ismember(f{:}, chng_error)
			fprintf('\tCorrecting error normalization by multiplying by sqrt(%d-1)\n', fz)
			tmp = tmp*sqrt(fz-1);
		end
		epsilon.(['fz' num2str(fz)]) = tmp;

		%%%% total
		epsilon.total = epsilon.total + tmp;
	end
	
end
