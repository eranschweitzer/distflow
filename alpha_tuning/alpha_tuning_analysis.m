clear variables;
close all;
%%
saveout = true;
results = load('alpha_parametrization.mat', 'results');
results = results.results;

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
results.tau = taumap(results.tau);
[~,I] = sort(results.tau);
% reorder
results.tau  = results.tau(I);
for f = {'amin', 'amax', 'err'}
  results.(f{:}) = results.(f{:})(:,I);
  results.total.(f{:}) = results.total.(f{:})(I);
end

%% total results
Ttot = array2table([results.tau.' results.total.amin results.total.amax results.total.err], ...
  'VariableNames', {'tau', 'amin', 'amax', 'err'});

%% Comparison to variable ranges
% The total results give a certain range. We want to explore what
% is the extra error when this range is used in all the different
% feeder sizes given that proably it is not always exactly the optimal.

epsilon = collect_alpha_epsilon();
[amin, amax] = meshgrid(results.alpha_range);
idx = find( (amin==Ttot.amin(1)) & (amax==Ttot.amax(1)) );

cmperr = zeros(length(results.n), 1);
ptr = 1;
for n = results.n.'
	% the error we are seeking is at index idx for saved method 7 (1 in the paper)
	% which is in the last column:
	cmperr(ptr) = epsilon.(['fz' num2str(n)])(idx,end);
	ptr = ptr + 1;
end
% cmp err is the error using the range from the total results, min err is the minimum error for the 
% specific feeder size.
Tcmp = array2table([results.n, results.err(:,1), cmperr], 'VariableNames', {'fz', 'minerr', 'cmperr'})

%% save
if saveout
	writetable(Ttot, 'alpha_tuning_total_error_analysis.csv')
	writetable(Tcmp, 'alpha_tuning_per_size_comparison.csv')
end
