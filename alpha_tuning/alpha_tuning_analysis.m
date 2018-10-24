clear variables;
close all;
%%
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

%% 