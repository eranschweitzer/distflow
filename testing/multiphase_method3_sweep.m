% exploring multiphase method 3: range of alpha linearlizy interpolated
% based on value of impedance.
%% load data
dssres = load('OpenDSSResults.mat');
load('Bus.mat');
load('Branch.mat');
opt = struct();
%% 
alpha_range = 0.47:0.001:0.5;
[amin, amax] = meshgrid(alpha_range);
%% run grid-search
res = cell(numel(amin),2);
for k = 1:numel(amin)
    opt.alpha = [amin(k), amax(k)];
    
    %% method 3
    opt.alpha_method = 3;
    r = distflow_multi(Bus, Branch, opt);
    res{k,1} = vertcat(r.vm);
    
    %% method 4 (quadratic)
    opt.alpha_method = 4;
    r = distflow_multi(Bus, Branch, opt);
    res{k,2} = vertcat(r.vm);
end

%% calculate error
err = cellfun(@(x) norm(dssres.x - x), res);