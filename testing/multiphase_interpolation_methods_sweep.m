% exploring multiphase methods for approximating losses
% method #    Description
% --------    -----------
%    3        alpha interpolated linearly based on branch impedance
%    4        alpha interpolated quadratically based on branch impedance
%    5        alpha interpolated linearly based on branch downstream mva magnitude
%    6        alpha interpolated quadratically based on branch downstream mva magnitude
%    7        alpha (per phase) interpolated linearly based on branch impedance
%    8        alpha (per phase) interpolated quadratically based on branch impedance
%    9        alpha (per phase) interpolated linearly based on branch downstream mva magnitude
%   10        alpha (per phase) interpolated quadratically based on branch downstream mva magnitude
% based on value of impedance.
%% load data
clear variables;
close all;
dssres = load('OpenDSSResults.mat');
load('Bus.mat');
load('Branch.mat');
warnstate = warning('off','all');
%% 
alpha_range = 0.45:0.001:0.55;
[amin, amax] = meshgrid(alpha_range);
%% run grid-search
mtds = 3:10;
res = cell(numel(amin),length(mtds));
if isempty(gcp('nocreate'))
    parpool(2)
end
parfor k = 1:numel(amin)
    opt = struct();
    opt.alpha = [amin(k), amax(k)];
    
    restmp = cell(1,length(mtds));
    for kk = 1:length(mtds)
		mtd = mtds(kk);
        opt.alpha_method = mtd;
        r = distflow_multi(Bus, Branch, opt);
        restmp{kk} = vertcat(r.vm);
    end
    res(k,:) = restmp;
end
delete(gcp('nocreate'))
warning(warnstate); %restore warnings
%% calculate error
err = cellfun(@(x) norm(dssres.x - x), res);
save('./testing/mutliphase_interpolation_results.mat', 'err', 'mtds')
