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
alpha_range = 0.40:0.001:0.6;
[amin, amax] = meshgrid(alpha_range);
%% run grid-search
mtds = 3:10;
res = cell(numel(amin),2*length(mtds));
if isempty(gcp('nocreate'))
    parpool(60)
end
parfor k = 1:numel(amin)
    
    restmp = cell(1,2*length(mtds));
    for kk = 1:length(mtds)
				mtd = mtds(kk);

    		opt = struct('suppress_warnings', 1);
    		opt.alpha = [amin(k), amax(k)];
        opt.alpha_method = mtd;
        r = distflow_multi(Bus, Branch, opt);
        restmp{kk} = vertcat(r.vm);
				
				opt.bustmpopt = opt;
				opt.calcmu = 1;
				opt.alpha = 0.5;
				opt.alpha_method=1;
				opt.gamma_method=2;
        r = distflow_multi(Bus, Branch, opt);
				restmp{kk + length(mtds)} = vertcat(r.vm);
    end
    res(k,:) = restmp;
end
delete(gcp('nocreate'))
warning(warnstate); %restore warnings
%% calculate error
err.norm = cellfun(@(x) norm(dssres.x - x), res);
err.max  = cellfun(@(x) max(abs(dssres.x-x)), res);
save('./testing/mutliphase_interpolation_results.mat', 'err', 'mtds', 'amin', 'amax')
