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
fname = {'IEEE_13.mat', 'IEEE_34.mat', 'IEEE_37.mat', 'IEEE_123.mat'};

%% 
alpha_range = 0.40:0.001:0.6;
[amin, amax] = meshgrid(alpha_range);
%% run grid-search
mtds = [11,12];
err = zeros(numel(amin),length(mtds));
if isempty(gcp('nocreate'))
    parpool(60)
end
parfor k = 1:numel(amin)
		if mod(k,100) == 0
			fprintf('running index %d\n', k)
		end
    errtmp = zeros(1,length(mtds));
    for f = 1:length(fname)
      data = load(fname{f});
      vnr   = vertcat(data.Bus.vm);
      sfnr  = cat(2,data.Branch.S).';
      for kk = 1:length(mtds)
        mtd = mtds(kk);

    %     	opt = struct('suppress_warnings', 1);
        opt = struct();
        opt.alpha = [amin(k), amax(k)];
        opt.alpha_method = mtd;
        [bnew, lnew] = distflow_multi(data.Bus, data.Branch, opt);

        v  = vertcat(bnew.vm);
        sf = vertcat(lnew.S);
        errtmp(kk) = (norm(vnr(4:end) - v(4:end),2) ...
          + norm(real(sfnr - sf),2) + norm(imag(sfnr - sf),2))/sqrt(length(sf));
  %       opt.bustmpopt = opt;
  %       opt.calcmu = 1;
  %       opt.alpha = 0.5;
  %       opt.alpha_method=1;
  %       opt.gamma_method=2;
  %       r = distflow_multi(Bus, Branch, opt);
  %       restmp{kk + length(mtds)} = vertcat(r.vm);
      end
      err(k,:) = err(k,:) + errtmp;
    end
end
delete(gcp('nocreate'))

%% save result
save('multiphase_parametrization_mtds11and12.mat', 'err', 'alpha_range', 'mtds');
%% calculate error
% err.norm = cellfun(@(x) norm(dssres.x - x), res);
% err.max  = cellfun(@(x) max(abs(dssres.x-x)), res);
% save('./testing/mutliphase_interpolation_results.mat', 'err', 'mtds', 'amin', 'amax')
