function sigma = getsigma(bus, t)
%%%
%%%  sigma = getsigma(bus)
%%%
%%%  Generates the sigma vector (essentially load power vector)
%%%  for the the multiphase distflow formulation
%%%
%%%  bus is an nx1 structure array (see distflow_multi)
%%%  t   (optional) is the time step (indexes into the profile structure)

if nargin < 2
	t = 0;
end
D = [1 -1 0; 0 1 -1; -1 0 1];
idx  = {1, [1,4].', [1,5,9].'};
sdflag = isfield(bus, 'sd');
ridx = cell(length(bus)-1,1);
vidx = cell(length(bus)-1,1);
ptr = 0;
for k = 2:length(bus)
    if sdflag && ~all(bus(k).sd == 0)
				profilesd = 1;
				if (t > 0) && isfield(bus(k), 'profile')
					if isfield(bus(k).profile, 'sd')
						profilesd = bus(k).profile.sd(t,:);
					end
				end
        tmpsd = ensure_col_vect(bus(k).sd).*profilesd.';
        sd = sqrt(3)*diag(D(:,bus(k).phase).'*diag(tmpsd)*D(:,bus(k).phase));
%         sd = diag(gamma_phi(bus(k).phase)*D(:,bus(k).phase).'*diag(tmpsd)*D(:,bus(k).phase));
    else
        sd = 0;
    end
        
    if any(bus(k).sy ~= 0) || any(sd~=0)
			  profilesy = 1;
				if (t > 0) && isfield(bus(k), 'profile')
					if isfield(bus(k).profile, 'sy')
						profilesy = bus(k).profile.sy(t,:);
					end
				end
        vidx{k} = ensure_col_vect(bus(k).sy).*profilesy.' + sd;
        ridx{k} = ptr + idx{length(bus(k).phase)};
        if length(vidx{k}) ~= length(ridx{k})
            error('getsigma: inconsistent sizes on bus %d between phase (%d x 1) and sy (%d x 1)', ...
                k, length(bus(k).phase), length(bus(k).sy))
        end
    end
    ptr = ptr + length(bus(k).phase)^2;
end
sigma = sparse(cell2mat(ridx), 1, cell2mat(vidx), ptr, 1);
