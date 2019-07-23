function [yl, ylc] = yload(bus)
%%%
%%%  [yl, ylc] = yload(bus)
%%%
%%%  generates the constant impedance load vectors that go into the K
%%%  matrix the multiphase distflow
%%%
%%%  bus is an nx1 structure array (see distflow_multi)

D = [1 -1 0; 0 1 -1; -1 0 1];
ydflag = isfield(bus,'yd');
yshflag = isfield(bus,'Ysh');
yyflag  = isfield(bus,'yy');
yl = cell(length(bus)-1, 1);
for k = 2:length(bus)
    %%% delta portion
    if ~ydflag
        yd = 0;
    elseif (length(bus(k).phase) < 2)
        if any(bus(k).yd ~= 0)
            warning('distflow_multi: Ignoring delta load on single phase bus.')
        end
        yd = 0;
    else
        yd = D(:,bus(k).phase)'*diag(conj(bus(k).yd))*D(:,bus(k).phase);
    end
    %%% shunt portion
    if ~yshflag
        ysh = 0;
    else
        ysh = bus(k).Ysh';
    end
    %% constant impedance laod
    if ~yyflag
        yy = 0;
    else
        yy = diag(bus(k).yy)';
    end
    yl{k-1} = (ysh + yy + yd).'; %note only transpose, NOT hermitian.
end
ylc = cellfun(@conj, yl, 'UniformOutput', false);
