function [v, Pf, Qf] = distflow_lossy(mpc, opt)

%% input check
define_constants;
nb = size(mpc.bus,1);
nl = size(mpc.branch,1);
if ~all(mpc.bus(:,BUS_I) == (1:nb)')
    error('distflow_lossy: bus number must be consequtive starting at 1.')
end
if (nb - (nl + 1)) ~= 0
    error('distflow_lossy: case must be radial but nb=%d and nl=%d', nb, nl)
end
if mpc.bus(1,BUS_TYPE) ~= REF
    error('distflow_lossy: source bus must be first bus.')
end
%% options structure
if nargin < 2
    opt = [];
end
if isscalar(opt)
    opt = struct('alpha', opt, 'alpha_method', 1);
end
opt = optdefaults(opt);
%% operator setup

NumberOfNodes=size(mpc.bus,1);
FromNode=mpc.branch(:,F_BUS);
ToNode=mpc.branch(:,T_BUS);
Branch = 1:length(FromNode);
NumberOfBranch=length(Branch);
F = sparse(Branch,FromNode,1,NumberOfBranch,NumberOfNodes);
T = sparse(Branch,ToNode,1,NumberOfBranch,NumberOfNodes);
M0 = F-T;
M = M0(:,2:end);

r = mpc.branch(:,BR_R);
x = mpc.branch(:,BR_X);

Rline = sparse(1:NumberOfBranch, 1:NumberOfBranch,r);
Xline = sparse(1:NumberOfBranch, 1:NumberOfBranch,x);
TFT   = T*F';
I     = speye(NumberOfBranch);

slack_voltage = mpc.bus(mpc.bus(:,BUS_TYPE) == REF, VM);

pc = mpc.bus(:,PD)/mpc.baseMVA;
qc = mpc.bus(:,QD)/mpc.baseMVA;

B    = (I-TFT)^(-1);
C    = 2*(Rline*B*Rline + Xline*B*Xline) - (Rline*Rline + Xline*Xline);


%% setup alpha
mtd = opt.alpha_method;
if isscalar(opt.alpha)
    if mtd ~= 1
        warning('distflow_lossy: scalar alpha (%0.4f) was passed but alpha_method=%d.\n\tScalar alpha only makes sense with method 1. Input method selection ignored.',...
            opt.alpha, mtd)
    end
    alpha = opt.alpha*I;
%     tau  =(Rline*(B - I)*Rline+Xline*(B - I)*Xline)*(Rline*Rline+Xline*Xline)^(-1);
%     K    = (opt.alpha*I - (1-2*opt.alpha)*tau);
%     Kinv = K^(-1);
else
% methods:
%   2: linear interpolation based on branch impedance
%   3: quadradic interpolation based on branch impedance
%   4: linear interpolation based on downstream power
%   5: quadratic interpolation based on downstream power
%   6: linear interpolation based on branch impedance TIMES downstream
%      power
%   7: quadratic interpolation based on branch impedance TIMES downstream
%      power
    
    switch mtd
        case {2,3}
            v = sqrt(r.^2 + x.^2);
        case {4,5}
            v = abs(B*T*(pc + 1i*qc));
        case {6,7}
            v = sqrt(r.^2 + x.^2).*abs(B*T*(pc + 1i*qc));
        otherwise
            error('distflow_lossy: method %d not implemented.', mtd)
    end
    if ismember(mtd, [3,5,7])
        v = v.^2;
    end
    vmax = max(v); vmin = min(v);
    slope = (opt.alpha(2) - opt.alpha(1))/(vmin - vmax);
    intercept = opt.alpha(2) - slope*vmin;
            
    alpha = sparse(1:nl,1:nl,slope*v + intercept, nl, nl);
end
beta = I - C*(Rline*Rline+Xline*Xline)^(-1)*(I - 2*alpha);
%% solution
Rlossy = M\(beta\(Rline*B*T)*2);
Xlossy = M\(beta\(Xline*B*T)*2);

v = [slack_voltage^2 ;slack_voltage^2+Rlossy*pc+Xlossy*qc];

c  = (Rline*Rline+Xline*Xline)^(-1)*(I - 2*alpha)*M0*v;
v = sqrt(v);
Pf = mpc.baseMVA*B*(T*pc + Rline*c);
Qf = mpc.baseMVA*B*(T*qc + Xline*c);


function opt = optdefaults(opt)
optd = struct('alpha', 0.5, 'alpha_method', 1);
if isempty(opt)
    opt = optd;
else
    opt = struct_compare(optd, opt);
end

function b = struct_compare(a, b)
% compares structure b to structure a.
% if b lacks a field in a it is added
% this is performed recursively, so if if a.x is a structure
% and b has field x, the the function is called on (a.x, b.x)
for f = fieldnames(a).'
	if ~isfield(b, f{:})
		b.(f{:}) = a.(f{:});
	elseif isstruct(a.(f{:}))
		b.(f{:}) = struct_compare(a.(f{:}), b.(f{:}));
	end
end