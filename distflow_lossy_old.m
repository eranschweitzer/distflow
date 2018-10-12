function [v, Pf, Qf] = distflow_lossy_old(mpc, alpha)

%% operator setup
define_constants;
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
tau  =(Rline*(B - I)*Rline+Xline*(B - I)*Xline)*(Rline*Rline+Xline*Xline)^(-1);
K    = (alpha*I - (1-2*alpha)*tau);
Kinv = K^(-1);

%% solution
Rlossy=(M\Kinv)*Rline*B*T;
Xlossy=(M\Kinv)*Xline*B*T;
v = [slack_voltage^2 ;slack_voltage^2+Rlossy*pc+Xlossy*qc];
v = sqrt(v);

Pf = mpc.baseMVA*B*(T*pc + ...
    (1-2*alpha)*Rline*(Rline*Rline + Xline*Xline)^(-1)*M0*v.^2);
Qf = mpc.baseMVA*B*(T*qc + ...
    (1-2*alpha)*Xline*(Rline*Rline + Xline*Xline)^(-1)*M0*v.^2); 