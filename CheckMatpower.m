% clc
clear variables
close all

%% MAT POWER CASE 
define_constants
% NumberOfNodes=6;
[n,e] = single_feeder_gen(250);
% save('6busmodel.mat','n','e');
% load('6busmodel.mat')
% forcing the system to have no transformers, forcing zero capacitance 
% e.funom(1)=e.funom(2);
% e.tunom(1)=e.tunom(2);
% n.unom(1)=n.unom(2);
% e.c=0*e.c;
mpc = matpower_fmt(n,e,60);
mpc = parallel_branch_join(mpc);
mpc.branch(:,BR_B) = 0;
mpc.branch(1,TAP) = 1;
mpc.branch(1,[BR_R,BR_X]) = mpc.branch(1,[BR_R,BR_X])/4;
results = runpf(mpc);
VoltageMatpower=results.bus(:,VM);

%% Variable Setup
NumberOfNodes=size(mpc.bus,1);
FromNode=mpc.branch(:,F_BUS);
ToNode=mpc.branch(:,T_BUS);
Branch=e.id;
NumberOfBranch=length(FromNode);
F=sparse(Branch,FromNode,1,NumberOfBranch,NumberOfNodes);
T=sparse(Branch,ToNode,1,NumberOfBranch,NumberOfNodes);
M0=F-T;
lengths=e.length;

r=mpc.branch(:,BR_R);
x=mpc.branch(:,BR_X);

Rline= sparse(1:NumberOfBranch, 1:NumberOfBranch,r);
Xline= sparse(1:NumberOfBranch, 1:NumberOfBranch,x);
TFT=T*F';
I=speye(NumberOfBranch);
M=M0(:,2:end);

slack_voltage=1.00;

pc=mpc.bus(:,PD)/mpc.baseMVA;
qc=mpc.bus(:,QD)/mpc.baseMVA;

B=(I-TFT)^(-1)-I;
tau=(Rline*B*Rline+Xline*B*Xline)*(Rline*Rline+Xline*Xline)^(-1);
D=(I+B)*T;

%% Lossless
% R=2*(M\Rline)*((I-TFT)\T);
% X=2*(M\Xline)*((I-TFT)\T);
R=2*(M\Rline)*(B+I)*T;
X=2*(M\Xline)*(B+I)*T;

Vlossless=[slack_voltage^2 ;slack_voltage^2+R*pc+X*qc];
Vlossless=sqrt(Vlossless);

Plossless = results.baseMVA*(B+I)*T*pc;
Qlossless = results.baseMVA*(B+I)*T*qc;


%% Alternate Lossy Case

% Rlossy=Rline*D;
% Xlossy=Xline*D;
% m0=M0(:,1);
% K=(I-tau)*M;
% K=(I-2*tau)*M;
% Kinv=K^(-1);
% Vlossy=[slack_voltage^2 ;slack_voltage^2+Kinv*(Rlossy*pc+Xlossy*qc)];
% Vlossy=sqrt(Vlossy);


% taufrac = 0.05;
alpha = 0.4975;
% taufrac = 1/(2*max(tau(:)));
% K=(I-taufrac*tau);
K = (alpha*I - (1-2*alpha)*tau);
Kinv=K^(-1);
% Rlossy=2*(M\Kinv)*Rline*(B+I)*T;
% Xlossy=2*(M\Kinv)*Xline*(B+I)*T;
Rlossy=(M\Kinv)*Rline*(B+I)*T;
Xlossy=(M\Kinv)*Xline*(B+I)*T;
Vlossy=[slack_voltage^2 ;slack_voltage^2+Rlossy*pc+Xlossy*qc];
Vlossy=sqrt(Vlossy);

Plossy = results.baseMVA*(B+I)*(T*pc + ...
    (1-2*alpha)*Rline*(Rline*Rline + Xline*Xline)^(-1)*M0*Vlossy.^2);
Qlossy = results.baseMVA*(B+I)*(T*qc + ...
    (1-2*alpha)*Xline*(Rline*Rline + Xline*Xline)^(-1)*M0*Vlossy.^2); 

%% errors

err = struct('V', VoltageMatpower - [Vlossless, Vlossy],...
             'P', results.branch(:,PF) - [Plossless, Plossy],...
             'Q', results.branch(:,QF) - [Qlossless, Qlossy]);

normerr = [norm(err.V(:,1),2), norm(err.P(:,1),2), norm(err.Q(:,1),2);
           norm(err.V(:,2),2), norm(err.P(:,2),2), norm(err.Q(:,2),2)];
maxerr  = [max(abs(err.V(:,1))), max(abs([err.P(:,1), err.Q(:,1)]),[], 1);
           max(abs(err.V(:,2))), max(abs([err.P(:,2), err.Q(:,2)]),[], 1)];
avgerr  = [mean(abs(err.V(:,1))), mean(abs([err.P(:,1), err.Q(:,1)]), 1);
           mean(abs(err.V(:,2))), mean(abs([err.P(:,2), err.Q(:,2)]), 1)];
stderr  = [std(abs(err.V(:,1))), std(abs([err.P(:,1), err.Q(:,1)]),0, 1);
           std(abs(err.V(:,2))), std(abs([err.P(:,2), err.Q(:,2)]),0, 1)];
%% Plotting
Nodes=n.id;
figure
% plot(Nodes,Vlossless,'r',Nodes,VoltageMatpower,'b',Nodes,Vlossy,'k','linewidth',1.5);
plot(Nodes,VoltageMatpower,'bo-','linewidth',1.5, 'MarkerSize', 10)
hold on;
plot(Nodes,Vlossless,'rx:','linewidth',1.5, 'MarkerSize', 10)
plot(Nodes,Vlossy,'k*--','linewidth',1.5, 'MarkerSize', 10);
legend('VMATPOWER','Vlossless','Vlossy')

figure
subplot(2,1,1)
% plot(Branch,Plossless,'r',...
%      Branch,results.branch(:,PF),'b',...
%      Branch,Plossy,'k','linewidth',1.5);
plot(Branch,results.branch(:,PF),'bo-','linewidth',1.5, 'MarkerSize', 10)
hold on;
plot(Branch,Plossless,'rx:','linewidth',1.5, 'MarkerSize', 10)
plot(Branch,Plossy,'k*--','linewidth',1.5, 'MarkerSize', 10);
legend('PMATPOWER','Plossless','Plossy')
subplot(2,1,2)
% plot(Branch,Qlossless,'r',...
%      Branch,results.branch(:,QF),'b',...
%      Branch,Qlossy,'k','linewidth',1.5);
plot(Branch,results.branch(:,QF),'bo-','linewidth',1.5, 'MarkerSize', 10)
hold on;
plot(Branch,Qlossless,'rx:','linewidth',1.5, 'MarkerSize', 10)
plot(Branch,Qlossy,'k*--','linewidth',1.5, 'MarkerSize', 10);
legend('QMATPOWER','Qlossless','Qlossy')

figure
subplot(2,2,1)
bar([(1:3).',(1:3).'], normerr.')
legend('lossless', 'lossy')
ylabel('||error||_2')
set(gca,'XTickLabels',{'V [p.u]', 'P [MW]', 'Q [MVAr]'})

subplot(2,2,2)
bar([(1:3).',(1:3).'], maxerr.')
legend('lossless', 'lossy')
ylabel('max|error|')
set(gca,'XTickLabels',{'V [p.u]', 'P [MW]', 'Q [MVAr]'})

subplot(2,2,3)
bar([(1:3).',(1:3).'], avgerr.')
legend('lossless', 'lossy')
ylabel('avg|error|')
set(gca,'XTickLabels',{'V [p.u]', 'P [MW]', 'Q [MVAr]'})

subplot(2,2,4)
bar([(1:3).',(1:3).'], stderr.')
legend('lossless', 'lossy')
ylabel('std|error|')
set(gca,'XTickLabels',{'V [p.u]', 'P [MW]', 'Q [MVAr]'})
% figure
% feeder_plot(n,e)


