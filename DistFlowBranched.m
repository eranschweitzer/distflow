clc
clear 
close all

FromNode=[1 2 3 4 5 3 7];
ToNode=[2 3 4 5 6 7 8];
% FromNode=[1 2 3 3 4 7 5];
% ToNode=  [2 3 4 7 5 8 6];
Branch=1:7;
NumberOfNodes=8;
NumberOfBranch=length(Branch);
% shuffle=[1 2 3 4 7 5 8 6];

%%

F=sparse(Branch,FromNode,1,NumberOfBranch,NumberOfNodes);
T=sparse(Branch,ToNode,1,NumberOfBranch,NumberOfNodes);
M=F-T;

slack_voltage=1.03;
r=0.01;
x=0.01;
pc=[0,0.3,0.4,0.2,0.5,0.3,0.3,0.3];
qc=[0,0.8,0.2,0.3,0.1,0.2,0.3,0.3];

% pc=pc(shuffle);
% qc=qc(shuffle);
Rline=r*speye(NumberOfBranch,NumberOfBranch);
Xline=x*speye(NumberOfBranch,NumberOfBranch);

TFT=T*F';
I=speye(size(TFT));
M=M(:,2:end);
R=2*(M\Rline)*((I-TFT)\T);
X=2*(M\Xline)*((I-TFT)\T);

V=[slack_voltage^2 ;slack_voltage^2+R*pc'+X*qc'];
disp(V')
    
   