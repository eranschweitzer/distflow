clc
clear 
close all

FromNode=[1 2 3 4 5];
ToNode=[2 3 4 5 6];
Branch=1:5;


%%
F=sparse(Branch,FromNode,ones(1,5),5,6);
T=sparse(Branch,ToNode,ones(1,5),5,6);

slack_voltage=1.03;
r=0.01;
x=0.01;
pc=[0,0.3,0.4,0.2,0.5,0.3];
qc=[0,0.8,0.2,0.3,0.1,0.2];

M=F-T;
NumberOfBranch=5;
NumberOfNodes=6;
Rline=r*speye(NumberOfBranch,NumberOfBranch);
Xline=x*speye(NumberOfBranch,NumberOfBranch);
TFT=T*F';
I=sparse(eye(size(TFT)));

M=M(:,2:end);
R=2*(M\Rline)*((I-TFT)\T);
X=2*(M\Xline)*((I-TFT)\T);

%%
V=[slack_voltage^2 ;slack_voltage^2+R*pc'+X*qc'];
disp(V')





    
   