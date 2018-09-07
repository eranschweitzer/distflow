%%
clc
clear 
close all

FromNode=[1 2 3 4 5];
ToNode=[2 3 4 5 6];
Branch=1:5;

%%
F=sparse(Branch,FromNode,ones(1,5),5,6);
T=sparse(Branch,ToNode,ones(1,5),5,6);

Vbase=4.16*1000;
Sbase=1000;
Zbase=Vbase*Vbase/Sbase;
UsePerUnit=0;

if (UsePerUnit==1)
    slack_voltage=4.16*1000/Vbase;
    r=1.993*10/Zbase;
    x=1.456*10/Zbase;
    pc=[0,0,0,0,0,1800]/Sbase;
    qc=[0,0,0,0,0,1800*.484]/Sbase;
    V2conversionfactor=1;
else
    r=0.0173;
    x=0.0173;
    slack_voltage=1.03*4.16/sqrt(3)*1000;
    pc=[0,900,270,600,0,1500]/3*1000;
    qc=[0,2400,600,900,0,1500]/3*1000;
    V2conversionfactor=4.16/sqrt(3)*1000;
end



M=F-T;
NumberOfBranch=5;
NumberOfNodes=6;
Rline=r*speye(NumberOfBranch,NumberOfBranch);
Xline=x*speye(NumberOfBranch,NumberOfBranch);

%% lossless case 
TFT=T*F';
I=sparse(eye(size(TFT)));

M=M(:,2:end);
R=2*(M\Rline)*((I-TFT)\T);
X=2*(M\Xline)*((I-TFT)\T);

iteration=300;
VoltageArray=[];
Pg=[];
Qg=[];
V=[slack_voltage^2 ;slack_voltage^2+R*pc'+X*qc'];
Vlossless=(sqrt(V/V2conversionfactor^2));

Plossless= (I-TFT)^(-1)*T*pc';
Qlossless= (I-TFT)^(-1)*T*qc';


%% Alternate Lossy
M=F-T;
NumberOfBranch=length(Branch);
NumberOfNodes=max([max(ToNode),max(ToNode)]);
Rline=r*speye(NumberOfBranch,NumberOfBranch);
Xline=x*speye(NumberOfBranch,NumberOfBranch);
TFT=T*F';

I=speye(size(TFT));

B=(I-TFT)^(-1)-I;
originalM=M;
m0=M(:,1);
M=M(:,2:end);

tau=(Rline*B*Rline+Xline*B*Xline)*(Rline*Rline+Xline*Xline)^(-1);
D=(I+B)*T;
Rlossy=Rline*D;
Xlossy=Xline*D;


K=(I-tau)*M;
Kinv=K^(-1);
Vlossy=[slack_voltage^2 ;slack_voltage^2+Kinv*(Rlossy*pc'+Xlossy*qc')];
Vlossy=(sqrt(Vlossy/V2conversionfactor^2));


%% OpenDSS simulation
[DSSObj, DSSText, gridpvpath] = DSSStartup;
DSSCircuit=DSSObj.ActiveCircuit;
DSSSolution=DSSCircuit.Solution;
DSSMon=DSSCircuit.Monitors;
DSSText.command = 'Clear';
DSSText.command = 'Compile C:\feeders\Radial5Busd\simulation21.dss';
% DSSCircuit.LineCodes.R0=r;
% DSSCircuit.LineCodes.R1=r;
% DSSCircuit.LineCodes.X0=x;
% DSSCircuit.LineCodes.X1=x;
DSSSolution.Solve;
BusInfo=getBusInfo(DSSObj,{'sourcebus','bus_02','bus_03','bus_04','bus_05','bus_06'});
VoltageOpenDSS=[BusInfo.voltagePU];
Nodes=1:NumberOfNodes;
plot(Nodes,Vlossless,'r',Nodes,VoltageOpenDSS,'b',Nodes,Vlossy,'k','linewidth',1.5);
legend('Vlossless','VoltageOpenDSS','Vlossy')
