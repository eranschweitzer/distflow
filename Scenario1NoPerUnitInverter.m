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
    Vsupply=4.16/sqrt(3)*1000; % Getting Line to neutral Voltage
    slack_voltage_pu=1.01;
    slack_voltage=slack_voltage_pu*Vsupply;
    pc=[0,900,270,600,0,1500]/3*1000; % To convert in watts
    qc=[0,2400,600,900,0,1500]/3*1000;
    V2conversionfactor=4.16/sqrt(3)*1000;
end



M=F-T;
NumberOfBranch=5;
NumberOfNodes=6;
Rline=r*speye(NumberOfBranch,NumberOfBranch);
Xline=x*speye(NumberOfBranch,NumberOfBranch);



%% Volt Var Volt Watt Curves
linear=[0.95,0.98,1.02,1.05,1.02,1.05];
%Stable
% no_deadband=[0.95,1.00,1.00,1.05,1.0,1.05];
% Unstable
no_deadband=[0.99,1.00,1.00,1.01,1.0,1.05];
curved=[0.95,0.98,1.02,1.05,1.05,1.1];
BusNames={'sourcebus','bus_02','bus_03','bus_04','bus_05','bus_06'};

%% Inverter Specifications

% Number_of_Inverters=4;
% stable case 
% InverterSmax=[20,50,40];
InverterNode=[3,4,5,6];
% InverterNode=[3,5,6];
InverterSmax=[80,50,50,90];
% droop_curves={'no_deadband','no_deadband'};
% droop_curves={'curved','linear','no_deadband','curved'};
% droop_curves={'curved','curved','curved','curved'};
droop_curves={'no_deadband','linear','curved','no_deadband'};
% droop_curves={'no_deadband','no_deadband','no_deadband','no_deadband'};
droop=containers.Map;
droop('curved')=curved;
droop('linear')=linear;
droop('no_deadband')=no_deadband;
Number_of_Inverters=length(InverterNode);
InverterArray(1:Number_of_Inverters)=InverterDistFlow;
for i = 1:Number_of_Inverters
    InverterArray(i).Smax=InverterSmax(i);
    InverterArray(i)=MaxCollection(InverterArray(i));
    InverterArray(i).Node=InverterNode(i);
    InverterArray(i).setpoints=droop(droop_curves{i});
end

%% lossless case 
TFT=T*F';
I=sparse(eye(size(TFT)));

M=M(:,2:end);
R=2*(M\Rline)*((I-TFT)\T);
X=2*(M\Xline)*((I-TFT)\T);

iteration=300;
tol=1e-4;
VoltageArray=[];
Pg=[];
Qg=[];
% V=[slack_voltage^2 ;slack_voltage^2+R*pc'+X*qc'];
% Vlossless=(sqrt(V/V2conversionfactor^2));

for itr=1:iteration
    if (itr==1)
        % Solving the first iteration without inverter
        V=sqrt([slack_voltage^2 ;slack_voltage^2+R*pc'+X*qc']);
        VoltageArray=V';
        Pg=zeros(1,NumberOfNodes);
        Qg=zeros(1,NumberOfNodes);
    else
        pg=zeros(1,NumberOfNodes);
        qg=zeros(1,NumberOfNodes);
        % Start Iteration for All inverters
        for i = 1:Number_of_Inverters
            pg(1,InverterArray(i).Node)=PCalc(InverterArray(i),V(InverterArray(i).Node)/Vsupply);
            qg(1,InverterArray(i).Node)=QCalc(InverterArray(i),V(InverterArray(i).Node)/Vsupply);
        end
        Pg=[Pg; pg];
        Qg=[Qg; qg];
        V=sqrt([slack_voltage^2 ; slack_voltage^2+R*(pc'-1000*pg')+X*(qc'-1000*qg')]);
        VoltageArray=[VoltageArray;V'];
%         VoltageArray=sqrt(VoltageArray);
        diffv=VoltageArray(itr,:)-VoltageArray(itr-1,:);
        diffv=abs(diffv);
        if max(diffv) < tol
            break;
        end
    end
end
Vlossless=VoltageArray/Vsupply;
Pglossless=Pg;
Qglossless=Qg;
disp(Vlossless)




%% OpenDSS simulation
[DSSObj, DSSText, gridpvpath] = DSSStartup;
DSSCircuit=DSSObj.ActiveCircuit;
DSSSolution=DSSCircuit.Solution;
DSSMon=DSSCircuit.Monitors;
DSSText.command = 'Clear';
DSSText.command = 'Compile C:\feeders\Radial5Busd\simulation21.dss';
setSourceInfo(DSSObj,{'source'},'pu',slack_voltage_pu);
% DSSCircuit.LineCodes.R0=r;
% DSSCircuit.LineCodes.R1=r;
% DSSCircuit.LineCodes.X0=x;
% DSSCircuit.LineCodes.X1=x;
for itr=1:iteration
    if (itr==1)
        % Solving the first iteration without inverter
        DSSSolution.Solve;
        BusInfo=getBusInfo(DSSObj,BusNames);
        V=[BusInfo.voltagePU];
        VoltageArray=V;
        Pg=zeros(1,NumberOfNodes);
        Qg=zeros(1,NumberOfNodes);
    else
        pg=zeros(1,NumberOfNodes);
        qg=zeros(1,NumberOfNodes);
        % Start Iteration for All inverters
        for i = 1:Number_of_Inverters
            pg(1,InverterArray(i).Node)=PCalc(InverterArray(i),V(InverterArray(i).Node));
            qg(1,InverterArray(i).Node)=QCalc(InverterArray(i),V(InverterArray(i).Node));
        end
        Pg=[Pg; pg];
        Qg=[Qg; qg];
        setLoadInfo(DSSObj,{'load_01','load_02','load_03','load_04','load_05','load_06'},'kw',...
            pc*3/1000-pg*3);  % Multiplication by 3 for three phase system, divide by 1000 to convert to kW
        setLoadInfo(DSSObj,{'load_01','load_02','load_03','load_04','load_05','load_06'},'kvar',...
            qc*3/1000-qg*3);
        DSSSolution.Solve;
        BusInfo=getBusInfo(DSSObj,BusNames);
        V=[BusInfo.voltagePU];
        VoltageArray=[VoltageArray;V];
%         VoltageArray=sqrt(VoltageArray);
        diffv=VoltageArray(itr,:)-VoltageArray(itr-1,:);
        diffv=abs(diffv);
        if max(diffv) < tol
            break;
        end
    end
end

VoltageOpenDSS=VoltageArray;
PgOpenDSS=Pg;
QgOpenDSS=Qg;
% Nodes=1:NumberOfNodes;
% plot(Nodes,Vlossless,'r',Nodes,VoltageOpenDSS,'b',Nodes,Vlossy,'k','linewidth',1.5);
% legend('Vlossless','VoltageOpenDSS','Vlossy')

%%
figure
BusNamesLegend={'sourcebus','bus 02','bus 03','bus 04','bus 05','bus 06'};
subplot(211)
plot(Vlossless,'linewidth',1.5)
title('Voltage From LossLessCase')
legend(BusNamesLegend)



subplot(212)
plot(VoltageOpenDSS,'linewidth',1.5)
title('Voltage From OpenDSS')
legend(BusNamesLegend)

%%
figure
BusNamesLegend={'sourcebus','bus 02','bus 03','bus 04','bus 05','bus 06'};
subplot(211)
plot(Pglossless,'linewidth',1.5)
title('Real Power of Inverters LossLessCase')
legend(BusNamesLegend)

subplot(212)
plot(PgOpenDSS,'linewidth',1.5)
title('Real Power of Inverters through OpenDSS')
legend(BusNamesLegend)

%%
figure
BusNamesLegend={'sourcebus','bus 02','bus 03','bus 04','bus 05','bus 06'};
subplot(211)
plot(Qglossless,'linewidth',1.5)
title('Reactive Power of Inverters LossLessCase')
legend(BusNamesLegend)

subplot(212)
plot(QgOpenDSS,'linewidth',1.5)
title('Reactive Power of Inverters through OpenDSS')
legend(BusNamesLegend)


%% making stem plots
figure
stem(Vlossless(end,:),'linewidth',1.5)
hold on
stem(VoltageOpenDSS(end,:),'linewidth',1.5)
legend('lossless','OpenDSS')
title('Node Voltage')

figure
stem(Pglossless(end,:),'linewidth',1.5)
hold on
stem(PgOpenDSS(end,:),'linewidth',1.5)
legend('lossless','OpenDSS')
title('Real Power from Inverters')

figure
stem(Qglossless(end,:),'linewidth',1.5)
hold on
stem(QgOpenDSS(end,:),'linewidth',1.5)
legend('lossless','OpenDSS')
title('Reactive Power from Inverters')

 %% Alternate Lossy
% M=F-T;
% NumberOfBranch=length(Branch);
% NumberOfNodes=max([max(ToNode),max(ToNode)]);
% Rline=r*speye(NumberOfBranch,NumberOfBranch);
% Xline=x*speye(NumberOfBranch,NumberOfBranch);
% TFT=T*F';
% 
% I=speye(size(TFT));
% 
% B=(I-TFT)^(-1)-I;
% originalM=M;
% m0=M(:,1);
% M=M(:,2:end);
% 
% tau=(Rline*B*Rline+Xline*B*Xline)*(Rline*Rline+Xline*Xline)^(-1);
% D=(I+B)*T;
% Rlossy=Rline*D;
% Xlossy=Xline*D;
% 
% 
% K=(I-tau)*M;
% Kinv=K^(-1);
% Vlossy=[slack_voltage^2 ;slack_voltage^2+Kinv*(Rlossy*pc'+Xlossy*qc')];
% Vlossy=(sqrt(Vlossy/V2conversionfactor^2));
