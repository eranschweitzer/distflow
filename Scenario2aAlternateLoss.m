clc
clear 
close all

FromNode=[1 2 3 4 5];
ToNode=[2 3 4 5 6];
Branch=1:5;

%% Volt Var Volt Watt Curves
linear=[0.95,0.98,1.02,1.05,1.02,1.05];
%Stable
% no_deadband=[0.95,1.00,1.00,1.05,1.0,1.05];
% Unstable
no_deadband=[0.99,1.00,1.00,1.01,1.0,1.05];
curved=[0.95,0.98,1.02,1.05,1.05,1.1];

%%
F=sparse(Branch,FromNode,ones(1,5),5,6);
T=sparse(Branch,ToNode,ones(1,5),5,6);

M=F-T;
r=0.005;
x=0.005;
r=r*[1 1 1 1 1];
x=x*[1 1 1 1 1];
NumberOfBranch=5;
NumberOfNodes=6;
Rline=r.*speye(NumberOfBranch,NumberOfBranch);
Xline=x.*speye(NumberOfBranch,NumberOfBranch);
TFT=T*F';
I=sparse(eye(size(TFT)));
B=(I-TFT)^(-1)-I;
originalM=M;
M=M(:,2:end);


tau=(Rline*B*Rline+Xline*B*Xline)*(Rline*Rline+Xline*Xline)^(-1);
D=(I+B)*T;
Rlossy=Rline*D;
Xlossy=Xline*D;


K=(I-tau)*M;
Kinv=K^(-1);
slack_voltage=1.03;
% load values
pc=[0,0.3,0.9,0.2,0,0.0];
qc=[0,0.8,0.2,0.3,0.1,0.0];
% pc=[0,0.3,0.8,1.2,0.0,0.0];
% qc=[0,0.8,0.9,0.7,0.0,0.0];

%% Inverter Specifications

% Number_of_Inverters=4;
% stable case 
% InverterSmax=[20,50,40];
InverterNode=[3,4,5,6];
% InverterNode=[3,5,6];
InverterSmax=[0,0,0,0];
% droop_curves={'no_deadband','no_deadband'};
droop_curves={'curved','linear','no_deadband','curved'};
% droop_curves={'curved','curved','curved','curved'};
% droop_curves={'no_deadband','linear','curved','no_deadband'};
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

%% Unstable Case

iteration=300;
tol=1e-4;
VoltageArray=[];
Pg=[];
Qg=[];
for itr=1:iteration
    if (itr==1)
        % Solving the first iteration without inverter
        V=sqrt([slack_voltage^2 ;Kinv*(Rlossy*pc'+Xlossy*qc')]);
        VoltageArray=V';
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
        V=sqrt([slack_voltage^2 ; slack_voltage^2+Kinv*(Rlossy*(pc'-0.01*pg')+Xlossy*(qc'-0.01*qg'))]);
        VoltageArray=[VoltageArray;V'];
%         VoltageArray=sqrt(VoltageArray);
        diffv=VoltageArray(itr,:)-VoltageArray(itr-1,:);
        diffv=abs(diffv);
        if max(diffv) < tol
            break;
        end
    end
end

disp(VoltageArray)
%%
plot(VoltageArray,'linewidth',1.2)
xlim([1 itr])





    
   