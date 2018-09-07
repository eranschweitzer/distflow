clc
clear 
close all

FromNode=[1 2 3 4 5];
ToNode=[2 3 4 5 6];
Branch=1:5;

%% Volt Var Volt Watt Curves
linear=[0.95,0.98,1.02,1.05,1.02,1.05];
% linear=linear';
no_deadband=[0.95,1.00,1.00,1.05,1.0,1.05];
curved=[0.95,0.98,1.02,1.05,1.05,1.1];

gen_control=[0,0,1,0,0,1];
% gen_control=gen_control*0;
Smaxa=[0,40,20,10,30,50];
% Pmax=MaxCollection(Smaxa);
%%
F=sparse(Branch,FromNode,ones(1,5),5,6);
T=sparse(Branch,ToNode,ones(1,5),5,6);

M=F-T;
r=0.01;
x=0.01;
r=r*[1 1.2 1.2 1.3 0.95];
x=x*[1 1.2 1.2 1.3 0.95];
NumberOfBranch=5;
NumberOfNodes=6;
Rline=r.*speye(NumberOfBranch,NumberOfBranch);
Xline=x.*speye(NumberOfBranch,NumberOfBranch);
TFT=T*F';
I=sparse(eye(size(TFT)));
originalM=M;
M=M(:,2:end);
R=2*(M\Rline)*((I-TFT)\T);
X=2*(M\Xline)*((I-TFT)\T);
slack_voltage=1.01;
% load values
% pc=[0,0.3,0.9,0.2,0,0.3];
% qc=[0,0.8,0.2,0.3,0.1,0.2];
pc=[0,0.3,0.0,0.2,0.9,0.3];
qc=[0,0.0,0.2,0,0.1,0.2];

%% Inverter Specifications

% Number_of_Inverters=4;
% stable case 
% InverterSmax=[20,50,40];
InverterNode=[2,3];
% InverterNode=[3,5,6];
InverterSmax=[40,40];
% droop_curves={'no_deadband','no_deadband'};
droop_curves={'no_deadband','no_deadband','curved'};
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

%%
iteration=300;
VoltageArray=[];
Pg=[];
Qg=[];
for itr=1:iteration
    if (itr==1)
        % Solving the first iteration without inverter
        V=[slack_voltage^2 ;slack_voltage^2+R*pc'+X*qc'];
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
        V=[slack_voltage^2 ; slack_voltage^2+R*(pc'-0.01*pg')+X*(qc'-0.01*qg')];
        VoltageArray=[VoltageArray;V'];
        diffv=VoltageArray(itr,:)-VoltageArray(itr-1,:);
        diffv=abs(diffv);
        if max(diffv) < 0.0001
            break;
        end
    end
end

% disp(sqrt(VoltageArray))





    
   