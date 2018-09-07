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

gen_control=[0,0,0,0,0,0];
Smaxa=[0,40,20,10,30,50];
Pmax=MaxCollection(Smaxa);
%%
F=sparse(Branch,FromNode,ones(1,5),5,6);
T=sparse(Branch,ToNode,ones(1,5),5,6);

Vbase=4.16*1000;
Sbase=1000;
Zbase=Vbase*Vbase/Sbase;

%     slack_voltage=1.03;
%     r=0.01;
%     x=0.02;
%     pc=[0,0.3,0.4,0.2,0.5,0.3];
%     qc=[0,0.8,0.2,0.3,0.1,0.2];

UsePerUnit=0;

if (UsePerUnit==1)
    slack_voltage=4.16*1000/Vbase;
    r=1.993*10/Zbase;
    x=1.456*10/Zbase;
    pc=[0,0,0,0,0,1800]/Sbase;
    qc=[0,0,0,0,0,1800*.484]/Sbase;
    V2conversionfactor=1;
else
    r=1.993*10;
    x=1.456*10;
    slack_voltage=4.16*1000;
    pc=[0,0,0,0,0,1800];
    qc=[0,0,0,0,0,1800*.484];
    V2conversionfactor=Vbase*Vbase;
end

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
iteration=300;
VoltageArray=[];
Pg=[];
Qg=[];
for itr=1:iteration
    if (itr==1)
        V=[slack_voltage^2 ;slack_voltage^2+R*pc'+X*qc'];
        VoltageArray=V';
        Pg=zeros(1,NumberOfNodes);
        Qg=zeros(1,NumberOfNodes);
    else
        pg=zeros(1,NumberOfNodes);
        qg=zeros(1,NumberOfNodes);
        for i=1:NumberOfNodes
            if (gen_control(i)>0)          
                pg(1,i)=Pcurve(V(i),Pmax(i),curved);
                qg(1,i)=Qcurve(V(i),Smaxa(i),Pmax(i),curved);
            end
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

disp(sqrt(VoltageArray/V2conversionfactor))
DistFlowAlternateLoss


function Pmax= MaxCollection(Smax)
    reactivepowercontribution=0.484;
    Pmax = Smax .^ 2 / (1 + reactivepowercontribution ^2 );
    Pmax = sqrt(Pmax);
end



    
   