% clc;
% clear ;
% close all;
% 
% 
% %% OpenDSS Integration
% [DSSObj, DSSText, gridpvpath] = DSSStartup;
% DSSCircuit=DSSObj.ActiveCircuit;
% DSSSolution=DSSCircuit.Solution;
% DSSMon=DSSCircuit.Monitors;
% DSSText.command = 'Clear';
% 
% DSSText.command = 'Compile C:\feeders\feeder13_B_R\feeder13BR.dss';
% % Easy way to change parameters for all the loads , making them ZIPV
% % The last parameter refers to the minimum voltage threshhold
% DSSText.Command= 'BatchEdit Load..* Model=1';
% % DSSText.Command='BatchEdit Load..* ZIPV=(0.2,0.05,0.75,0.2,0.05,0.75,0.6)';
% 
% DSSSolution.Solve();
% LoadInfo=getLoadInfo(DSSObj);
% LoadBusNames={LoadInfo.name};
% BusInfo=getBusInfo(DSSObj);
% BusNames={BusInfo.name};
% BusMap=containers.Map;
% BusMap('sourcebus')=1;
% counter=2;
% for bus=2:length(BusNames)
%     test=BusNames(bus);
%     test=test{1};
%     if (sum(ismember(LoadBusNames,test))==0)
%         BusMap(test)=counter;
%         counter=counter+1;
%     end
% end

%% OpenDSS simulation
[DSSObj, DSSText, gridpvpath] = DSSStartup;
DSSCircuit=DSSObj.ActiveCircuit;
DSSSolution=DSSCircuit.Solution;
DSSMon=DSSCircuit.Monitors;
% DSSText.command = 'Clear';
DSSText.command = 'Compile C:\feeders\Radial5Busd\simulation.dss';

BusInfo=getBusInfo(DSSObj);
VoltageOpenDSS=[BusInfo.voltagePU];
