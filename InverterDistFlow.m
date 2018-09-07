classdef InverterDistFlow
   properties  (Access=public)
     Smax=0
     setpoints=zeros(1,5);
     power_factor=0.9;
     Pmax=0;
     Node=0;
   end
   
   methods
       
   
   function obj= MaxCollection(obj)
        reactivepowercontribution=tan(acos(obj.power_factor));
        pmax = sqrt(obj.Smax * obj.Smax / (1 + reactivepowercontribution ^2 ));
        obj.Pmax=pmax;
    end   
       
    function [Q]= QCalc(obj,v)
%     defines the Qcurve at various voltage points v given user specified inputs
%     VOLT-VAR POINTS
        setpoint=obj.setpoints;
    
        vminq = setpoint(1);
        vdead1 = setpoint(2);
        vdead2 = setpoint(3);
        vmaxq = setpoint(4);

%      points below vminq,
        if v <= vminq
            Q = Qmax(obj.Smax,obj.Pmax);
%     linearly decrease Qinj between vminq and vdead1
        elseif ((v > vminq) && (v <= vdead1))
            Q = (vdead1 - v)/(vdead1-vminq)*Qmax(obj.Smax,Pcurve(v,obj.Pmax,setpoint));
    %     zero in the dead-band
        elseif (vdead1==vdead2) || ((v > vdead1) && (v <= vdead2))
            Q = 0.0;
%       linearly decrease Qinj between vdead2 and vmaxq+
        elseif ( (v > vdead2) && (v <= vmaxq))
            Q = (vdead2 - v)/(vmaxq - vdead2)*Qmax(obj.Smax,Pcurve(v,obj.Pmax,setpoint));
%             Q = (vdead2 - v)/(vmaxq - vdead2)*Qmax(obj.Smax,PCalc(obj,v));
    
%     maintain the maximum (negative) injection given the watt injection at the given voltage
        elseif  v > vmaxq
            Q = -Qmax(obj.Smax,Pcurve(v,obj.Pmax,setpoint));
        end
    end
    
    function P= PCalc(obj,v)
%     defines the Pcurve at various voltage points v given user specified inputs

%    VOLT-WATT POINTS
        setpoint=obj.setpoints;
        vbreakp = setpoint(5);
        vmaxp = setpoint(6);

       if v <= vbreakp
            P = obj.Pmax;
       elseif ((v > vbreakp) && (v <= vmaxp))
            P = (vmaxp - v)/(vmaxp - vbreakp)*obj.Pmax;
       else
            P = 0.0;
       end
    end
    
  end
end

