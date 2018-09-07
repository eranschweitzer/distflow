function Q= Qcurve(v,Smax,Pmax,setpoint)
%     defines the Qcurve at various voltage points v given user specified inputs
%     VOLT-VAR POINTS
    vminq = setpoint(1);
    vdead1 = setpoint(2);
    vdead2 = setpoint(3);
    vmaxq = setpoint(4);
%      points below vminq,
    if v <= vminq
        Q = Qmax(Smax,Pmax);
%     linearly decrease Qinj between vminq and vdead1
    elseif ((v > vminq) && (v <= vdead1))
        Q = (vdead1 - v)/(vdead1-vminq)*Qmax(Smax,Pcurve(v,Pmax,setpoint));
%     zero in the dead-band
    elseif (vdead1==vdead2) || ((v > vdead1) && (v <= vdead2))
        Q = 0.0;
%       linearly decrease Qinj between vdead2 and vmaxq+
    elseif ( (v > vdead2) && (v <= vmaxq))
        Q = (vdead2 - v)/(vmaxq - vdead2)*Qmax(Smax,Pcurve(v,Pmax,setpoint));
    
%     maintain the maximum (negative) injection given the watt injection at the given voltage
    elseif  v > vmaxq
        Q = -Qmax(Smax,Pcurve(v,Pmax,setpoint));
    end
end  

