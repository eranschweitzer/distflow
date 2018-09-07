function P= Pcurve(v,Pmax,setpoint)
%     defines the Pcurve at various voltage points v given user specified inputs

%    VOLT-WATT POINTS
    vbreakp = setpoint(5);
    vmaxp = setpoint(6);

    if v <= vbreakp
        P = Pmax;
      
    elseif ((v > vbreakp) && (v <= vmaxp))
        P = (vmaxp - v)/(vmaxp - vbreakp)*Pmax;
    
    else
        P = 0.0;
    end
end