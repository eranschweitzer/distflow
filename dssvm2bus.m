function bus = dssvm2bus(bus, vm)
%%% add the opendss vm solution to the bus structure

ptr1 = 0;
for k = 1:length(bus)
    ptr2 = ptr1 + length(bus(k).phase);
    bus(k).vm = vm(ptr1+1:ptr2);
    ptr1 = ptr2;
end