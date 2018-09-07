  # -*- coding: utf-8 -*-
"""
Created on Sun Jul 22 20:08:01 2018
 
@author: shamm
"""


import numpy as np
import matplotlib.pyplot as plt
import math



linear=np.array([0.95,0.98,1.02,1.05,1.02,1.05])
no_deadband=np.array([0.95,1.00,1.00,1.05,1.0,1.05])
curved=np.array([0.95,0.98,1.02,1.05,1.05,1.1])
power_factor=0.9
reactivepowercontribution = math.tan(math.acos(power_factor))
### Change the setpoint variable if other options are needed to be seen according to the three arrays, change the array you want to simulate
setpoint=np.copy(no_deadband)
# setpoint=np.copy(linear)
# setpoint=np.copy(curved)

breakpointlist=[linear,no_deadband,curved,linear,no_deadband,curved]
breakpointlist=[linear,no_deadband,curved,linear,curved,curved]
# breakpointlist=[linear,no_deadband,no_deadband,curved,no_deadband,no_deadband]
breakpointlist=[curved,curved,curved,curved,curved,curved,curved,curved]

# ##### VOLT-VAR POINTS
# vminq = setpoint[0]
# vdead1= setpoint[1]
# vdead2= setpoint[2]
# vmaxq = setpoint[3]
#
# ##### VOLT-WATT POINTS
# vbreakp = setpoint[4]
# vmaxp = setpoint[5]

Qmax = lambda Sout,Pout: np.sqrt(Sout**2 - Pout**2)

def MaxCollection(Smax):
    
    Pmax = Smax ** 2 / (1 + reactivepowercontribution ** 2)
    Pmax = np.sqrt(Pmax)
    return Pmax

def Pcurve(v,Pmax,setpoint):
    ### defines the Pcurve at various voltage points v given user specified inputs

    ##### VOLT-WATT POINTS
    vbreakp = setpoint[4]
    vmaxp = setpoint[5]

    v = np.array(v)
    P = np.zeros(v.shape)
    
    tmp = v <= vbreakp
    P[tmp] = Pmax
    
    tmp = (v > vbreakp) & (v <= vmaxp)
    P[tmp] = (vmaxp - v[tmp])/(vmaxp - vbreakp)*Pmax
    
    tmp = (v > vmaxp)
    P[tmp] = 0.0
    return P

def Qcurve(v,Smax,Pmax,setpoint):
    ### defines the Qcurve at various voltage points v given user specified inputs
    v = np.array(v)
    Q = np.zeros(v.shape)

    ##### VOLT-VAR POINTS
    vminq = setpoint[0]
    vdead1 = setpoint[1]
    vdead2 = setpoint[2]
    vmaxq = setpoint[3]

    ## points below vminq,
    tmp = v <= vminq
    Q[tmp] = Qmax(Smax,Pmax)
    
    ## linearly decrease Qinj between vminq and vdead1
    tmp = (v > vminq) & (v <= vdead1)
    Q[tmp] = (vdead1 - v[tmp])/(vdead1-vminq)*Qmax(Smax,Pcurve(v[tmp],Pmax,setpoint))
    
    ## zero in the dead-band
    if (vdead1==vdead2):
        tmp = (v >= vdead1) & (v <= vdead2)
        Q[tmp] = 0.0
    else:
        tmp = (v > vdead1) & (v <= vdead2)
        Q[tmp] = 0.0
    ## linearly decrease Qinj between vdead2 and vmaxq+
    tmp = (v > vdead2) & (v <= vmaxq)
    Q[tmp] = (vdead2 - v[tmp])/(vmaxq - vdead2)*Qmax(Smax,Pcurve(v[tmp],Pmax,setpoint))
    
   
    
    ## maintain the maximum (negative) injection given the watt injection at the given voltage
    tmp = v > vmaxq
    Q[tmp] = -Qmax(Smax,Pcurve(v[tmp],Pmax,setpoint))
    return Q





# Network Information
lineinfo=[0]
networkset=set(lineinfo)
networklist=[networkset]
for i in range(1,6):
    lineinfo.append(i)
    networkset = set(lineinfo)
    networklist.append(networkset)

lineinfo=[0,1,2,6]
networkset = set(lineinfo)
networklist.append(networkset)
lineinfo=[0,1,2,6,7]
networkset = set(lineinfo)
networklist.append(networkset)


ratio=1
r=0.01
x=r*ratio

r=r*2
x=x*2
# First Node is the 0th node
nodes=np.linspace(0,7,8)
number_of_nodes=len(nodes)
vslack=1.03


rarray={}
rarray['01']=r
rarray['12']=r
rarray['23']=r
rarray['34']=r
rarray['45']=r
rarray['26']=r
rarray['67']=r


xarray={}
xarray['01']=x
xarray['12']=x
xarray['23']=x
xarray['34']=x
xarray['45']=x
xarray['26']=x
xarray['67']=x


v=np.zeros(shape=(number_of_nodes))
pc=1*np.array([0,0.3,0.4,0.2,0.5,0.3,0,0])
# pc=np.array([0,0,0.0,0,0,4])
qc=1*np.array([0,0.8,0.2,0.3,0.1,0.2,0,0])
gen_control=np.array([0,0,1,0,0,1,0,0])
gen_control=0*gen_control
Smaxa=np.array([0,40,20,10,30,50,20,20])

tol=1e-4
iteration=300
basecase=0
V=np.zeros(shape=(iteration,number_of_nodes))
Power=np.zeros(shape=(iteration,2))
inverterloc=np.where(gen_control>0.0)[0]


V=np.zeros(shape=(iteration,number_of_nodes))
Power=np.zeros(shape=(iteration,2))
RealPower=np.zeros(shape=(iteration,number_of_nodes))
ReactivePower=np.zeros(shape=(iteration,number_of_nodes))
print('\n')
for itr in range(0,iteration):
    #doing a flat run analysis without any PV penetration, only the first time
    if (itr>0):
        basecase=1
    for i in range(1, number_of_nodes):
        # print('\n')
        # print('i:' + str(i)+'\n')
        sumpq = 0
        
        for j in range(1, number_of_nodes):
            # print('j:' + str(j)+'\n')
            Pmax = MaxCollection(Smaxa[j])
            Smax=Smaxa[j]

            sumr=0
            sumx=0
            p = list(networklist[i] & networklist[j])
            for k in range(len(p)-1):
                sumr += rarray[str(k)+str(k+1)]
                sumx += xarray[str(k) + str(k + 1)]
            sumpq=sumpq-pc[j]*sumr-qc[j]*sumx +basecase*gen_control[j]*Pcurve(v[j],Pmax,breakpointlist[j])/100*sumr+basecase*gen_control[j]*Qcurve(v[j],Smax,Pmax,breakpointlist[j])/100*sumx
            RealPower[itr,j]=gen_control[j]*Pcurve(v[j],Pmax,breakpointlist[j])
            ReactivePower[itr,j]=gen_control[j]*Qcurve(v[j],Smax,Pmax,breakpointlist[j])

        v[i] = vslack ** 2 + sumpq
    v[0] = vslack ** 2
    V[itr, :] = v
    print('Iteration:', str(itr), '  ', str(v))
    if (itr>0):
        diffv=abs(V[itr,:]-V[itr-1,:])
        #diffv = abs(V[itr, inverterloc] - V[itr - 1, inverterloc])
        if (max(diffv)<tol):
            break
#print(RealPower[0:itr,inverterloc])
print('\n')
#print(ReactivePower[0:itr,inverterloc])