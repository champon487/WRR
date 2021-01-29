# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 16:13:50 2019

@author: river801
"""

from numba import njit, jit, f8
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import time

#@jit(nopython=True)
#@jit("f8[:](f8[:],f8[:],f8[:],f8[:],f8,f8,f8,f8,f8)", nopython=True)
@jit
def momentum(h,u,wu,q,g,rdx,ib,dt,snm):
    
    h_b = np.zeros(h.shape)
    press = np.zeros(h.shape)
    rough = np.zeros(h.shape)
    advec = np.zeros(h.shape)
    
    h_b[1:nx] = (h[1:nx]+h[2:nx+1])*0.5
    press[1:nx] = -g*(-h[1:nx]+h[2:nx+1])*rdx+g*ib
    rough[1:nx] = g*snm**2.0*u[1:nx]/h_b[1:nx]**(1.333333)
    advec[1:nx] = u[1:nx]*(-u[0:nx-1]+u[1:nx])*rdx
    wu[1:nx] = (u[1:nx]+(-advec[1:nx]+press[1:nx])*dt)/(1.0+rough[1:nx]*dt)
    q[1:nx] = wu[1:nx]*h_b[1:nx]
    
    return wu

#@jit("f8[:](f8[:])", nopython=True, parallel=True)
@jit
def continuity(h,q,rdx,dt):
    
    h[1:nx] = h[1:nx]-(-q[0:nx-1]+q[1:nx])*rdx*dt
    
    return h

@jit
def momentum2(h,u,wu,q,g,rdx,ib,dt,snm,nx):
    
    for i in np.arange(1,nx):
        h_bar = (h[i]+h[i+1])*0.5
        pressure = -g*(-h[i]+h[i+1])*rdx+g*ib
        roughness = g*snm**2.0*u[i]/h_bar**(4.0/3.0)
        advection = u[i]*(-u[i-1]+u[i])*rdx
        wu[i] = (u[i]+(-advection+pressure)*dt)/(1.0+roughness*dt)
        q[i] = wu[i]*h_bar
    
    return wu, q

#@jit("f8[:](f8[:])", nopython=True, parallel=True)
@jit
def continuity2(h,q,rdx,dt,nx):
    
    for i in np.arange(1,nx):
        h[i] = h[i]-(-q[i-1]+q[i])*rdx*dt
    
    return h

# set workiing directly and file name

#path    = "C:/Users/river801/Dropbox/temp/WRR/"
#fname   = "qr-t"

# read data: Time(h), Discharge(m3/s), Precipitation(mm/h)

#with open(path+fname+".txt","r") as fpin:
with open("../statistical_analysis/qr-t.txt","r") as fpin:
    lines = fpin.readlines()
    tt = []
    qq = []
    rr = []
    for line in lines:
        value = line.replace('\n', '').split('\t')
        tt.append(float(value[0]))
        qq.append(float(value[1]))
        rr.append(float(value[2]))
    
fpin.close()

nt = len(tt)-1
tuk   = 3600.
etime = tt[nt]*3600.


# Setting the river geometry and model parameters

g = 9.81

chlen  = 100*100       # length of river reach (m)
wid    = 200            # channel width (m)
snm    = 0.03           # mannings roughness coefs
ib     = 1.0/2000.0      # bed slope

nx  = 500    # number of grid
dx  = chlen/float(nx)    # size of grid
rdx = 1.0/dx
dt  = 0.2    # computational time step (sec)

x   = np.zeros(nx+1)
h   = np.zeros(nx+1)
u   = np.zeros(nx+1)
wu  = np.zeros(nx+1)
q   = np.zeros(nx+1)

h_b   = np.zeros(nx+1)
press = np.zeros(nx+1)
rough = np.zeros(nx+1)
advec = np.zeros(nx+1)


x[0] = 0.

x[1:] = x[0:-1]+dx

#for i in np.arange(1,nx+1):
#    x[i] = x[i-1]+dx
    
q0 = qq[0]
h0 = (snm*q0/(wid*ib**0.5))**0.6

h[:] = h0
u[:] = q0/(h0*wid)
q[:] = q0/wid

#for i in np.arange(nx+1):
#    h[i]  = h0
#    u[i]  = q0/(h0*wid)
#    q[i]  = q0/wid

tt = 0.
optime = 0.

t_cal = []
q_obs = []
q_cal = []

t0 = time.time()
    
while tt<etime:
    nn = int(tt/3600)
    
    q0 = qq[nn]+(qq[nn+1]-qq[nn])/3600.*(tt-nn*3600.)
    h0 = (snm*q0/(wid*ib**0.5))**0.6
    
    h[0] = h0
    u[0] = q0/(h[0]*wid)
    q[0] = q0/wid
    
#    momentum(h,u,wu,q,g,rdx,ib,dt,snm)
    momentum2(h,u,wu,q,g,rdx,ib,dt,snm,nx)
    
#    h_b[1:nx] = (h[1:nx]+h[2:nx+1])*0.5
#    press[1:nx] = -g*(-h[1:nx]+h[2:nx+1])/dx+g*ib
#    rough[1:nx] = g*snm**2.0*u[1:nx]/h_b[1:nx]**(4.0/3.0)
#    advec[1:nx] = u[1:nx]*(-u[0:nx-1]+u[1:nx])/dx
#    wu[1:nx] = (u[1:nx]+(-advec[1:nx]+press[1:nx])*dt)/(1.0+rough[1:nx]*dt)
#    q[1:nx] = wu[1:nx]*h_b[1:nx]
    
       
#    for i in np.arange(1,nx):
#        h_bar = (h[i]+h[i+1])*0.5
#        pressure = -g*(-h[i]+h[i+1])/dx+g*ib
#        roughness = g*snm**2.0*u[i]/h_bar**(4.0/3.0)
#        advection = u[i]*(-u[i-1]+u[i])/dx
#        wu[i] = (u[i]+(-advection+pressure)*dt)/(1.0+roughness*dt)
#        q[i] = wu[i]*h_bar
            
    wu[nx] = wu[nx-1]
    q[nx] = q[nx-1]
    
#    for i in np.arange(1,nx):
#        h[i] = h[i]-(-q[i-1]+q[i])/dx*dt
    
#    continuity(h,q,rdx,dt)
    continuity2(h,q,rdx,dt,nx)
    
#    h[1:nx] = h[1:nx]-(-q[0:nx-1]+q[1:nx])/dx*dt
        
    h[nx] = h[nx-1]
    
    u = wu
    
#    for i in np.arange(1,nx+1):
#        u[i] = wu[i]
    
    if optime>tuk:
        print(tt/3600,q0,q[nx]*wid)
        t_cal.append(tt/3600.)
        q_obs.append(q0)
        q_cal.append(q[nx]*wid)
        optime = optime-tuk
        
#        plt.plot(t_cal,q_obs,color='k', label="Given hydrograph")
#        plt.plot(t_cal,q_cal,color='b', label="Downstream end")
        
    
    optime+=dt
    tt+=dt
    
t1 = time.time()
print("Computational time= ",t1-t0)

# nx = 100
# Original 692 sec   slice 35 sec  jit 14 sec
# jit+型指定  13 sec  +parallelはだめ

# nx = 500
# jit+型指定 34.5 sec
# slice 61.9 sec
# original 3601 sec
    
plt.plot(t_cal,q_obs,color='k', label="Given hydrograph")
plt.plot(t_cal,q_cal,color='b', label="Downstream end")

plt.xlabel( "Time (h)" )
plt.ylabel( "Discharge (m$^3$/s)" )
#plt.set.xlim([0,72])
#plt.set.ylim([0,3000])

plt.legend()

plt.savefig("dynamic.jpg",dpi=300)
plt.close()

