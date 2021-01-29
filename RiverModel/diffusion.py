# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 16:13:50 2019

@author: river801
"""

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

# set workiing directly and file name

path    = "C:/Users/river801/Dropbox/temp/WRR/"
fname   = "qr-t"

# read data: Time(h), Discharge(m3/s), Precipitation(mm/h)

with open(path+fname+".txt","r") as fpin:
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
tuk   = 200.
etime = tt[nt]*3600.


# Setting the river geometry and model parameters

chlen  = 50*1000       # length of river reach (m)
wid    = 200            # channel width (m)
snm    = 0.03           # mannings roughness coefs
ib     = 1.0/1000.0      # bed slope

nx = 50    # number of grid
dx = chlen/float(nx)    # size of grid
dt = 1.0    # computational time step (sec)

dx2 = dx*dx
dx3 = dx*dx2

x   = np.zeros(nx+1)
h   = np.zeros(nx+1)
wh  = np.zeros(nx+1)
gh  = np.zeros(nx+1)
wgh = np.zeros(nx+1)
u   = np.zeros(nx+1)
wu  = np.zeros(nx+1)
q   = np.zeros(nx+1)
wq  = np.zeros(nx+1)

x[0] = 0.

for i in np.arange(1,nx+1):
    x[i] = x[i-1]+dx
    
q0 = qq[0]
h0 = (snm*q0/(wid*ib**0.5))**0.6

for i in np.arange(nx+1):
    h[i]  = h0
    u[i]  = q0/(h0*wid)
    gh[i] = 0.
    q[i]  = q0

time = 0.
optime = 0.

t_cal = []
q_obs = []
q_cal = []

fpout = open(path+"out.txt",'w')
    
while time<etime:
    nn = int(time/3600)
    
    q0 = qq[nn]+(qq[nn+1]-qq[nn])/3600.*(time-nn*3600.)
    h0 = (snm*q0/(wid*ib**0.5))**0.6
    
    h[0] = h0
    gh[0] = 0.
    u[0] = q0/(h[0]*wid)
    qq[0] = q0
    
    for i in np.arange(1,nx):
        dhdx = ib-(-h[i-1]+h[i+1])/dx*0.5
        u[i] = h[i]**(2./3.)*dhdx**0.5/snm
        
    u[nx] = u[nx-1]
    h[nx] = h[nx-1]
    gh[nx] = gh[nx-1]
    
    for i in np.arange(1,nx):
        dhdx = ib-(-h[i-1]+h[i+1])/dx*0.5
        diff = 0.5*h[i]**(5.0/3.0)/dhdx**0.5*(h[i-1]-2.0*h[i]+h[i+1])/dx2/snm
        wh[i]	= h[i]+diff*dt
        
    wh[0] = h0
    wh[nx] = wh[nx-1]
    
    for i in np.arange(1,nx):
        wgh[i] = gh[i]+(-wh[i-1]+wh[i+1]+h[i-1]-h[i+1])*0.5/dx
        
    wgh[0] = 0.
    wgh[nx] = wgh[nx-1]
    
    for i in np.arange(1,nx):
        dhdx = ib-(-wh[i-1]+wh[i+1])/dx*0.5
        u[i] = wh[i]**(2./3.)*dhdx**0.5/snm
        
    
    for i in np.arange(1,nx):
        xx = -u[i]*dt*5./3.
        
        if u[i]>0:
            isn = 1.
        else:
            isn = -1.
            
        im1 = i-int(isn)
        
        a1 = ((wgh[im1]+wgh[i])*dx*float(isn)-2.0*(wh[i]-wh[im1]))/(dx3*float(isn))
        b1 = (3.0*(wh[im1]-wh[i])+(wgh[im1]+2.0*wgh[i])*dx*float(isn))/dx2
        h[i]	= ((a1*xx+b1)*xx+wgh[i])*xx+wh[i]
        gh[i]	= (3.0*a1*xx+2.0*b1)*xx+wgh[i]
        
    for i in np.arange(1,nx):
        gh[i] = gh[i]-(gh[i]*(-u[i-1]+u[i+1]))*0.5*dt/dx
    
    for i in np.arange(nx+1):
        q[i] = h[i]*u[i]*wid
    
    if optime>tuk:
        print(time/3600,q0,q[nx])
        
        output = str(time/3600)+"\t"+str(q0)+"\t"+str(q[nx])+"\n"
        
        fpout.writelines(output)
        
        t_cal.append(time/3600.)
        q_obs.append(q0)
        q_cal.append(q[nx])
        optime = optime-tuk
    
    optime+=dt
    time+=dt
    
fpout.close()
    
plt.plot(t_cal,q_obs,color='k', label="Given hydrograph")
plt.plot(t_cal,q_cal,color='b', label="Downstream end")

plt.xlabel( "Time (h)" )
plt.ylabel( "Discharge (m$^3$/s)" )
#plt.set.xlim([0,72])
#plt.set.ylim([0,3000])

plt.legend()

plt.savefig(path+fname+"-diffusion.jpg",dpi=300)
plt.close()

