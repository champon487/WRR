# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 16:13:50 2019

@author: river801
"""

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

# set workiing directly and file name

path    = "C:/Users/tiwas/Dropbox/temp/WRR/"
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

etime = tt[len(tt)-1]


# Storage routine model

p1 = 0.6
p2 = 0.4648
k1 = 40
k2 = 120

A  = 950    # watershed area (km2)
dt = 0.1    # Calculation time step (h)

t = []
q = []

time = 0

y1 = 0
y2 = 0

while time<etime:
    nn = int(time)
    
    wy2 = y2+(rr[nn]-y1**(1/p2)-k1*p1/p2*y1**(p1/p2-1)*y2)*dt/k2
    wy1 = y1+y2*dt
    wq  = y1**(1./p2)*A/3.6
    
    t.append(time)
    q.append(wq)
    
    y1 = wy1
    y2 = wy2
    
    print("Time=",time,"(h)")
    time+=dt


# Show the observed precipitation and discharge

ax11 = plt.subplot()
ax11.plot( tt, qq, color='k', linewidth=2, label="Observation" )
ax11.plot( t, q, color='r', linewidth=3, label="Calculation" )
ax12 = ax11.twinx()
ax12.bar( tt, rr, color='b' )

ax11.set_xlabel( "Time (h)" )
ax11.set_ylabel( "Discharge (m$^3$/s)" )
ax12.set_ylabel( "Precipitation (mm/h)" )
ax11.set_xlim([0,72])
ax11.set_ylim([0,3000])
ax12.set_ylim([100,0])

ax11.legend()

plt.savefig(path+fname+"_runoff.jpg",dpi=300)
plt.close()
