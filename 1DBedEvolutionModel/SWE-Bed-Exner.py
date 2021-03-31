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
import yaml

@jit
def momentum2(h,eta,u,wu,q,g,rdx,dt,snm,nx):
    
    for i in np.arange(1,nx):
        h_bar = (h[i]+h[i+1])*0.5
        pressure = -g*(-h[i]+h[i+1])*rdx-g*(-eta[i]+eta[i+1])*rdx
        roughness = g*snm**2.0*u[i]/h_bar**(4.0/3.0)
        advection = u[i]*(-u[i-1]+u[i])*rdx
        wu[i] = (u[i]+(-advection+pressure)*dt)/(1.0+roughness*dt)
        q[i] = wu[i]*h_bar
        
        wu[nx] = wu[nx-1]
        q[nx] = q[nx-1]
    
    return wu, q

@jit
def continuity2(h,q,rdx,dt,nx):
    
    for i in np.arange(1,nx):
        h[i] = h[i]-(-q[i-1]+q[i])*rdx*dt
        
    h[nx] = h[nx-1]
    
    return h

@jit
def bedload(h,u,qb,qbin,snm,spec,diam,tsc,g,nx):
        
    for i in np.arange(1,nx):
        h_bar = (h[i]+h[i+1])*0.5
        ts = (snm*u[i])**2./(spec*diam*h_bar**(1./3.))
        if ts>tsc:
            qb[i] = 8.0*(ts-tsc)**1.5*(spec*g*diam**3.)**0.5
        else:
            qb[i] = 0.
            
    qb[ 0] = qbin
    qb[nx] = qb[nx-1]
         
    return qb

@jit
def exner(eta,qb,rdx,dt,poro,nx):
        
    for i in np.arange(1,nx+1):
        eta[i] = eta[i]-(-qb[i-1]+qb[i])*rdx*dt/(1.-poro)
        
    return eta

def main(args):

    #set config
    with open(args[1], 'r') as yml:
        config = yaml.load(yml, Loader=yaml.FullLoader)

    # Setting the river geometry and model parameters

    g  = config['g']
    # nu = config['nu']

    q_min  = config['q_min']          # minimun discharge m3/s
    q_max  = config['q_max']         # maximum discharge m3/s
    t_q    = config['t_q']*3600.   # time from min discharge to maximum discharge (sec)
    chlen  = config['chlen']       # length of river reach (m)
    wid    = config['wid']            # channel width (m)
    snm    = config['snm']           # mannings roughness coefs
    ib     = config['ib']      # bed slope
    spec   = config['spec']
    diam   = config['diam']   # sediment diameter
    tsc    = config['tsc']
    poro   = config['poro']
    qbin   = config['qbin']   #sediment supply m2/s

    i_vis = config['i_vis']

    tuk   = 3600.*config['output_interval']
    etime = t_q*2.*config['hyd_cycle']

    nx = config['nx']    # number of grid
    dx = chlen/float(nx)    # size of grid
    dt = config['dt']    # computational time step (sec)

    rdx = 1./dx

    x   = np.zeros(nx+1)
    z   = np.zeros(nx+1)
    h   = np.zeros(nx+1)
    u   = np.zeros(nx+1)
    wu  = np.zeros(nx+1)
    q   = np.zeros(nx+1)
    hh  = np.zeros(nx+1)
    xc  = np.zeros(nx+1)
    eta = np.zeros(nx+1)
    qb  = np.zeros(nx+1)

    dzdx = np.zeros(nx+1)

    x[0] = 0.

    for i in np.arange(1,nx+1):
        x[i] = x[i-1]+dx
        

    z00  = chlen*ib
    z[0] = z00

    for i in np.arange(1,nx+1):
        z[i] = z[i-1]-ib*dx
            
    for i in np.arange(1,nx+1):
        eta[i] = (z[i-1]+z[i])*0.5
        xc[i] = (x[i-1]+x[i])*0.5
        
    h0 = (snm*q_min/(wid*ib**0.5))**0.6

    hmax = (snm*q_max/(wid*ib**0.5))**0.6

    h[:] = h0
    u[:] = q_min/(h0*wid)
    q[:] = q_min/wid

    tt = 0.
    t_hyd = 0.
    optime = 0.

    n = 0

    t0 = time.time()
        
    while tt<etime:
        
        if t_hyd<t_q:
            qt = q_min+(q_max-q_min)/t_q*t_hyd
        else:
            qt = q_max-(q_max-q_min)/t_q*(t_hyd-t_q)
        
        h[0] = (snm*qt/(wid*ib**0.5))**0.6
        u[0] = qt/(h[0]*wid)
        q[0] = qt/wid
        
        momentum2(h,eta,u,wu,q,g,rdx,dt,snm,nx)
        
        continuity2(h,q,rdx,dt,nx)
        
        u = wu
        
        bedload(h,u,qb,qbin,snm,spec,diam,tsc,g,nx)
            
        exner(eta,qb,rdx,dt,poro,nx)
            
        if optime>tuk:
            
            print("Time= ",tt/3600)
            optime = optime-tuk
            
            hh = h+eta
            
            plt.xlim([0,chlen])
            plt.xlabel( "Downstream distance (m)" ) 
            
            if i_vis==0:
            
                plt.ylim([0,z00+hmax*10.])
                
                plt.plot(xc[1:nx],eta[1:nx],color='k',label='Bed surface')
                plt.plot(xc[1:nx],hh[1:nx],color='b',label='Water surface')
                plt.ylabel( "Elevation (m)" )
            
            else:
            
                plt.ylim([0.5,1.5])
                
                dzdx[1:nx] = -(-eta[1:nx]+eta[2:nx+1])*rdx/ib
                
                plt.plot(xc[1:nx],dzdx[1:nx],color='b',label='slope')
                
                plt.ylabel( "slope/initial slope" )
            
            plt.legend()
            
            nnn = str(n)
            
            plt.savefig('fig/Figure' + nnn.zfill(4) +".jpg", dpi=300)
            
            plt.close()
            
            n += 1
        
        optime+=dt
        tt+=dt
        t_hyd+=dt
        
        if t_hyd>2.*t_q:
            t_hyd = t_hyd-2.*t_q
        
    t1 = time.time()
    print("Computational time= ",t1-t0)

    # nx=100 t=100
    # jit 19 sec
    # normal 270 sec
    # normal 319 sec

    # nx=100 t=300
    # jit 53 sec
    # normal 806 sec

#--------------------------------------------------
# root
#--------------------------------------------------
if __name__ == '__main__':
    import sys

    # set workiing directly and file name
    args = sys.argv
    main(args)