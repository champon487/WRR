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
import os
import shutil

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
    
    return h

@jit
def bedload(h,u,qb,qbin,snm,spec,diam,tsc,g,nx):
        
    for i in np.arange(1,nx):
        h_bar = (h[i]+h[i+1])*0.5
        ts = (snm*u[i])**2./(spec*diam*h_bar**(1./3.))
        if ts>tsc:
            qb[i] = 4.0*(ts-tsc)**1.5*(spec*g*diam**3.)**0.5
        else:
            qb[i] = 0.
         
    return qb

@jit
def SscTransport(c,dd,ee,wc,u,h,q,g,snm,wf,rep,dx,dt,nx):
    
    for i in np.arange(1,nx+1):
        uc = (u[i]+u[i-1])*0.5
        usta = (g*(snm*uc)**2./(h[i]**(1./3.)))**0.5
        beta = 15.*wf/usta
        zet  = usta/wf*rep**0.6
        
        dd[i] = wc[i]*alf(beta)*wf
        
        if usta>wf:
            ee[i] = ent_GarciaParker(zet)*wf
        else:
            ee[i] = 0.
        
        c[i] = wc[i]+(-(q[i]*wc[i]-q[i-1]*wc[i-1])/dx+ee[i]-dd[i])*dt/h[i]
        
    return c,dd,ee

@jit
def exner(eta,qb,dd,ee,rdx,dt,poro,nx):
        
    for i in np.arange(1,nx+1):
        eta[i] = eta[i]-(-qb[i-1]+qb[i])*rdx*dt/(1.-poro) + (dd[i]-ee[i])*dt/(1.-poro)
        
    return eta

@jit
def alf(bet):
    if bet>20:
        alfx = bet
    else:
        alfx = bet/(1.-np.exp(-bet))
    
    return alfx

@jit
def ent_GarciaParker(zet):
    
    aa = 1.3*10**(-7)
    
    ent = aa*zet**5./(1.+aa*zet**5./0.3)
    
    return ent

def wf_Dietrich(spec,g,diam,nu):
    
    rep = (spec*g*diam**3.)**0.5/nu
    b1 = 2.891394
    b2 = 0.95296
    b3 = 0.056835
    b4 = 0.002892
    b5 = 0.000245
    
    wf = np.exp(-b1+b2*np.log(rep)-b3*(np.log(rep))**2.-b4*(np.log(rep))**3.+b5*(np.log(rep))**4.)*(spec*g*diam)**0.5

    return wf


def CriticalShieldsNumber_Iwagaki(spec,g,diam,nu):
    
    rep = (spec*g*diam**3.)**0.5/nu
    
    if rep<=2.14:
        ustac2 = 0.14*spec*g*diam
    elif rep>2.14:
        ustac2 = (0.1235*spec*g)**(0.78125)*nu**(0.4375)*diam**(0.34375)
    elif rep>54.2:
        ustac2 = 0.034*spec*g*diam
    elif rep>162.7:
        ustac2 = (0.01505*spec*g)**(1.136364)*nu**(-0.272727273)*diam**(1.40909091)
    elif rep>671:
        ustac2 = 0.05*spec*g*diam
    
    taustac = ustac2/(spec*g*diam)
    
    return taustac

def main(args):

    #set config
#    with open(args[1], 'r') as yml:     # command line?
    with open("config.yml",'r') as yml:
        config = yaml.load(yml, Loader=yaml.FullLoader)

    # Setting the river geometry and model parameters

    g  = config['g']
    nu = config['nu']

    q_min  = config['q_min']          # minimun discharge m3/s
    q_max  = config['q_max']         # maximum discharge m3/s
    t_q    = config['t_q']*3600.   # time from min discharge to maximum discharge (sec)
    chlen  = config['chlen']       # length of river reach (m)
    wid    = config['wid']            # channel width (m)
    snm    = config['snm']           # mannings roughness coefs
    ib     = config['ib']      # bed slope
    spec   = config['spec']
    diam   = config['diam']   # sediment diameter
    poro   = config['poro']
    qbin   = config['qbin']   #sediment supply m2/s

    i_qbin = config['i_qbin']   # 0: equilibrium, 1: specified by qbin

    tuk   = 3600.*config['output_interval']
    etime = 3600.*config['etime']
    morpho_start = 3600.*config['morpho_start']

    nx = config['nx']    # number of grid
    dx = chlen/float(nx)    # size of grid
    dt = config['dt']    # computational time step (sec)
    
    tsc = CriticalShieldsNumber_Iwagaki(spec,g,diam,nu)
    wf = wf_Dietrich(spec,g,diam,nu)
    
    rep = (spec*g*diam**3.)**0.5/nu
    
    path_op = config['path_op']
    
    if os.path.exists(path_op):
        shutil.rmtree(path_op)

    os.mkdir(path_op)

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
    c   = np.zeros(nx+1)
    wc  = np.zeros(nx+1)
    dd  = np.zeros(nx+1)
    ee  = np.zeros(nx+1)

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
            
            # Boundary condition
        
        h[0] = (snm*qt/(wid*ib**0.5))**0.6
        u[0] = qt/(h[0]*wid)
        q[0] = qt/wid
        
        wu[0] = u[0]
        
        h[nx] = h[nx-1]
        
        momentum2(h,eta,u,wu,q,g,rdx,dt,snm,nx)
        
        continuity2(h,q,rdx,dt,nx)
        
        u = wu
        
        wc = c
        
        # boundary condition
        
        usta = (g*(snm*u[0])**2./(h[0]**(1./3.)))**0.5
        
        if usta>wf:
            beta = 15.*wf/usta
            zet  = usta/wf*rep**0.6
            wc[0] = ent_GarciaParker(zet)/alf(beta)
        else:
            wc[0] = 0.
        
        if i_qbin==0:
            qb[0] = qb[1]
        else:       
            qb[0] = qbin
            
        qb[nx] = qb[nx-1]
                
        SscTransport(c,dd,ee,wc,u,h,q,g,snm,wf,rep,dx,dt,nx)
        
        bedload(h,u,qb,qbin,snm,spec,diam,tsc,g,nx)
        
        if tt>morpho_start:
            
            exner(eta,qb,dd,ee,rdx,dt,poro,nx)
            
        if optime>tuk:
            
            print("Time= ",tt/3600)
            optime = optime-tuk
            
            hh = h+eta
            
            plt.title("Time= {0:.2f} hours".format(tt/3600))            
            plt.xlim([0,chlen])
            plt.xlabel( "Downstream distance (m)" ) 
            
            plt.ylim([0,z00+hmax*2.])
                
            plt.plot(xc[1:nx],eta[1:nx],color='k',label='Bed surface')
            plt.plot(xc[1:nx],hh[1:nx],color='b',label='Water surface')
            plt.ylabel( "Elevation (m)" )
            
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
    print("Computational time= ",t1-t0, "  sec")


#--------------------------------------------------
# root
#--------------------------------------------------
if __name__ == '__main__':
    import sys

    # set workiing directly and file name
    args = sys.argv
    main(args)