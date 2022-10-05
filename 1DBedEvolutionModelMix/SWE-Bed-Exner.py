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
def momentum2(h,eta,u,wu,q,g,rdx,dt,snmm,nx):
    
    for i in np.arange(1,nx):
        h_bar = (h[i]+h[i+1])*0.5
        snm_bar = (snmm[i]+snmm[i+1])*0.5
        pressure = -g*(-h[i]+h[i+1])*rdx-g*(-eta[i]+eta[i+1])*rdx
        roughness = g*snm_bar**2.0*u[i]/h_bar**(4.0/3.0)
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
def bedload(h,u,qbk,qbin,snmm,spec,dk,tsck,g,nx,nk):
    
    for k in np.arange(0,nk+1):
        for i in np.arange(1,nx):
            h_bar = (h[i]+h[i+1])*0.5
            snm_bar = (snmm[i]+snmm[i+1])*0.5
            ts = (snm_bar*u[i])**2./(spec*dk[k]*h_bar**(1./3.))
            if ts>tsck[k][i]:
                qbk[k][i] = 4.0*(ts-tsck[k][i])**1.5*(spec*g*dk[k]**3.)**0.5
            else:
                qbk[k][i] = 0.
         
    return qbk

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
def exner(eta,qbk,dz,dzk,rdx,dt,poro,nx,nk):
    
    dz[:] = 0.

    for k in np.arange(0,nk+1):    
        for i in np.arange(1,nx+1):
            dzk[k][i] = -(-qbk[k][i-1]+qbk[k][i])*rdx*dt/(1.-poro)
            dz[i] = dz[i]+dzk[k][i]
    
    for i in np.arange(1,nx+1):
        eta[i] = eta[i]+dz[i]
        
    return eta, dz, dzk

@jit
def sorting(fak,ftk,fdk,dz,dzk,e_t,e_d,nb,fak_new,ftk_new,fdk_new,nx,nk,poro):
    
    for i in np.arange(1,nx+1):
        if dz[i]<0.:
            if (e_t[i]+dz[i])>0.:
                e_t_new = e_t[i]+dz[i]
                nb_new = nb[i]
                for k in np.arange(0,nk+1):
                    fak_new[k] = fak[k][i]-ftk[k][i]*dz[i]/(1.-poro)+dzk[k][i]/(1.-poro)
                    ftk_new[k] = ftk[k][i]
                    fdk_new[k] = fdk[int(nb[i])][k][i]
            else:
                e_t_new = e_t[i]+dz[i]+e_d
                nb_new = nb[i]-1
                for k in np.arange(0,nk+1):
                    fak_new[k] = fak[k][i]+ftk[k][i]*e_t[i]/(1.-poro)-(e_t[i]+dz[i])*fdk[int(nb[i])][k][i]/(1.-poro)+dzk[k][i]/(1.-poro)
                    ftk_new[k] = fdk[int(nb[i])][k][i]
                    fdk_new[k] = fdk[int(nb[i])-1][k][i]
        else:
            if (e_t[i]+dz[i])<e_d:
                e_t_new = e_t[i]+dz[i]
                nb_new = nb[i]
                for k in np.arange(0,nk+1):
                    fak_new[k] = fak[k][i]-fak[k][i]*dz[i]/(1.-poro)+dzk[k][i]/(1.-poro)
                    ftk_new[k] = (e_t[i]*ftk[k][i]+dz[i]*fak[k][i])/e_t_new
                    fdk_new[k] = fdk[int(nb[i])][k][i]
            else:
                e_t_new = e_t[i]+dz[i]-e_d
                nb_new = nb[i]+1
                for k in np.arange(0,nk+1):
                    fak_new[k] = fak[k][i]*(1.-dz[i]/(1.-poro))+dzk[k][i]/(1.-poro)
                    ftk_new[k] = fak[k][i]
                    fdk_new[k] = ftk[k][i]*e_t[i]/e_d+(1.-e_t[i]/e_d)*fak[k][i]
                    
        e_t[i] = e_t_new
        nb[i]  = nb_new
        
        fak_new[fak_new<0] = 0.
        f_tot = np.sum(fak_new)
        for k in np.arange(0,nk+1):
            fak_new[k] = fak_new[k]/f_tot  
        
        ftk_new[ftk_new<0] = 0.
        f_tot = np.sum(ftk_new)
        for k in np.arange(0,nk+1):
            ftk_new[k] = ftk_new[k]/f_tot
        
        fdk_new[fdk_new<0] = 0.
        f_tot = np.sum(fdk_new)
        for k in np.arange(0,nk+1):
            fdk_new[k] = fdk_new[k]/f_tot
        
        for k in np.arange(0,nk+1):
            fak[k][i] = fak_new[k]
            ftk[k][i] = ftk_new[k]
            fdk[int(nb[i])][k][i] = fdk_new[k]
            
    return fak, ftk, fdk
            

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
    if rep>2.14:
        ustac2 = (0.1235*spec*g)**(0.78125)*nu**(0.4375)*diam**(0.34375)
    if rep>54.2:
        ustac2 = 0.034*spec*g*diam
    if rep>162.7:
        ustac2 = (0.01505*spec*g)**(1.136364)*nu**(-0.272727273)*diam**(1.40909091)
    if rep>671:
        ustac2 = 0.05*spec*g*diam
    
    taustac = ustac2/(spec*g*diam)
    
    return taustac

@jit
def RepresentativeDiameter(dm,dk,fak,nx,nk):
    
    dm[:] = 0.
    
    for k in np.arange(0,nk+1):
        for i in np.arange(1,nx+1):
            dm[i] = dm[i]+dk[k]*fak[k][i]
            
    return dm

@jit
def ManningN(snmm,dm,g,nx):
    
    for i in np.arange(1,nx+1):
        snmm[i] = (2.5*dm[i])**(1./6.)/(7.66*g**0.5)
            
    return snmm

@jit
def hidingexposure(tsck,tscm,dm,dk,spec,g,nu,nx,nk):
    
    for i in np.arange(1,nx+1):
        tscm[i] = 0.05  #CriticalShieldsNumber_Iwagaki(spec,g,dm[i],nu)
        
    for k in np.arange(0,nk+1):
        for i in np.arange(1,nx+1):
            xi_eh = (np.log10(23.)/np.log10(21.*dk[k]/dm[i]+2.))**2.
            tsck[k][i] = xi_eh*tscm[i]
            
    return tscm,tsck

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
    
    nk = 2
    
    dk   = np.zeros(nk+1)
    fk   = np.zeros(nk+1)
    wfk  = np.zeros(nk+1)
    rep  = np.zeros(nk+1)
    tsck = np.zeros(nk+1)
    fak_new = np.zeros(nk+1)
    ftk_new = np.zeros(nk+1)
    fdk_new = np.zeros(nk+1)
    
    dk[0] = 0.0001
    dk[1] = 0.001
    dk[2] = 0.005
    
    fk[0] = 1./3.
    fk[1] = 1./3.
    fk[2] = 1./3.
    
    for k in np.arange(0,nk+1):
        wfk[k]  = wf_Dietrich(spec,g,dk[k],nu)
        rep[k]  = (spec*g*dk[k]**3.)**0.5/nu
        
    e_m = dk[nk]
    e_d = 0.1
    nm  = 20
    
    
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
    dz  = np.zeros(nx+1)
    eta = np.zeros(nx+1)
    
    nb = np.zeros(nx+1)
    e_t = np.zeros(nx+1)
    eta_base = np.zeros(nx+1)
    dm = np.zeros(nx+1)
    snmm = np.zeros(nx+1)
    tscm = np.zeros(nx+1)
    
    dzk  = np.zeros((nk+1,nx+1))
    fak  = np.zeros((nk+1,nx+1))
    ftk  = np.zeros((nk+1,nx+1))
    qbk  = np.zeros((nk+1,nx+1))
    tsck = np.zeros((nk+1,nx+1))
    
    fdk  = np.zeros((nm,nk+1,nx+1))
    
    # c   = np.zeros(nx+1)
    # wc  = np.zeros(nx+1)
    # dd  = np.zeros(nx+1)
    # ee  = np.zeros(nx+1)

    x[0] = 0.

    for i in np.arange(1,nx+1):
        x[i] = x[i-1]+dx
        
   
    z00  = chlen*ib
    z[0] = z00
    
    for i in np.arange(1,nx+1):
        if x[i]<x[int(nx*0.5)]:
            z[i] = z[i-1]-ib*dx
        else:
            z[i] = z[i-1]-ib*0.2*dx
            
    for i in np.arange(1,nx+1):
        eta[i] = (z[i-1]+z[i])*0.5
        xc[i] = (x[i-1]+x[i])*0.5
        eta_base[i] = eta[i] - e_d*10.
        
    for i in np.arange(1,nx+1):
        nb[i] = 1
        e_t[i] = eta[i]-(eta_base[i]+e_m)-e_d*nb[i]
        while e_t[i]>=e_d:
            nb[i] = nb[i] + 1
            e_t[i] = eta[i]-(eta_base[i]+e_m)-e_d*nb[i]
            
    
    for k in np.arange(0,nk+1):
        for i in np.arange(1,nx+1):
            fak[k][i] = fk[k]
            ftk[k][i] = fk[k]
            for m in np.arange(0,int(nb[i])+1):
                fdk[m][k][i] = fk[k]
                
    RepresentativeDiameter(dm,dk,fak,nx,nk)
    
    hidingexposure(tsck,tscm,dm,dk,spec,g,nu,nx,nk)
    
    ManningN(snmm,dm,g,nx)
    
    snm = snmm[1]
        
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
        
        h[0] = (snmm[1]*qt/(wid*ib**0.5))**0.6
        u[0] = qt/(h[0]*wid)
        q[0] = qt/wid
        
        wu[0] = u[0]
        
        h[nx] = h[nx-1]
        
        momentum2(h,eta,u,wu,q,g,rdx,dt,snmm,nx)
        
        continuity2(h,q,rdx,dt,nx)
        
        u = wu
        
      #  wc = c
        
        # boundary condition
        
      
      #  usta = (g*(snmm[1]*u[0])**2./(h[0]**(1./3.)))**0.5
        
      #  if usta>wf:
      #      beta = 15.*wf/usta
      #      zet  = usta/wf*rep**0.6
      #      wc[0] = ent_GarciaParker(zet)/alf(beta)
      #  else:
      #      wc[0] = 0.
            
      
              
        # SscTransport(c,dd,ee,wc,u,h,q,g,snm,wf,rep,dx,dt,nx)
        
        bedload(h,u,qbk,qbin,snmm,spec,dk,tsck,g,nx,nk)
        
        for k in np.arange(0,nk+1):
            qbk[k][ 0] = qbk[k][1]
            qbk[k][nx] = qbk[k][nx-1] 
            
      
          
        
        if tt>morpho_start:
            
            exner(eta,qbk,dz,dzk,rdx,dt,poro,nx,nk)
            
            sorting(fak,ftk,fdk,dz,dzk,e_t,e_d,nb,fak_new,ftk_new,fdk_new,nx,nk,poro)
            
            RepresentativeDiameter(dm,dk,fak,nx,nk)
    
            hidingexposure(tsck,tscm,dm,dk,spec,g,nu,nx,nk)
            
            ManningN(snmm,dm,g,nx)
            
        if optime>tuk:
            
            print("Time= ",tt/3600, "hours")
            optime = optime-tuk
            
            hh = h+eta
            
            
            '''
            
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
            
            '''
            
            fig, ax1 = plt.subplots()
            ax2 = ax1.twinx()
            
            ax1.set_title("Time= {0:.2f} hours".format(tt/3600))
            ax1.set_xlim([0,chlen])
            ax1.set_xlabel( "Downstream distance (m)" ) 
            
            ax1.set_ylim([0,z00+hmax*2.])
            
            ax1.plot(xc[1:nx],eta[1:nx],color='k',label='Bed surface')
            ax1.plot(xc[1:nx],hh[1:nx],color='b',label='Water surface')
            ax1.set_ylabel( "Elevation (m)" )
            
         #   ax1.legend()
            
            ax2.plot(xc[1:nx],dm[1:nx]*1000,color='r',label='Mean diameter')
        #    ax2.plot(x[1:nx],qbk[0][1:nx],label='Mean diameter')
        #    ax2.plot(x[1:nx],qbk[1][1:nx],label='Mean diameter')
        #    ax2.plot(x[1:nx],qbk[2][1:nx],label='Mean diameter')
        #    ax2.set_ylim([dk[0]*1000,dk[nk]*1000])
            ax2.set_ylim([1.5,2.5])
            ax2.set_ylabel( "Diameter (mm)" )
        #    ax2.legend()
        
            handler1, label1 = ax1.get_legend_handles_labels()
            handler2, label2 = ax2.get_legend_handles_labels()
            
            ax1.legend(handler1 + handler2, label1 + label2)
            
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