# -*- coding: utf-8 -*-
"""
Analysis code for hydrological statistics

@author: Toshiki Iwasai @ Hokkaido University
"""

# Include some libraries 

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm # normal distribution
from scipy import optimize
from statistics import mean, median, variance, stdev

def func(parameter,x,y):
    a = parameter[0]
    b = parameter[1]
    residual = y-(a*x+b)
    return residual

# define work directory and name of input text file 

path    = "C:/work/class/2019/wrr/rainfall/"
fname   = "kadai"

# open input file, which includes time-precipitation (mm/h) 

with open(path+fname+".txt","r") as fpin:
    lines = fpin.readlines()
    tt = []
    rr = []
    for line in lines:
        value = line.replace('\n', '').split('\t')
        tt.append(float(value[0]))
        rr.append(float(value[1]))
    
fpin.close()

# plot time series of year-maximum hourly-precipitation rate 

plt.plot( tt, rr )
plt.xlabel( "Year" )
plt.ylabel( "Precipitation (mm/h)" )
plt.xlim([1970,2020])
plt.ylim([0,70])
plt.savefig(path+fname+"_tt-rr.jpg",dpi=300)
plt.close()

# make the x as log, and calculate mean, variance, and standard deviation 

rrl = np.log10(rr)
ave_rrl = median(rrl)
var_rrl = variance(rrl)
std_rrl = stdev(rrl)

# plot histgram of the data as log scale 

plt.hist(rrl,bins=10,range=[0.9,1.8],normed=True,histtype='step')
plt.xlabel( "Precipitation (mm/h)" )
plt.ylabel( "Probability" )
#plt.xlim([0,70])
plt.savefig(path+fname+"_hist-log.jpg",dpi=300)
plt.close()

# Find eq of log-normal probability distribution 
 
pdf = []

p1 = list(range(1,70))

for i in p1:
    wpdf = np.exp(-0.5*((np.log10(float(i))-ave_rrl)/std_rrl)**2)/(std_rrl*float(i)*(2*3.14)**0.5)
    pdf.append(wpdf)

# plot histgram on normal x coordinate wit log-normal PDF 
    
plt.hist(rr,bins=15,range=[0,70],normed=True,histtype='step')
plt.plot(p1,pdf)
plt.xlabel( "Precipitation (mm/h)" )
plt.ylabel( "Probability" )
plt.xlim([0,70])
plt.savefig(path+fname+"_hist.jpg",dpi=300)
plt.close()

# ploting position 

rs = sorted(rr,reverse=True)
nn = len(rs)+1
p = list(range(1,nn))

pp = []
ff = []

for i in p:
    wp = (2*float(i)-1)/float(2*nn)     # Hazen plot
    wf = 1-wp
    pp.append(wp)
    ff.append(wf)
    
plt.figure(figsize=(3.5,4))

ppp = norm.ppf( ff, loc=0, scale=1 )
qqq = np.log10(rs)

# regression analysis 

parameter0 = [0.0,0.0]
result = optimize.leastsq(func,parameter0,args=(ppp,qqq))
aa=result[0][0]
bb=result[0][1]

# plot the data on the log-normal probability paper

xmin=np.log10(1)
xmax=np.log10(100)
ymin=norm.ppf(0.001, loc=0, scale=1)
ymax=norm.ppf(0.999, loc=0, scale=1)
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.tick_params(labelbottom=False)
plt.tick_params(labelleft=False)
plt.tick_params(which='both', width=0)

plt.plot(qqq,ppp,'o',color='k',markersize=2)

plt.plot([xmin, xmax], [(xmin-bb)/aa, (xmax-bb)/aa], color='k', linestyle='-', linewidth=1)

# setting of x-y axes
_dy=np.array([0.001,0.01,0.1,0.5,0.9,0.99,0.999])
dy=norm.ppf(_dy, loc=0, scale=1)
plt.hlines(dy, xmin, xmax, color='grey',linewidth=0.5)
#_dx=np.array([1,10,100])
_dx=np.array([1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100])
dx=np.log10(_dx)
plt.vlines(dx, ymin, ymax, color='grey',linewidth=0.5)

fs=10
for i in range(1,4):
    plt.text(float(i), ymin-0.1, str(10**i), ha = 'center', va = 'top', fontsize=fs)

for i in range(0,7):
    plt.text(xmin-0.01, dy[i], str(_dy[i]), ha = 'right', va = 'center', fontsize=fs)
    
plt.text(0.5*(xmin+xmax), ymin-0.5, 'Rainfall (mm/h)', ha = 'center', va = 'center', fontsize=fs)
plt.text(xmin-0.25,0.5*(ymin+ymax),'Non-exceedance probability', ha = 'center', va = 'center', fontsize=fs, rotation=90)


plt.savefig(path+fname+"_pp-rr.jpg",dpi=300)
plt.close()

