# -*- coding: utf-8 -*-
"""
Created on Tue May 11 10:51:46 2021

@author: Zsuzsanna Koczor-Benda, UCL
"""
from __future__ import print_function
import math
import numpy as np
import matplotlib.pyplot as plt

# functions for calculating broadened spectrum for all input files
def lorentz0(res,gamma):
    rmax=gamma*30
    npoints=int(round(rmax/res))
    g=np.zeros((2*npoints))
    xs=[k*res for k in range(npoints)]
    for i,x in enumerate(xs):
        g[i]= 1/(math.pow((x),2) + math.pow(1/2*gamma, 2))
    return g

def displ_lorentz0(l0,x0,y,numpoints,xmin,res):
    disp=np.zeros((numpoints))
    lenl=len(l0)
    p0=int(round((x0-xmin)/res))
    if p0>numpoints+lenl or p0<-lenl:
        return disp
    if p0>numpoints:
        for p in range(numpoints+lenl-p0): 
            if numpoints-p-1 <0:
                break
            disp[numpoints-p-1]=y*l0[p0-numpoints+p]
    elif p0<0:
        for p in range(lenl+p0):  
            if p==numpoints:
                break
            disp[p]=y*l0[p-p0]
    else:
        for p in range(min(numpoints-p0,lenl)):  
            disp[p0+p]=y*l0[p]
        for p in range(min(p0,lenl)):            
            disp[p0-p-1]=y*l0[p]

    return disp
    
def calc_broadened_spectrum(freqs,rawints,xmin, xmax, res, gamma):
    l0=lorentz0(res,gamma)
    numpoints=int(round((xmax-xmin)/res))
    wn=np.zeros(numpoints)
    sp_ints=np.zeros((numpoints)) 
    x=xmin
    for i in range(0,numpoints):
        wn[i]=x
        x+=res      

    spectrum=np.zeros((numpoints))
    for state in range(0,len(freqs)):
        spectrum+=displ_lorentz0(l0,freqs[state],rawints[state],numpoints,xmin,res)
    sp_ints=spectrum/math.pi * 1/2 * gamma
    return wn, sp_ints

def calc_broadened_conv_spectrum(freqs,rawints1,xmin, xmax, res, gamma1,gamma2):
    l0=lorentz0(res,gamma1)
    l02=lorentz0(res,gamma2)
    numpoints=int(round((xmax-xmin)/res))
    wn=np.zeros(numpoints)
    sp_ints=np.zeros((numpoints)) 
    x=xmin
    for i in range(0,numpoints):
        wn[i]=x
        x+=res      

    spectrum=np.zeros((numpoints))
    for state in range(0,len(freqs)):
        ones=1 #np.full((np.shape(rawints1[f,state])),1)
            
        spectrum+=rawints1[state]*displ_lorentz0(l0,freqs[state],ones,numpoints,xmin,res)*displ_lorentz0(l02,freqs[state],ones,numpoints,xmin,res)
        sp_ints=spectrum/(math.pi**2) * 1/4 * gamma1*gamma2  
    return wn, sp_ints

# create broadened spectrum for single molecule
def create_average_spec_single(freqs,  IRints, Rints, Pints,xmin=30,xmax=1000,res=0.5,gammaIR=5,gammaR=5,sclf=0.98):
    
    # scale frequencies
    freqs=sclf*freqs
    fmin=xmin-200
    fmax=xmax+200 
    
    maxdof=len(freqs)
    R_fre=np.zeros((maxdof))
    IR_fre=np.zeros((maxdof))
    P_fre=np.zeros((maxdof))

    for l in range(0,maxdof):
        if(fmin<=freqs[l]<fmax):
                R_fre[l]=Rints[l]
                IR_fre[l]=IRints[l]
                P_fre[l]=Pints[l]
                
    prod_ints=np.zeros_like(IR_fre)
    R_ints=np.zeros_like(IR_fre)
    IR_ints=np.zeros_like(IR_fre)
    for f,fr in enumerate(freqs):
        if fr<xmin or fr>xmax:
            continue
        prod_ints[f]=P_fre[f] 
        R_ints[f]=R_fre[f]
        IR_ints[f]=IR_fre[f]
            
    # average all orientations
    wn, IR_spec = calc_broadened_spectrum(freqs,IR_fre,xmin, xmax, res, gammaIR)
    wn, R_spec = calc_broadened_spectrum(freqs,R_fre,xmin, xmax, res, gammaR)
    wn, conv_spec=calc_broadened_conv_spectrum(freqs,P_fre,xmin, xmax, res, gammaIR,gammaR)    
    return wn,R_spec,IR_spec,conv_spec,freqs,prod_ints,R_ints,IR_ints 

# plot polar plot of IR, Raman, conv value for initial molecular orientation
def polar_plot(theta,r):
    plt.rcParams.update({'font.size': 10})
    plt.close("Orientation dependence of intensity")
    fig, ax = plt.subplots(num="Orientation dependence of intensity",subplot_kw={'projection': 'polar'})
    fig.set_size_inches(3, 3)
    line=ax.plot(theta,r,linewidth=4.0,c='k')
    #ax.set_rmax(max(r))
    ax.set_rticks([])
    #ax.set_rlabel_position(45)  
    ax.grid(True)
    #ax.set_title("Polar plot IR, R, conv", va='bottom')
    #plt.show()
    return line
    