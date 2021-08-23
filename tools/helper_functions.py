# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 11:34:25 2021

@author: Zsuzsanna Koczor-Benda, UCL
"""

import numpy as np
import math

phys_constants=dict(amu = 1.66053878e-27,    # kg
                    h   = 6.62606896e-34,    # Js
                    c   = 2.99792458e+8,     # m/s
                    k   = 1.3806504e-23      # J/K)
                    )
scaling_factors=dict(    
                    D=2.541746473,  # for dipole moment: 1 ea0 =2.54 Debye
                    A=0.529177210903, # 1 bohr = 0.529 Angstrom
                    IRfac=126.8 # factor from Philippe's paper for D2/A2amu to km/mol
                    )

def calc_scaling(T):
    # calculate scaling factors 
    #  R[cm^4/kg] = 10^-32/amu[kg] * R[Angstrom^4/amu]
    #  c[cm/s]=10^2 * c[m/s] 
    #  factor of 10^4 needed to agree with Principles of SERS book Eq. A.10
    amu=phys_constants['amu']
    h=phys_constants['h']
    c=phys_constants['c']
    k=phys_constants['k']
    scalingR=math.pow(10, 4)*h*math.pow(math.pi, 2)/((math.pow(10, 32)*amu)*(math.pow(10, 2)*c)*22.5)
    scalingR=45*math.pow(scaling_factors['A'],4)*scalingR
    scalingpolar=math.pow(scaling_factors['A'],4)
    scalingexp=-h*(math.pow(10, 2)*c)/(k*T)
    scalingIR=math.pow(scaling_factors['D']/scaling_factors['A'],2)*scaling_factors['IRfac']
    return scalingIR,scalingR,scalingexp, scalingpolar

def single_rot(p,a=0,b=0,c=0):
    Rz=np.array([[math.cos(a), -math.sin(a),0],[math.sin(a),math.cos(a),0],[0,0,1]])
    Rx=np.array([[1,0,0],[0,math.cos(b), -math.sin(b)],[0,math.sin(b),math.cos(b)]])
    Rz2=np.array([[math.cos(c), -math.sin(c),0],[math.sin(c),math.cos(c),0],[0,0,1]])
    R=np.matmul(Rz2,np.matmul(Rx,Rz))
    rot=np.matmul(p,np.transpose(R))
    return rot

def single_rot_T(p,a=0,b=0,c=0):
    Rz=np.array([[math.cos(a), -math.sin(a),0],[math.sin(a),math.cos(a),0],[0,0,1]])
    Rx=np.array([[1,0,0],[0,math.cos(b), -math.sin(b)],[0,math.sin(b),math.cos(b)]])
    Rz2=np.array([[math.cos(c), -math.sin(c),0],[math.sin(c),math.cos(c),0],[0,0,1]])
    R=np.matmul(Rz2,np.matmul(Rx,Rz))
    rot=np.matmul(p,R)
    return rot

def unit_rot(axis,rad):
    if axis==1 :
        R=np.array([[math.cos(rad),-math.sin(rad),0],[math.sin(rad),math.cos(rad),0],[0,0,1]]) 
    elif axis==2 :
        R=np.array([[math.cos(rad),0,math.sin(rad)],[0,1,0],[-math.sin(rad),0,math.cos(rad)]])
    elif axis==3 :
        R=np.array([[1,0,0],[0,math.cos(rad),-math.sin(rad)],[0,math.sin(rad),math.cos(rad)]])
    else :
        R=np.array([[1,0,0],[0,1,0],[0,0,1]])
    return R
    
def get_oriented_box(coords,vdw):
        coormin=np.zeros_like(coords)
        coormax=np.zeros_like(coords)
        for i in range(len(vdw)): 
            coormin[i,:]=coords[i,:]-vdw[i]
            coormax[i,:]=coords[i,:]+vdw[i]
        cell=np.max(coormax,axis=0)-np.min(coormin,axis=0)
        return cell    

VdwRad={'H': 1.1 ,
'He': 1.4 ,
'Li': 1.81 ,
'Be': 1.53 ,
'B': 1.92 ,
'C': 1.7 ,
'N': 1.55 ,
'O': 1.52 ,
'F': 1.47 ,
'Ne': 1.54 ,
'Na': 2.27 ,
'Mg': 1.73 ,
'Al': 1.84 ,
'Si': 2.1 ,
'P': 1.8 ,
'S': 1.8 ,
'Cl': 1.75 ,
'Ar': 1.88 ,
'K': 2.75 ,
'Ca': 2.31 ,
'Sc': 2.3 ,
'Ti': 2.15 ,
'V': 2.05 ,
'Cr': 2.05 ,
'Mn': 2.05 ,
'Fe': 2.05 ,
'Co': 2.0 ,
'Ni': 2.0 ,
'Cu': 2.0 ,
'Zn': 2.1 ,
'Ga': 1.87 ,
'Ge': 2.11 ,
'As': 1.85 ,
'Se': 1.9 ,
'Br': 1.83 ,
'Kr': 2.02 ,
'Rb': 3.03 ,
'Sr': 2.49 ,
'Y': 2.4 ,
'Zr': 2.3 ,
'Nb': 2.15 ,
'Mo': 2.1 ,
'Tc': 2.05 ,
'Ru': 2.05 ,
'Rh': 2.0 ,
'Pd': 2.05 ,
'Ag': 2.1 ,
'Cd': 2.2 ,
'In': 2.2 ,
'Sn': 1.93 ,
'Sb': 2.17 ,
'Te': 2.06 ,
'I': 1.98 ,
'Xe': 2.16 ,
'Cs': 3.43 ,
'Ba': 2.68 ,
'La': 2.5 ,
'Ce': 2.48 ,
'Pr': 2.47 ,
'Nd': 2.45 ,
'Pm': 2.43 ,
'Sm': 2.42 ,
'Eu': 2.4 ,
'Gd': 2.38 ,
'Tb': 2.37 ,
'Dy': 2.35 ,
'Ho': 2.33 ,
'Er': 2.32 ,
'Tm': 2.3 ,
'Yb': 2.28 ,
'Lu': 2.27 ,
'Hf': 2.25 ,
'Ta': 2.2 ,
'W': 2.1 ,
'Re': 2.05 ,
'Os': 2.0 ,
'Ir': 2.0 ,
'Pt': 2.05 ,
'Au': 2.1 }

def rotate_molecule(atoms,refcoords,filename,phi=0,theta=0,xi=0):  
    
    # get Van der Waals radii
    vdw=np.zeros_like(atoms,dtype=float)
    for i,a in enumerate(atoms): 
        vdw[i]= VdwRad[a]

    # do custom rotation with phi and theta
    rotated=single_rot(refcoords,a=phi,b=theta,c=xi)

    # get dimensions of cell encapsulating the molecule
    rotc=np.zeros_like(refcoords)
    for s,c in enumerate(refcoords):
        if atoms[s]=='Au':
            vdw[s]=0     # do not count gold atom for cell size
            rotc[s]=np.array([0,0,0])
            continue 
        rotc[s]=rotated[s]
    cell=get_oriented_box(rotc,vdw) 

    # rotate axes instead of mol
    axes=np.eye(3)
    rot_axes=single_rot_T(axes,a=phi,b=theta,c=xi)

    return rot_axes,cell
