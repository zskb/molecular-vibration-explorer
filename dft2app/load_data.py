# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 11:30:47 2021

@author: Zsuzsanna Koczor-Benda, UCL
"""
import numpy as np
import math


def symmetrize(a):
    return a + a.T - np.diag(a.diagonal())

# load Hessian from g09 fchk
# to be diagonalized to get normal modes
def load_fchk(filename):
    Z=[]
    pol=np.zeros((3,3))
    with open(filename) as inpfile:
        line = inpfile.readline()
        while line:
            if "Atomic numbers" in line:
                natoms=int(line.split()[4])
                line = inpfile.readline()
                while "Nuclear charges" not in line:
                    spl=line.split()
                    for charge in spl:
                        Z.append(charge) 
                    line = inpfile.readline()
            if "Real atomic weights" in line:
                W=np.zeros((natoms))
                line = inpfile.readline()
                a=0
                while "Atom fragment" not in line:
                    spl=line.split()
                    for weight in spl:
                        W[a]=weight 
                        a+=1
                    line = inpfile.readline()
                break
            line = inpfile.readline()
        H=np.zeros((natoms*3,natoms*3))
        D=np.zeros((natoms,3,3))
        P=np.zeros((natoms,3,3,3))
        i=0
        j=0
        found_der=0
        while line:
            if "Cartesian Force Constants" in line:
                found_der=1
                line = inpfile.readline()             
                while "Dipole Moment" not in line:
                    spl=line.split()
                    for h in spl:
                        H[i,j]=float(h)
                        j+=1
                        if j>i:
                            j=0
                            i+=1
                    line = inpfile.readline()   
                    if "Nonadiabatic" in line:
                        while "Dipole Moment" not in line:
                            line = inpfile.readline() 
                        break
                break
            line = inpfile.readline()
        if found_der==0:
            return 0
        line = inpfile.readline() 
        line = inpfile.readline() 
        i=0
        j=0
        a=0
     #       if "Dipole Derivatives" in line:
        line = inpfile.readline()             
        while "Polarizability" not in line:
                    spl=line.split()
                    for h in spl:
                        D[i,a,j]=float(h)
                        j+=1
                        if j>=3:
                            j=0
                            a+=1
                        if a>=3:
                            a=0
                            i+=1
                    line = inpfile.readline()  
                    
        
        line = inpfile.readline() 
        spl=line.split()
        pol[0,0]=spl[0]
        pol[0,1]=spl[1]
        pol[1,1]=spl[2]
        pol[0,2]=spl[3]
        pol[1,2]=spl[4]
        line = inpfile.readline()
        spl=line.split()
        pol[2,2]=spl[0]
        anisotr=math.pow(2,-0.5)*math.pow((pol[0,0]-pol[1,1])**2+(pol[1,1]-pol[2,2])**2+(pol[2,2]-pol[0,0])**2+
                         6*(pol[0,1]**2+pol[1,2]**2+pol[0,2]**2),0.5) 
        line = inpfile.readline() 
        line = inpfile.readline() 
        i=0
        j=0
        k=0
        l=0
        while "HyperPolarizability" not in line:
                    spl=line.split()
                    for h in spl:
                        P[i,j,k,l]=float(h)
                        l+=1
                        if l>k:
                            l=0
                            k+=1
                        if k>=3:
                            k=0
                            j+=1
                        if j>=3:
                            j=0
                            i+=1
                    line = inpfile.readline()      
    for i in range(3*natoms):
        for j in range(i):
            H[j,i]=H[i,j]
    
    return H,Z,W,D,P,natoms,anisotr


def load_data_low_acc(filename):     
    Z=[]
    with open(filename) as inpfile:
        line = inpfile.readline()
        while line:
            if "NAtoms= " in line:
                natoms=int(line.split()[1])
                nmode=3*natoms-6
                
            if "Input orientation:" in line:
                line = inpfile.readline()
                line = inpfile.readline()
                line = inpfile.readline()
                line = inpfile.readline()
                line = inpfile.readline()
                while "---" not in line:
                    spl=line.split()
                    Z.append(spl[1])
                    line = inpfile.readline()
                break
            line = inpfile.readline()
    freqs=[]
    red=[]
    displ=np.zeros((nmode,natoms,3))
    
    with open(filename) as inpfile:
        line = inpfile.readline()
        pmode=0
        while line:           
            if " Frequencies" in line:  
                spl=line.split()
                nfreq=len(spl)-2
                for i in range(2,nfreq+2):
                    freqs.append(float(spl[i]))
                line = inpfile.readline()
                spl=line.split()  
                for i in range(3,nfreq+3):
                    red.append(float(spl[i])) #reduced masses
                while "   1  " not in line:
                    line = inpfile.readline()
                spl=line.split()
                at=0
                while line:
                    for a in range(3):
                        displ[pmode,at,a]=float(spl[a+2])
                        displ[pmode+1,at,a]=float(spl[a+5])
                        displ[pmode+2,at,a]=float(spl[a+8])
                    at+=1
                    line = inpfile.readline()
                    spl=line.split()
                    if len(spl)<=3 :
                        break
                pmode=pmode+3
            line = inpfile.readline()
    return np.array(freqs),displ,Z,np.array(red)

def load_deriv_data(filename,natoms):     
    foot=""
    P=""
    with open(filename) as inpfile:
        line = inpfile.readline()
        while line:
            if "\\Freq\\RB3LYP\\" in line:
                foot=foot+line.strip()
                while line:
                    line = inpfile.readline()
                    foot=foot+line.strip()
                    if "\@" in line:
                        break
                
                break
            line = inpfile.readline()
    spl=foot.split("\\")
    for s in spl:
        if "DipoleDeriv" in s:
            DD=(s.split("=")[1]).split(",")
        elif "Polar=" in s and "Hyper" not in s:
            PP=(s.split("=")[1]).split(",")
        elif "PolarDeriv" in s:
            PD=(s.split("=")[1]).split(",")
    i=0
    j=0
    a=0
    D=np.zeros((natoms,3,3))
    P=np.zeros((natoms,3,3,3))
    pol=np.zeros((3,3))
    aniso=0
    for d in DD:
        D[i,a,j]=float(d)
        j+=1
        if j>=3:
            j=0
            a+=1
        if a>=3:
            a=0
            i+=1
    i=0
    j=0
    k=0
    l=0
    for p in PD:
        P[i,j,k,l]=float(p)
        l+=1
        if l>k:
             l=0
             k+=1
        if k>=3:
             k=0
             j+=1
        if j>=3:
             j=0
             i+=1
    i=0
    j=0
    for p in PP:
        pol[j,i]=float(p)
        j+=1
        if j>i:
            j=0
            i+=1
    aniso=math.pow(2,-0.5)*math.pow((pol[0,0]-pol[1,1])**2+(pol[1,1]-pol[2,2])**2+(pol[2,2]-pol[0,0])**2+
                         6*(pol[0,1]**2+pol[1,2]**2+pol[0,2]**2),0.5)
    return D,P,aniso

# load Hessian from Gaussian fchk and diagonalize                
def load_from_Hessian(filename):
    try:
        H,Z,W,D,P,nat,aniso=load_fchk(filename)  
    #    print(D)
    except:
        print("Error loading fchk")
        return 0
    Hw=np.zeros_like(H)
    # in fchk: lower triangle of H
    for i in range(nat*3):
        for j in range(nat*3):
            a=i//3
            b=j//3
            Hw[i,j]=H[i,j]/np.sqrt(W[a]*W[b]) # do mass-weighting
    w, v = np.linalg.eig(Hw)
    ws=np.sort(w)
    sort=np.argsort(w)
    vs=v[:,sort]
    nm=nat*3-6 
    vw=np.zeros((nm,nat,3))
    for mode in range(0,nat*3): 
        if mode<6: # remove first 6 modes: translation+rotation
            continue
        for i in range(nat):
            for alpha in range(3):
                vw[mode-6,i,alpha]=vs[i*3+alpha,mode] # these are already mass weighted, G09 output is not  
    ww=np.sqrt(ws[6:]) # remove first 6 modes
    
    norm=np.zeros((nm))
    redm=np.zeros((nm))
    Qcart=np.zeros_like(vw)
    for m in range(nm):
        for at in range(nat):
            for alpha in range(3):
                Qcart[m,at,alpha]=vw[m,at,alpha]/np.sqrt(W[at])
                norm[m]+=np.square(Qcart[m,at,alpha])
        redm[m]=1/norm[m]
        norm[m]=1/np.sqrt(norm[m])
        Qcart[m,:,:]=norm[m]*Qcart[m,:,:]
  #  print("redm ",redm)  
   # print(Qcart)
    D_normal=np.zeros((nm,3))
    P_normal=np.zeros((nm,3,3))
  #  print("P ",P)
    Psym=np.zeros((nm,3,3))
    for m in range(nm):
        for i in range(nat):
            for al in range(3):
                    D_normal[m,:]+=D[i,al,:]*Qcart[m,i,al]/np.sqrt(redm[m])#/np.sqrt(W[i])
                    P_normal[m,:,:]+=P[i,al,:,:]*Qcart[m,i,al]/np.sqrt(redm[m])
        Psym[m]=symmetrize(P_normal[m,:,:])
    sclf=5140.487127157137 # get frequencies in cm-1 
    return ww*sclf,vw,Z,W,D_normal,Psym,nat,aniso    
    
# load data from Gaussian output             
def load_from_output(filename):
    try:
        fr,displ,Z,redm=load_data_low_acc(filename)
        #  print("loaded .out")
        nat=len(Z)
        D,P,aniso=load_deriv_data(filename,nat)  
  #      print(aniso)
    except:
        print("Error loading output")
        return 0
    nm=len(fr)
 #   print("modes ",nm)
    Qcart=displ
  #  print("redm ",redm)        
    D_normal=np.zeros((nm,3))
    P_normal=np.zeros((nm,3,3))
  #  print("P ",P)
    Psym=np.zeros((nm,3,3))
    for m in range(nm):
        for i in range(nat):
            for al in range(3):
                    D_normal[m,:]+=D[i,al,:]*Qcart[m,i,al]/np.sqrt(redm[m])#/np.sqrt(W[i])
                    P_normal[m,:,:]+=P[i,al,:,:]*Qcart[m,i,al]/np.sqrt(redm[m])
    
        Psym[m]=symmetrize(P_normal[m,:,:])
    return fr,Z,Qcart,D_normal,Psym,nat,aniso 

def load_file(filename):
    # read from fchk more accurate
    if filename.endswith(".fchk"):
        try:
            fr,Q,Z,W,D,P,nat,aniso=load_from_Hessian(filename)
            print("Loaded ",filename)
            return fr,Z,Q,D,P,nat,aniso
        except:
            print("Error with loading file ",filename)
            return 
    # read from output less accurate            
    if filename.endswith(".out"):
        try:
            fr,Z,Q,D,P,nat,aniso =load_from_output(filename) 
            print("Loaded ",filename)
            return fr,Z,Q,D,P,nat,aniso
        except:
            print("Error with loading file ",filename)
            return
        


