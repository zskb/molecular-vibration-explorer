# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 11:30:47 2021

@author: Zsuzsanna Koczor-Benda, UCL
"""
import numpy as np
                
def load_from_dat(filename):
    with open(filename) as inpfile:
        line = inpfile.readline()
        if "Number of atoms" in line:
            line = inpfile.readline()
            natoms=int(line.split()[0])
            line = inpfile.readline()
            line = inpfile.readline()
        if("Anisotropy" in line):
            line = inpfile.readline()
            aniso=float(line.split()[0])
            line = inpfile.readline()
        while line:       
            if "Frequencies" in line:
                fr_=[]
                line = inpfile.readline()
                m=0
                while "--" not in line:
                    fr_.append(float(line.split()[0]))
                    m+=1
                    line = inpfile.readline()
                fr=np.array(fr_)
                break
            line = inpfile.readline()
            
        nmodes=np.shape(fr)[0]
        Q=np.zeros((nmodes,natoms,3))
        while line:       
            if "Displacements" in line:
                
                line = inpfile.readline()
                line = inpfile.readline()
                for m in range(nmodes):
                    for a in range(natoms):
                        spl=line.split()
                        for al in range(3):
                            Q[m,a,al]=float(spl[al])
                        line = inpfile.readline()
                    line = inpfile.readline()
                break
            line = inpfile.readline() 
        D=np.zeros((nmodes,3))
        while line:       
            if "Dipole derivatives" in line:
                
                line = inpfile.readline()
                line = inpfile.readline()
                for m in range(nmodes):
                    spl=line.split()
                    for al in range(3):
                        D[m,al]=float(spl[al])
                    line = inpfile.readline()
                    line = inpfile.readline()
                break
            line = inpfile.readline()     
#              
        P=np.zeros((nmodes,3,3))
        while line: 
            if "Polarizability derivatives" in line:
                
                line = inpfile.readline()
                line = inpfile.readline()
                for m in range(nmodes):
                    for a in range(3):
                        spl=line.split()
                        for al in range(3):
                            P[m,a,al]=float(spl[al])
                        line = inpfile.readline()
                    line = inpfile.readline()
                break
            line = inpfile.readline()     
    return fr,Q,D,P,natoms,aniso

def load_data_file(filename):
    if filename.endswith(".dat"):
        try:
            fr,Q,D,P,nat,aniso=load_from_dat(filename)
       #     print("Loaded ",filename)
            return fr,Q,D,P,nat,aniso
        except:
            print("Error with loading file ",filename)
            return 


