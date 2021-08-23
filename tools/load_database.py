# -*- coding: utf-8 -*-
"""
Created on Tue May 11 10:51:46 2021

@author: Zsuzsanna Koczor-Benda, UCL
"""
from __future__ import print_function
import numpy as np

def read_in_average_ints_all(filename):
    fname=[]
    smiles=[]
    ints=[]
    with open(filename) as inpfile:
        line=inpfile.readline()
        line=inpfile.readline()
        while line:
            spl=line.split()
            fname.append(spl[0])
            smiles.append(spl[1]) 
            ints.append([float(vib) for vib in spl[2:]])
            line=inpfile.readline()
    ints_np=np.array(ints)
    return ints_np,fname,smiles

# load intensities and get target properties
def get_intensities_all(sclf,SH):

    if SH:
        freqs,fname,smiles=read_in_average_ints_all("data_SH/frequencies_SH.txt")
        IR,fname,smiles=read_in_average_ints_all("data_SH/IR_intensities_SH.txt")
        Raa,fname,smiles=read_in_average_ints_all("data_SH/Raman_Stokes_intensities_aa_SH.txt")
        Rab,fname,smiles=read_in_average_ints_all("data_SH/Raman_Stokes_intensities_ab_SH.txt")
        Paaa,fname,smiles=read_in_average_ints_all("data_SH/Conversion_anti-Stokes_intensities_aaa_SH.txt")
        Pbaa,fname,smiles=read_in_average_ints_all("data_SH/Conversion_anti-Stokes_intensities_baa_SH.txt")
        Paab,fname,smiles=read_in_average_ints_all("data_SH/Conversion_anti-Stokes_intensities_aab_SH.txt")
        Pabc,fname,smiles=read_in_average_ints_all("data_SH/Conversion_anti-Stokes_intensities_abc_SH.txt")
    else:
        freqs,fname,smiles=read_in_average_ints_all("data_SAu/frequencies.txt")
        IR,fname,smiles=read_in_average_ints_all("data_SAu/IR_intensities.txt")
        Raa,fname,smiles=read_in_average_ints_all("data_SAu/Raman_Stokes_intensities_aa.txt")
        Rab,fname,smiles=read_in_average_ints_all("data_SAu/Raman_Stokes_intensities_ab.txt")
        Paaa,fname,smiles=read_in_average_ints_all("data_SAu/Conversion_anti-Stokes_intensities_aaa.txt")
        Pbaa,fname,smiles=read_in_average_ints_all("data_SAu/Conversion_anti-Stokes_intensities_baa.txt")
        Paab,fname,smiles=read_in_average_ints_all("data_SAu/Conversion_anti-Stokes_intensities_aab.txt")
        Pabc,fname,smiles=read_in_average_ints_all("data_SAu/Conversion_anti-Stokes_intensities_abc.txt")
                
    # scale frequencies
    freqs=sclf*freqs
    
    return freqs,IR,Raa,Rab,Paaa,Pbaa,Paab,Pabc,smiles,fname

# get intensities in custom frequency range
def get_intens_range(freqs,intens,fmin,fmax):
    intens_range=np.zeros_like(freqs)
    for f in range(0,len(freqs)):
        for mode in range(0,len(freqs[f])):
            if (freqs[f,mode]>= fmin) and (freqs[f,mode]<=fmax):
                intens_range[f,mode]=intens[f,mode]
    return intens_range

# calculate target quantities P, A, or R  
# means and stds from 1.3k randomly selected thiols
def get_target(intens,target_type="P"):    
    if target_type=="P": # conversion
        Tmean=-28.367672289689146 
        Tstd=0.3809242121851108
    elif target_type=="R": # Raman Stokes
        Tmean=-28.91883480924085
        Tstd=0.2625275063027046  
    else: # IR absorption
        Tmean=2.290063042859388
        Tstd=0.22266072753602675
    sint=np.sum(intens,axis=1)
    target_=np.log10(sint, out=np.full_like(sint,np.NINF),where=sint!=0)    
    target=(target_-Tmean)/Tstd
    return target




