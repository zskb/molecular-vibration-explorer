# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 18:32:37 2022

@author: Zsuzsanna Koczor-Benda

Extract data from Gaussian formatted checkpoint files (.fchk) or output files (.out)
for use in Molecular Vibration Explorer

"""

###########################{libraries}##################################

import math
import os
import numpy as np
from timeit import default_timer as timer
import copy
from load_data import load_file
from cairosvg import svg2png
from rdkit import Chem # optional, only for having a unified SMILES format across the database and creating 2D molecule depictions
from rdkit.Chem import Draw


#######################{user editable variables}########################


# reading from formatted checkpoint file (.fchk) is more accurate than from output file (.out)
# to generate .fchk file from .chk use the formchk command of Gaussian (see https://gaussian.com/formchk/)
read_from_fchk=False

file_ext=".fchk"
if read_from_fchk==False:
    file_ext=".out" 
    
# directory where .fchk or .out files are located
# format for Windows:
dirname="C:\\\\your-local-path\\\\test-molecules"
# filenames are presumed to have format "freq-0003921.out" or "freq-0003921.fchk"
# where the number is an identifier of the molecule (can contain numbers/letters/symbols valid for filenames, except '.')

# we need a file called "smiles_data.txt" that contains SMILES code for each filename, see example
# use online databases/tools or molecule editors to get SMILES, generating SMILES directly from output is not reliable
smilesfile="test-molecules\\\\smiles_data.txt" 

maxdof=501 #need to set the maximum degrees of freedom in the database


#######################{script constants}###############################


phys_constants=dict(amu = 1.66053878e-27,    # kg
                    h   = 6.62606896e-34,    # Js
                    c   = 2.99792458e+8,     # m/s
                    k   = 1.3806504e-23      # J/K)
                    )
scaling_factors=dict(    # J/K)
                    D=2.541746473,  # for dipole moment: 1 ea0 =2.54 Debye
                    A=0.529177210903, # 1 bohr = 0.529 Angstrom
                    IRfac=126.8 # for D2/A2amu to km/mol
                    )

phys_param=dict(laser =785, # wavelength in nm 
                T = 298.15 
                )


###########################{functions}##################################


def calc_scaling(): 
    amu=phys_constants['amu']
    h=phys_constants['h']
    c=phys_constants['c']
    k=phys_constants['k']
    T=phys_param['T']
    scaling=math.pow(10, 4)*h*math.pow(math.pi, 2)/((math.pow(10, 32)*amu)*(math.pow(10, 2)*c)*22.5)
    scalingexp=-h*(math.pow(10, 2)*c)/(k*T)
    return scaling,scalingexp


def write_int_file(filnam, smiles, filename,ints,found,maxdof):
    ints_np=np.zeros((len(ints),maxdof))
    for i in range(0,len(ints)):
        for j in range(0,len(ints[i])):
            ints_np[i,j]=ints[i][j]
    with open(filnam,"w") as ofile:
        ofile.write("code smiles normal_mode_intensities\n".format())
        for m in range(len(filename)):
            if found[m]==0:
                continue
            ofile.write("{} {} ".format(filename[m],smiles[m]))
            for d in range(np.shape(ints_np)[1]):
                ofile.write("{:.5e} ".format(ints_np[m,d]))
            ofile.write("\n")
            
def write_files_pol(fr, ir_av,r_aa,r_ab,p_aaa,p_baa,p_aab,p_abc, smiles, filename,found):
    write_int_file("frequencies.txt", smiles, filename,fr,found,maxdof)
    write_int_file("IR_intensities.txt", smiles, filename,ir_av,found,maxdof)
    write_int_file("Raman_Stokes_intensities_aa.txt", smiles, filename,r_aa,found,maxdof)
    write_int_file("Raman_Stokes_intensities_ab.txt", smiles, filename,r_ab,found,maxdof)
    write_int_file("Conversion_anti-Stokes_intensities_aaa.txt", smiles, filename,p_aaa,found,maxdof)
    write_int_file("Conversion_anti-Stokes_intensities_baa.txt", smiles, filename,p_baa,found,maxdof)
    write_int_file("Conversion_anti-Stokes_intensities_aab.txt", smiles, filename,p_aab,found,maxdof)
    write_int_file("Conversion_anti-Stokes_intensities_abc.txt", smiles, filename,p_abc,found,maxdof)
    print("Data files written")
    
def symmetrize(a):
    return a + a.T - np.diag(a.diagonal())

def full_average_IR(d):
    full_av=(math.pow(d[0],2) + math.pow(d[1],2) + math.pow(d[2],2))/3
    return full_av

def full_average_R_pol(p,e_in,e_out):
    w1,w2,w3=e_in
    x1,x2,x3=e_out
    a11=p[0,0]
    a12=p[0,1]
    a13=p[0,2]
    a22=p[1,1]
    a23=p[1,2]
    a33=p[2,2]
    aux0=(6.*(w1*(x1*((w2*x2)+(w3*x3)))))+((w1**2)*(((2.*(x1**2))-(x3**2))-(x2**2)));
    aux1=((6.*(w2*(w3*(x2*x3))))+aux0)-((w2**2)*((x1**2)+((-2.*(x2**2))+(x3**2))));
    aux2=((w2**2)*((x1**2)+((3.*(x2**2))+(x3**2))))+((w3**2)*((x1**2)+((x2**2)+(3.*(x3**2)))));
    aux3=(4.*(w1*(x1*((w2*x2)+(w3*x3)))))+(((w1**2)*((3.*(x1**2))+((x2**2)+(x3**2))))+aux2);
    aux4=((w3**2)*((3.*(x1**2))+((3.*(x2**2))+(4.*(x3**2)))))+((w1**2)*((4.*(x1**2))+(3.*((x2**2)+(x3**2)))));
    aux5=(2.*(w1*(x1*((w2*x2)+(w3*x3)))))+(((w2**2)*((3.*(x1**2))+((4.*(x2**2))+(3.*(x3**2)))))+aux4);
    aux6=((a11**2)*((4.*(w2*(w3*(x2*x3))))+aux3))+((a12**2)*((2.*(w2*(w3*(x2*x3))))+aux5));
    aux7=(a11*((a22+a33)*(aux1-((w3**2)*((x1**2)+((x2**2)+(-2.*(x3**2))))))))+aux6;
    aux8=(2.*(a22*(a33*((w3**2)*(x3**2)))))+((3.*((a33**2)*((w3**2)*(x3**2))))+aux7);
    aux9=(3.*((a22**2)*((w3**2)*(x3**2))))+((4.*((a23**2)*((w3**2)*(x3**2))))+aux8);
    aux10=(3.*((a23**2)*((w2**2)*(x3**2))))+(((a33**2)*((w2**2)*(x3**2)))+((4.*((a13**2)*((w3**2)*(x3**2))))+aux9));
    aux11=((a33**2)*((w1**2)*(x3**2)))+((3.*((a13**2)*((w2**2)*(x3**2))))+(((a22**2)*((w2**2)*(x3**2)))+aux10));
    aux12=(3.*((a13**2)*((w1**2)*(x3**2))))+(((a22**2)*((w1**2)*(x3**2)))+((3.*((a23**2)*((w1**2)*(x3**2))))+aux11));
    aux13=(6.*(a22*(a33*(w2*(w3*(x2*x3))))))+((4.*((a33**2)*(w2*(w3*(x2*x3)))))+aux12);
    aux14=(4.*((a22**2)*(w2*(w3*(x2*x3)))))+((2.*((a23**2)*(w2*(w3*(x2*x3)))))+aux13);
    aux15=(4.*((a33**2)*(w1*(w3*(x1*x3)))))+((2.*((a13**2)*(w2*(w3*(x2*x3)))))+aux14);
    aux16=(2.*((a23**2)*(w1*(w3*(x1*x3)))))+((6.*(a22*(a33*(w1*(w3*(x1*x3))))))+aux15);
    aux17=(2.*((a13**2)*(w1*(w3*(x1*x3)))))+((4.*((a22**2)*(w1*(w3*(x1*x3)))))+aux16);
    aux18=((a22**2)*((w3**2)*(x2**2)))+((3.*((a23**2)*((w3**2)*(x2**2))))+(((a33**2)*((w3**2)*(x2**2)))+aux17));
    aux19=(3.*((a33**2)*((w2**2)*(x2**2))))+((3.*((a13**2)*((w3**2)*(x2**2))))+aux18);
    aux20=(4.*((a23**2)*((w2**2)*(x2**2))))+((2.*(a22*(a33*((w2**2)*(x2**2)))))+aux19);
    aux21=(4.*((a13**2)*((w2**2)*(x2**2))))+((3.*((a22**2)*((w2**2)*(x2**2))))+aux20);
    aux22=((a22**2)*((w1**2)*(x2**2)))+((3.*((a23**2)*((w1**2)*(x2**2))))+(((a33**2)*((w1**2)*(x2**2)))+aux21));
    aux23=(4.*((a33**2)*(w1*(w2*(x1*x2)))))+((3.*((a13**2)*((w1**2)*(x2**2))))+aux22);
    aux24=(2.*((a23**2)*(w1*(w2*(x1*x2)))))+((6.*(a22*(a33*(w1*(w2*(x1*x2))))))+aux23);
    aux25=(2.*((a13**2)*(w1*(w2*(x1*x2)))))+((4.*((a22**2)*(w1*(w2*(x1*x2)))))+aux24);
    aux26=((a22**2)*((w3**2)*(x1**2)))+((3.*((a23**2)*((w3**2)*(x1**2))))+(((a33**2)*((w3**2)*(x1**2)))+aux25));
    aux27=(3.*((a23**2)*((w2**2)*(x1**2))))+(((a33**2)*((w2**2)*(x1**2)))+((3.*((a13**2)*((w3**2)*(x1**2))))+aux26));
    aux28=(3.*((a33**2)*((w1**2)*(x1**2))))+((3.*((a13**2)*((w2**2)*(x1**2))))+(((a22**2)*((w2**2)*(x1**2)))+aux27));
    aux29=(4.*((a23**2)*((w1**2)*(x1**2))))+((2.*(a22*(a33*((w1**2)*(x1**2)))))+aux28);
    aux30=(4.*((a13**2)*((w1**2)*(x1**2))))+((3.*((a22**2)*((w1**2)*(x1**2))))+aux29);
    aux31=((aux30-(a22*(a33*((w2**2)*(x3**2)))))-(a22*(a33*((w1**2)*(x3**2)))))-(a22*(a33*((w3**2)*(x2**2))));
    aux32=((aux31-(a22*(a33*((w1**2)*(x2**2)))))-(a22*(a33*((w3**2)*(x1**2)))))-(a22*(a33*((w2**2)*(x1**2))));
    output=0.0666667*aux32;
    return output

def full_average_conv_aaa(d,p):
    m1=d[0]
    m2=d[1]
    m3=d[2]
    a11=p[0,0]
    a12=p[0,1]
    a13=p[0,2]
    a22=p[1,1]
    a23=p[1,2]
    a33=p[2,2]
    aux0=(4.*(a13**2))+((a22**2)+((4.*(a23**2))+((2.*(a22*a33))+(5.*(a33**2)))));
    aux1=(a22*((3.*((m1**2)+(m2**2)))+(m3**2)))+(a33*((3.*(m1**2))+((m2**2)+(3.*(m3**2)))));
    aux2=a11*((4.*((3.*(a12*(m1*m2)))+((3.*(a13*(m1*m3)))+(a23*(m2*m3)))))+aux1);
    aux3=(3.*((a11**2)*((5.*(m1**2))+((m2**2)+(m3**2)))))+((4.*((a12**2)*((3.*((m1**2)+(m2**2)))+(m3**2))))+(2.*aux2));
    aux4=(8.*(a12*((((3.*a22)+a33)*(m1*m2))+(2.*(((a23*m1)+(a13*m2))*m3)))))+aux3;
    aux5=(8.*(((a13*((a22+(3.*a33))*m1))+(3.*(a23*((a22+a33)*m2))))*m3))+((3.*(aux0*(m3**2)))+aux4);
    aux6=(12.*((a23**2)*(m2**2)))+((6.*(a22*(a33*(m2**2))))+((3.*((a33**2)*(m2**2)))+aux5));
    aux7=(16.*(a13*(a23*(m1*m2))))+((4.*((a13**2)*(m2**2)))+((15.*((a22**2)*(m2**2)))+aux6));
    aux8=(4.*((a23**2)*(m1**2)))+((2.*(a22*(a33*(m1**2))))+((3.*((a33**2)*(m1**2)))+aux7));
    output=0.00952381*((12.*((a13**2)*(m1**2)))+((3.*((a22**2)*(m1**2)))+aux8));
    return output

def full_average_conv_baa(d,p):
    m1=d[0]
    m2=d[1]
    m3=d[2]
    a11=p[0,0]
    a12=p[0,1]
    a13=p[0,2]
    a22=p[1,1]
    a23=p[1,2]
    a33=p[2,2]
    aux0=(a13*(((3.*a11)+(a22+(3.*a33)))*m1))+((2.*(a12*(a13*m2)))+(a23*((a11+(3.*(a22+a33)))*m2)));
    aux1=(6.*(a11*a22))+((9.*(a22**2))+((8.*(a23**2))+((4.*((a11+a22)*a33))+(3.*(a33**2)))));
    aux2=(-4.*(((2.*(a12*(a23*m1)))+aux0)*m3))+(((9.*(a11**2))+((12.*(a12**2))+((8.*(a13**2))+aux1)))*(m3**2));
    aux3=(6.*(a11*(a33*(m2**2))))+((4.*(a22*(a33*(m2**2))))+((9.*((a33**2)*(m2**2)))+aux2));
    aux4=(4.*(a11*(a22*(m2**2))))+((3.*((a22**2)*(m2**2)))+((8.*((a23**2)*(m2**2)))+aux3));
    aux5=(9.*((a11**2)*(m2**2)))+((8.*((a12**2)*(m2**2)))+((12.*((a13**2)*(m2**2)))+aux4));
    aux6=(-12.*(a12*(a22*(m1*m2))))+((-8.*(a13*(a23*(m1*m2))))+((-4.*(a12*(a33*(m1*m2))))+aux5));
    aux7=(6.*(a22*(a33*(m1**2))))+((9.*((a33**2)*(m1**2)))+((-12.*(a11*(a12*(m1*m2))))+aux6));
    aux8=(9.*((a22**2)*(m1**2)))+((12.*((a23**2)*(m1**2)))+((4.*(a11*(a33*(m1**2))))+aux7));
    aux9=(8.*((a12**2)*(m1**2)))+((8.*((a13**2)*(m1**2)))+((4.*(a11*(a22*(m1**2))))+aux8));
    output=0.00952381*((3.*((a11**2)*(m1**2)))+aux9);
    return output

def full_average_conv_aab(d,p):
    m1=d[0]
    m2=d[1]
    m3=d[2]
    a11=p[0,0]
    a12=p[0,1]
    a13=p[0,2]
    a22=p[1,1]
    a23=p[1,2]
    a33=p[2,2]
    aux0=(8.*(a13**2))+((2.*(a22**2))+((8.*(a23**2))+((-3.*(a22*a33))+(3.*(a33**2)))));
    aux1=((a12**2)*((8.*((m1**2)+(m2**2)))+(5.*(m3**2))))+((a11**2)*((3.*(m1**2))+(2.*((m2**2)+(m3**2)))));
    aux2=(2.*(((a13*(((-2.*a22)+a33)*m1))+(a23*((a22+a33)*m2)))*m3))+((aux0*(m3**2))+aux1);
    aux3=(-3.*(a22*(a33*(m2**2))))+((2.*((a33**2)*(m2**2)))+((6.*(a12*(((a23*m1)+(a13*m2))*m3)))+aux2));
    aux4=(5.*((a13**2)*(m2**2)))+((3.*((a22**2)*(m2**2)))+((8.*((a23**2)*(m2**2)))+aux3));
    aux5=(2.*((a33**2)*(m1**2)))+((6.*(a13*(a23*(m1*m2))))+((2.*(a12*((a22+(-2.*a33))*(m1*m2))))+aux4));
    aux6=(8.*((a13**2)*(m1**2)))+((2.*((a22**2)*(m1**2)))+((5.*((a23**2)*(m1**2)))+aux5));
    aux7=(4.*(a23*(m2*m3)))+((3.*(a33*(m3**2)))+(a22*((3.*((m1**2)+(m2**2)))+(m3**2))));
    aux8=(3.*(a33*(m1**2)))+((-2.*(a12*(m1*m2)))+((a33*(m2**2))+((-2.*(a13*(m1*m3)))+aux7)));
    output=0.00952381*((aux6-(a11*aux8))-(a22*(a33*(m1**2))));
    return output

def full_average_conv_abc(d,p):
    m1=d[0]
    m2=d[1]
    m3=d[2]
    a11=p[0,0]
    a12=p[0,1]
    a13=p[0,2]
    a22=p[1,1]
    a23=p[1,2]
    a33=p[2,2]
    aux0=(((5.*(a13**2))+((3.*(a22**2))+((5.*(a23**2))+(a33**2))))-(a22*a33))*(m3**2);
    aux1=((a12**2)*((5.*((m1**2)+(m2**2)))+(11.*(m3**2))))+((a11**2)*((m1**2)+(3.*((m2**2)+(m3**2)))));
    aux2=(-4.*(a12*(((a22+(-2.*a33))*(m1*m2))+(3.*(((a23*m1)+(a13*m2))*m3)))))+aux1;
    aux3=(4.*(a13*(((2.*a22)-a33)*(m1*m3))))+((-4.*(a23*((a22+a33)*(m2*m3))))+(aux0+aux2));
    aux4=((a22**2)*(m2**2))+((5.*((a23**2)*(m2**2)))+((3.*((a33**2)*(m2**2)))+aux3));
    aux5=(3.*((a33**2)*(m1**2)))+((-12.*(a13*(a23*(m1*m2))))+((11.*((a13**2)*(m2**2)))+aux4));
    aux6=(3.*((a22**2)*(m1**2)))+((11.*((a23**2)*(m1**2)))+((-5.*(a22*(a33*(m1**2))))+aux5));
    aux7=(a33*((m1**2)+((5.*(m2**2))+(m3**2))))+(a22*((m1**2)+((m2**2)+(5.*(m3**2)))));
    aux8=((5.*((a13**2)*(m1**2)))+aux6)-(a11*((4.*((a12*(m1*m2))+((a13*(m1*m3))+(-2.*(a23*(m2*m3))))))+aux7));
    output=0.00952381*(aux8-(a22*(a33*(m2**2))));
    return output


# build molecules from smiles using RDKit
def build_molecules(smiles):
    mols=[]
    errs=[]
    for i in range(0,len(smiles)):
        if Chem.MolFromSmiles(smiles[i]) is None: 
            errs.append(i)
            mols.append(None)
            continue
        else:
            mols.append(Chem.MolFromSmiles(smiles[i])) 
    return mols, errs


def generate_molecule_pictures(filenames,smiles):
    mols,errs =  build_molecules(smiles)   
    for m,mm in enumerate(mols):
        Draw.MolToFile(mm,"./{}.svg".format(filenames[m]),ImgSize=(100,100))
    for filename in os.listdir("./"):
        if filename.endswith(".svg"):
                svg_path=filename
                png_name=filename[:-3]+"png"
                svg2png(url=svg_path,write_to=png_name)
                print("Saved ",png_name)


#############################{main}#####################################


os.chdir(dirname)
print(dirname)

v0= math.pow(10, 7)/phys_param['laser']
scalingIR=math.pow(scaling_factors['D']/scaling_factors['A'],2)*scaling_factors['IRfac']
scaling,scalingexp= calc_scaling()
scaling=45*math.pow(scaling_factors['A'],4)*scaling

mol=0
start = timer()
Fr=[]
IR=[]
R_aa=[]
R_ab=[]
P_aaa=[]
P_baa=[]
P_aab=[]
P_abc=[]
filenames=[]

for filename in os.listdir("./"):
        if filename.endswith(file_ext):
                fr,Z,Q,D,P0,nat,aniso=load_file(filename)
                nmodes=len(fr)
                prod_aaa=np.zeros_like(fr)
                prod_baa=np.zeros_like(fr)
                prod_aab=np.zeros_like(fr)
                prod_abc=np.zeros_like(fr)
                r_aa=np.zeros_like(fr)
                r_ab=np.zeros_like(fr)
                ir_av=np.zeros_like(fr)

                for m in range(nmodes):

                    P=symmetrize(P0[m,:,:])
                  
                    scalingR=scaling* math.pow(v0 - fr[m], 4) / (
                                fr[m] * (1 - math.exp(scalingexp * fr[m]))) # Stokes for thermal population
                    
                    scalingTHOR=scaling* math.pow(v0 + fr[m], 4) / fr[m] # anti-Stokes without population
                    scalingaR=scaling* math.pow(v0 + fr[m], 4) / fr[m] *(
                           1/(-1+math.exp(-scalingexp * fr[m]))) # # anti-Stokes for thermal population
  
                    prod_aaa[m]=scalingIR *scalingTHOR*full_average_conv_aaa(D[m,:],P)
                    prod_baa[m]=scalingIR *scalingTHOR*full_average_conv_baa(D[m,:],P)
                    prod_aab[m]=scalingIR *scalingTHOR*full_average_conv_aab(D[m,:],P)
                    prod_abc[m]=scalingIR *scalingTHOR*full_average_conv_abc(D[m,:],P)
                    ir_av[m]=scalingIR*full_average_IR(D[m,:]) 
                    r_aa[m]=scalingR*full_average_R_pol(P,np.array([0,0,1]),np.array([0,0,1]))
                    r_ab[m]=scalingR*full_average_R_pol(P,np.array([0,0,1]),np.array([0,1,0]))
                Fr.append(fr)
                IR.append(ir_av)
                R_aa.append(r_aa)
                R_ab.append(r_ab)
                P_aaa.append(prod_aaa) 
                P_baa.append(prod_baa) 
                P_aab.append(prod_aab) 
                P_abc.append(prod_abc) 
                filenames.append(filename.split(".")[0])
                mol+=1

smiles=copy.deepcopy(filenames)
found=np.zeros((np.shape(smiles)))
with open(smilesfile) as inpfile:
    line=inpfile.readline()
    while(line):
        spl=line.split()
        for n,fn in enumerate(filenames):
            if fn in spl[1]: 
        #        smiles[n]=Chem.MolToSmiles(Chem.MolFromSmiles(spl[1])) # convert to RDkit SMILES format so duplicates can be identified in database
                smiles[n]=spl[0]
                found[n]=1
        line=inpfile.readline()
    
print("found ", int(np.sum(found)), " molecules in smiles file")

write_files_pol(Fr,IR, R_aa,R_ab,P_aaa,P_baa,P_aab,P_abc, smiles, filenames,found)

generate_molecule_pictures(filenames,smiles)

end = timer()
print("Time {:.2e} sec".format(end - start), " molecules: ",mol)
