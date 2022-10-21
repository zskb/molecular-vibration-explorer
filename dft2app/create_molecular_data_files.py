# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 09:15:34 2022

@author: Zsuzsanna Koczor-Benda
"""

##########################{libraries}###################################

import numpy as np
import openbabel as ob # works with Openbabel version 2.4.1
from load_data import load_file 
import os


#######################{user editable variables}########################


# reading from formatted checkpoint file (.fchk) is more accurate than from output file (.out)
# to generate .fchk file from .chk use the formchk command of Gaussian (see https://gaussian.com/formchk/)
# NOTE: .out files are ALWAYS needed (for rotating molecule in reference orientation), .fchk files are optional but recommended
read_from_fchk=False

file_ext=".fchk"
if read_from_fchk==False:
    file_ext=".out" 
    
# directory where .fchk and .out files are located
dirname="C:\\\\your-local-path\\\\test-molecules"
# filenames are presumed to have format "freq-0003921.out" or "freq-0003921.fchk"
# where the number is an identifier of the molecule (can contain numbers/letters/symbols valid for filenames, except '.')

# we need a file called "smiles_data.txt" that contains SMILES code for each filename, see example
# use online databases/tools or molecule editors to get SMILES, generating SMILES directly from output is not reliable
smilesfile="test-molecules\\\\smiles_data.txt" 

# define a tag for the dft .out and .fchk filenames (ex. "freq-0003921.out" or "freq-0003921.fchk")
tagname="freq-"

# decide whether or not to orient the molecule with respect to a specific link
apply_rotation=False

# decide on the link to consider for orientation based on atomic number
# in our database we consider single S (16) atoms and links with Au (79) and H (1) atoms
n_ref=16
n_atnum=1


###########################{functions}##################################


def convert_to_mol(filename,molfilename):
    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "mol") 
    mol = ob.OBMol()
    obConversion.ReadFile(mol, filename)
    outMDL = obConversion.WriteString(mol)
    obConversion.WriteFile(mol,molfilename)  
    
def rotate_molecule(mol,filename,Q,D,P,phi=0,theta=0):  
    # rotate to reference orientation first, then do custom rotation with phi, theta
    
    # for rotating to reference position
    # locate thiol S atom and its first- and second-level neighbours
    S_coords=np.zeros(3)
    neighb2_idxs=np.zeros(4,dtype=int)
    neighb2_num=np.zeros(4,dtype=int)
    # find S atom with Au/H neighbour (in case there are more than one S atoms in molecule)
    # note there is only one thiol per molecule in the database
    thiol=False
    n1_indx=-1
    for a in ob.OBMolAtomIter(mol): 
        if a.GetAtomicNum()==n_ref: 
            thiol=False
            n=0
            for neighbour in ob.OBAtomAtomIter(a):
                if neighbour.GetAtomicNum()==n_atnum: 
                      thiol=True
            if thiol==True:
                S_indx=a.GetIdx()
                S_coords[0]=a.GetX()
                S_coords[1]=a.GetY()
                S_coords[2]=a.GetZ()
                for neighbour in ob.OBAtomAtomIter(a):   
                    if neighbour.GetAtomicNum()!=n_atnum:
                        n1_indx=neighbour.GetIdx()
                        for neighbour2 in ob.OBAtomAtomIter(neighbour):
                            if neighbour2.GetIdx()==S_indx:
                                n+=1
                                continue
                            neighb2_idxs[n]=neighbour2.GetIdx()                     
                            neighb2_num[n]=neighbour2.GetAtomicNum()
                            n+=1
                        break     
    if n1_indx==-1:
        print("Error: thiol group not found in molecule")
        print([[atom.GetX(),atom.GetY(),atom.GetZ()] for atom in ob.OBMolAtomIter(mol)] )
        return 0
    
    n2_indx=neighb2_idxs[np.argmax(neighb2_num)]
    
    coords_orig = [[atom.GetX(),atom.GetY(),atom.GetZ()] for atom in ob.OBMolAtomIter(mol)] 
    
    # translate S to origin
    translate=np.array([-S_coords[0],-S_coords[1],-S_coords[2]] )
    atnums= [atom.GetAtomicNum() for atom in ob.OBMolAtomIter(mol)] 
    coords=np.zeros(np.shape(coords_orig))
    for i in range(0,len(coords_orig)):
        coords[i]=coords_orig[i]+translate
    
    # rotate n1 to z axis, n2 to x axis
    n1_coords=coords[n1_indx-1]
    n1_norm=np.linalg.norm(n1_coords)
    n1_normed=n1_coords/n1_norm
    cosphi=n1_normed[2]/np.sqrt(1-np.power(n1_normed[1],2))
    sinphi=n1_normed[0]/np.sqrt(1-np.power(n1_normed[1],2))

    rot_matrix=np.zeros((3,3))
    rot_matrix[0,0]=cosphi
    rot_matrix[0,2]=-sinphi
    rot_matrix[2,0]=sinphi
    rot_matrix[1,1]=1
    rot_matrix[2,2]=cosphi

    rotcoords=np.transpose(np.matmul(rot_matrix,np.transpose(coords)))
    
    n1_coords2=rotcoords[n1_indx-1]
    n1_norm2=np.linalg.norm(n1_coords2)
    n1_normed2=n1_coords2/n1_norm2
    cosrho=n1_normed2[2]/np.sqrt(1-np.power(n1_normed2[0],2))
    sinrho=n1_normed2[1]/np.sqrt(1-np.power(n1_normed2[0],2))

    rot_matrix2=np.zeros((3,3))
    rot_matrix2[0,0]=1
    rot_matrix2[1,1]=cosrho
    rot_matrix2[1,2]=-sinrho
    rot_matrix2[2,1]=sinrho
    rot_matrix2[2,2]=cosrho
    
    rotcoords2=np.transpose(np.matmul(rot_matrix2,np.transpose(rotcoords)))
    
    n2_coords3=rotcoords2[n2_indx-1]
    n2_norm3=np.linalg.norm(n2_coords3)
    n2_normed3=n2_coords3/n2_norm3
    costheta=n2_normed3[0]/np.sqrt(1-np.power(n2_normed3[2],2))
    sintheta=n2_normed3[1]/np.sqrt(1-np.power(n2_normed3[2],2))

    rot_matrix3=np.zeros((3,3))
    rot_matrix3[0,0]=costheta
    rot_matrix3[0,1]=sintheta
    rot_matrix3[1,0]=-sintheta
    rot_matrix3[1,1]=costheta
    rot_matrix3[2,2]=1
    
    # rotate original coordinates, dipole and polarizability derivatives to reference 
    rot=np.matmul(rot_matrix3,np.matmul(rot_matrix2,rot_matrix)) 
    Drot=np.zeros_like(D)
    Prot=np.zeros_like(P)
    Qrot=np.zeros_like(Q)
    for m in range(np.shape(D)[0]):
        Drot[m]=np.matmul(np.transpose(rot),D[m]) 
        Prot[m]=np.matmul(np.transpose(rot),np.matmul(P[m],rot)) 
        Qrot[m]=np.transpose(np.matmul(rot,np.transpose(Q[m])))
    refcoords=np.transpose(np.matmul(rot,np.transpose(coords)))
        
#   save in .xyz and .mol formats for displaying later
    code=(filename[:-4]).split(tagname,1)[1]
    filename="{}.xyz".format(code)
    xyz="{} \n\n".format(len(atnums))
    molfile = open(filename, "w+")   
    molfile.write(xyz)
    for m in range(0,len(atnums)):
        molfile.write("  {}    {:.8f} {:.8f} {:.8f} \n".format(atnums[m],refcoords[m,0],refcoords[m,1],refcoords[m,2]))
    molfile.close()
    molfilename=filename[:-3]+"mol"
    convert_to_mol(filename,molfilename)
    print("Saved ",molfilename)
    return Drot,Prot,Qrot

def write_data_file(molcode,fr,Q,D,P,nat,aniso):
    filename=molcode+".dat"
    with open(filename,'w+') as ofile:
        ofile.write("Number of atoms \n{}\n--\n".format(nat))
        ofile.write("Anisotropy \n{:.3f}\n--\n".format(aniso))
        ofile.write("Frequencies \n")
        for f in fr:
            ofile.write("{:.12f}\n".format(f))
        ofile.write("--\n")
        ofile.write("Displacements \n")
        for n,f in enumerate(Q):
            ofile.write("Mode {}\n".format(n))
            for at in f:
                ofile.write("{:.12e} {:.12e} {:.12e}\n".format(at[0],at[1],at[2]))
        ofile.write("--\n")
        ofile.write("Dipole derivatives \n")
        for n,d in enumerate(D):
            ofile.write("Mode {}\n".format(n))
            ofile.write("{:.12e} {:.12e} {:.12e}\n".format(d[0],d[1],d[2]))
        ofile.write("--\n")
        ofile.write("Polarizability derivatives \n")
        for n,f in enumerate(P):
            ofile.write("Mode {}\n".format(n))
            for at in f:
                ofile.write("{:.12e} {:.12e} {:.12e}\n".format(at[0],at[1],at[2]))
        ofile.write("--\n")
        print("Written ",filename)
        print("-----------------------")
        
def manipulate_molfile(smilesfile):
    # write smiles code in .mol file 
    found=0
    with open(smilesfile) as inpfile:
        line=inpfile.readline()
        while(line):
            spl=line.split()
            for filename in os.listdir("./"):
                if ".mol" in filename:
                    fn=filename.split(".")[0]

                    if fn in spl[1]: 
                        smiles=spl[0]
                        with open(filename, 'r') as file:
                            # read a list of lines into data
                            data = file.readlines()
                        
                        # now change the 2nd line
                        data[2] = smiles+'\n'
                        
                        # and write everything back
                        with open(filename, 'w') as file:
                            file.writelines( data )
                        
                        found+=1
                        break
            line=inpfile.readline()       


#############################{main}#####################################
        

os.chdir(dirname)
print(dirname)
nm=0
for filename in os.listdir("./"):
    if filename.endswith(file_ext): 
        fname=filename.split(".")[0]
        molcode=fname.split(tagname,1)[1]
        try:
            fr,Z,Q,D,P,nat,aniso=load_file(filename)
            outfilename=dirname+"\\\\"+fname+ ".out"  # we need the .out file for OpenBabel
            obConversion = ob.OBConversion()
            obConversion.SetInAndOutFormats("out", "can") 
            obConversion.SetOptions(str("h"), ob.OBConversion.OUTOPTIONS)
            mol0 = ob.OBMol()
            obConversion.ReadFile(mol0, outfilename)
            print("Generated SMILES: ",(obConversion.WriteString(mol0)).split()[0])
            if apply_rotation==True:
                Drot,Prot,Qrot=rotate_molecule(mol0,outfilename,Q,D,P,phi=0,theta=0)
                write_data_file(molcode,fr,Qrot,Drot,Prot,nat,aniso)
            else:
                coords = [[atom.GetX(),atom.GetY(),atom.GetZ()] for atom in ob.OBMolAtomIter(mol0)]
                atnums= [atom.GetAtomicNum() for atom in ob.OBMolAtomIter(mol0)]
                #   save in .xyz and .mol formats for displaying later
                code=(outfilename[:-4]).split(tagname,1)[1]
                filename="{}.xyz".format(code)
                xyz="{} \n\n".format(len(atnums))
                molfile = open(filename, "w+")   
                molfile.write(xyz)
                for m in range(0,len(atnums)):
                    molfile.write("  {}    {:.8f} {:.8f} {:.8f} \n".format(atnums[m],coords[m][0],coords[m][1],coords[m][2]))
                    #molfile.write("  {}     \n".format(atnums[m]))
                molfile.close()
                molfilename=filename[:-3]+"mol"
                convert_to_mol(filename,molfilename)
                print("Saved ",molfilename)
                write_data_file(molcode,fr,Q,D,P,nat,aniso)
            nm+=1
        except:
            print("Error with molecule")

manipulate_molfile(smilesfile)
print("Saved ",nm," molecules")



