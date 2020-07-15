#!/usr/bin/python

from water import *
import Bio.PDB as bpdb
#import numpy
import subprocess,os
import sys,getopt
from shutil import copyfile

"""Main program for optimizing water orientations in neutron structures"""


try:
   opts, args = getopt.getopt(sys.argv[1:],"hp:z:o:rnet:c:",["help","pdb=","mtz=","out=","refine","norefine","omit","threshold=","cycles="])
except getopt.GetoptError:
   print 'nwo.py -p <pdbfile> -o <outputpdbfile> -z <reflectionfile> -t <rscc-threshold> -c <max-number-of-cycles> -r (refine at each step) OR -n (no refinement) OR -e (use omit map)'
   sys.exit(2)
for opt, arg in opts:
   if opt == '-h':
      print 'nwo.py -p <pdbfile> -o <outputpdbfile> -z <reflectionfile> -t <rscc-threshold> -c <max-number-of-cycles> -r (refine at each step) OR -n (no refinement) OR -e (use omit map)'
      sys.exit()
   elif opt in ("-p", "--pdb"):
      pdb = arg
   elif opt in ("-o", "--out"):
      outfile = arg
   elif opt in ("-z", "--mtz"):
      mtz = arg
   elif opt in ("-t", "--threshold"):
      threshold = float(arg)
   elif opt in ("-c", "--cycles"):
      max_no_of_cycles = int(arg)
   elif opt in ("-r", "--threshold"):
      flag = "refine"
   elif opt in ("-n", "--norefine"):
      flag = "norefine"
   elif opt in ("-e", "--omit"):
      flag = "omit"


############# FIX  to copy alternate conformations
############# From github biopython issue 455

copy=bpdb.Atom.copy
def myCopy(self):
    shallow = copy.copy(self)
    for child in self.child_dict.values():
        shallow.disordered_add(child.copy())
    return shallow
bpdb.Atom.DisorderedAtom.copy=myCopy

############ START


#Read CRYST1 line from pdb

with open(pdb, "r") as input_pdb:
     for line in input_pdb:
         if(line[0:6]=="CRYST1"):
             cryst1 = line

with open("cryst1", "w") as cryst1_file:
    cryst1_file.write(cryst1)


ini_structure = bpdb.PDBParser().get_structure('ini', pdb) # load the initial file (from ReadySet/Maestro/random guess) 

#read waters from PDB


best_structure = ini_structure.copy()
best_residue_list = bpdb.Selection.unfold_entities(best_structure[0], 'R')
best_waters_list = [w for w in  best_residue_list if w.get_full_id()[3][0] == "W"]


s = best_structure.copy()
residue_list = bpdb.Selection.unfold_entities(s[0], 'R')
waters_list = [w for w in  residue_list if w.get_full_id()[3][0] == "W"]

#run a first refinement / make maps without waters

class NoWaters(bpdb.Select):
    def accept_residue(self,residue):
        if (residue.get_resname()=="HOH"):
            return False
        else:
            return True


copyfile(pdb,"x.pdb")
copyfile(mtz,"x.mtz")
io = bpdb.PDBIO()

if(flag=="omit"):
    io.set_structure(ini_structure)
    io.save("no-waters.pdb",NoWaters())
    return_code = subprocess.call("phenix.maps no-waters.pdb x.mtz scattering_table=neutron labels=F-obs,SIGF-obs",shell=True)
    copyfile("no-waters_map_coeffs.mtz","x.mtz")

if(flag=="refine"):
    return_code = subprocess.call(" phenix.refine refine.def --overwrite --unused_ok > log",shell=True)
    refined_structure = bpdb.PDBParser().get_structure('r', 'x-refined_001.pdb')



#run real space correlation

if(flag=="refine"):
    return_code = subprocess.call("phenix.real_space_correlation x-refined_001.pdb x-refined_001.mtz use_hydrogens=True scattering_table=neutron data_labels=F-obs,SIGF-obs > rscc",shell=True)
elif(flag=="norefine"):
    return_code = subprocess.call("phenix.real_space_correlation x.pdb x.mtz use_hydrogens=True scattering_table=neutron data_labels=F-obs,SIGF-obs > rscc",shell=True)
elif(flag=="omit"):
    return_code = subprocess.call("phenix.real_space_correlation x.pdb x.mtz use_hydrogens=True scattering_table=neutron map_coefficients_label=2FOFCWT,PH2FOFCWT > rscc",shell=True)

#read rscc output
start = "<----id string---->  occ     ADP      CC   Rho1   Rho2"
found = False

with open('rscc') as rscc_file:
       for line in rscc_file:
           if found:
               chain = line[1]
               resid = int(line[9:13])
               atomn = line[17:19].strip()
               rscc = float(line[35:41])
               resname = line[5:8]
               if (resname == "HOH"):
                   s[0][chain][('W',resid,' ')][atomn].bfactor = rscc
           if line.strip() == start:
               found = True

# Make a first estimation of waters that are in wrong orientations

bad_waters = 0

for w in best_waters_list:
    w['D1'].set_bfactor(waters_list[best_waters_list.index(w)]['D1'].bfactor)
    w['D2'].set_bfactor(waters_list[best_waters_list.index(w)]['D2'].bfactor)


for w in waters_list:
    if(w['O'].bfactor < 0.5):
        print("The oxygen atom in water " , w.id[1], "does not seem to be visible in the neutron map. Consider deleting it.")
    if((w['D1'].bfactor < threshold) or (w['D2'].bfactor  < threshold)):
        bad_waters+=1
    if((w['D1'].bfactor + w['D2'].bfactor) > (best_waters_list[waters_list.index(w)]['D1'].bfactor + best_waters_list[waters_list.index(w)]['D2'].bfactor)):
        best_waters_list[waters_list.index(w)]['D1'].set_coord(w['D1'].coord)
        best_waters_list[waters_list.index(w)]['D2'].set_coord(w['D2'].coord)
        best_waters_list[waters_list.index(w)]['D1'].set_bfactor(w['D1'].bfactor)
        best_waters_list[waters_list.index(w)]['D2'].set_bfactor(w['D2'].bfactor)

        
#Start generating new water orientations

no_of_cycles = 0

while(no_of_cycles < max_no_of_cycles and bad_waters > 0):
    for w in waters_list:
        generate_new_water_orientation(w, threshold)
        for a in w:
            if(flag=="refine"):
                a.set_bfactor(refined_structure[a.full_id[1]][a.full_id[2]][a.full_id[3]][a.full_id[4][0]].bfactor)
            else:
                a.set_bfactor(ini_structure[a.full_id[1]][a.full_id[2]][a.full_id[3]][a.full_id[4][0]].bfactor)
    io = bpdb.PDBIO()
    io.set_structure(s)
    io.save("temp")
    filenames = ['cryst1', 'temp']
    #add CRYST1 header to pdb
    with open('x.pdb', 'w') as pdb:   
        for name in filenames:
            with open(name) as inp:
                for line in inp:
                    pdb.write(line)
    if(flag=="refine"):
        return_code = subprocess.call(" phenix.refine refine.def --overwrite  --unused_ok > log",shell=True)
        return_code = subprocess.call(" phenix.real_space_correlation x-refined_001.pdb x-refined_001.mtz detail=atom use_hydrogens=True scattering_table=neutron data_labels=F-obs,SIGF-obs > rscc",shell=True)
    elif(flag=="norefine"):
        return_code = subprocess.call(" phenix.real_space_correlation x.pdb x.mtz use_hydrogens=True detail=atom scattering_table=neutron data_labels=F-obs,SIGF-obs > rscc",shell=True)
    elif(flag=="omit"):
        return_code = subprocess.call("phenix.real_space_correlation x.pdb x.mtz use_hydrogens=True detail=atom scattering_table=neutron map_coefficients_label=2FOFCWT,PH2FOFCWT > rscc",shell=True)
    found = False
    with open('rscc') as rscc_file:
        for line in rscc_file:
            if found:
                chain = line[1]
                resid = int(line[9:13])
                atomn = line[17:19].strip()
                rscc = float(line[35:41])
                resname = line[5:8]
                if (resname == "HOH"):
                   s[0][chain][('W',resid,' ')][atomn].bfactor = rscc
            if line.strip() == start:
                found = True
    bad_waters = 0
    for w in waters_list:
        if((w['D1'].bfactor < threshold) or (w['D2'].bfactor  < threshold)):
            bad_waters+=1
        if((w['D1'].bfactor + w['D2'].bfactor) > (best_waters_list[waters_list.index(w)]['D1'].bfactor + best_waters_list[waters_list.index(w)]['D2'].bfactor)):
            best_waters_list[waters_list.index(w)]['D1'].set_coord(w['D1'].coord)
            best_waters_list[waters_list.index(w)]['D2'].set_coord(w['D2'].coord)
            best_waters_list[waters_list.index(w)]['D1'].set_bfactor(w['D1'].bfactor)
            best_waters_list[waters_list.index(w)]['D2'].set_bfactor(w['D2'].bfactor)


    no_of_cycles+=1
    if(flag=="refine"):
        refined_structure = bpdb.PDBParser().get_structure('r', 'x-refined_001.pdb')

#Delete deuterium atoms that are still not correct and write out final pdb file


class GoodRSCC(bpdb.Select):
    def accept_atom(self,atom):
        if (atom.occupancy == 0.0):
            return False
        else:
            return True

for w in best_waters_list:
    for a in w:
        if(a.id != "O" and a.bfactor < (threshold - 0.2)):
           a.occupancy = 0.0
        if(flag=="refine"):
           a.bfactor = refined_structure[a.full_id[1]][a.full_id[2]][a.full_id[3]][a.full_id[4][0]].bfactor
        else:
           a.bfactor = ini_structure[a.full_id[1]][a.full_id[2]][a.full_id[3]][a.full_id[4][0]].bfactor

io.set_structure(best_structure)
io.save(outfile, GoodRSCC())

#Remove unneeded files

os.remove("x.pdb")
os.remove("x.mtz")
os.remove("cryst1")
if(flag=="omit"):
   os.remove("no-waters.pdb")
   os.remove("no-waters_map_coeffs.mtz")
   os.remove("no-waters_2mFo-DFc_map.ccp4")
#'''       



        
    




       


   




