#!/usr/bin/python

################################################################################
# prepare_pdb.py                                                               #
# Version: 040909                                                              #
# Goal:    the script calls the program 'pdb2gmx' from the GROMACS suite to    #
#          produce a topology and a 'gro' file and merges the two files into   #
#          into one "pdb-like" file suitable for easymifs                      #
# Usage:   prepare_pdb.py PDB [-k]                                             #
#          The PDB file cannot contain missing atoms or residues or HETATM     #
#          records. The -k option will produce a clean PDB file (by stripping  #
#          all the HETATM records)                                             #
################################################################################

import sys
import os
import re

################################################################################
# CONSTANTS                                                                    #
################################################################################

ff = "G43b1" # in vacuo force field

################################################################################
# MAIN PROGRAM                                                                 #
################################################################################

## parse the parameters
if len(sys.argv) < 2:
  print("Usage: process_gromacs_input.py PDB [-k]")
  sys.exit(1)
pdbName = sys.argv[1]
stripHET = False
if len(sys.argv) == 3:
  if sys.argv[2] == "-k":
    stripHET = True

## read the PDB file
pdbStem = pdbName.split(".pdb")
if len(pdbStem) == 1:
  pdbStem = pdbName.split(".ent")
  if len(pdbStem) == 1:
    print("The file should have extension .pdb")
    sys.exit(1)
pdbStem = pdbStem[0]
pdb = open(pdbName, "r")
if stripHET: # strip the HETATMS
  strippedPdb = open(pdbStem + "_stripped.pdb", "w")
for line in pdb:
  if line[:6] != "HETATM" and stripHET: 
    strippedPdb.write(line)
  elif line[:6] == "HETATM" and not stripHET:
    print("HETATM record found in the PDB file...\nPlease consider using the -k option or manually remove the HETATM records\n")
    pdb.close()
    sys.exit(2)
  
if stripHET:
  strippedPdb.close()
  pdbName = pdbStem + "_stripped.pdb" # make pdb2gmx read this file
pdb.close()

## call the pdb2gmx program
pathname = os.path.dirname(sys.argv[0])
gmxlib = pathname + os.sep + "pdb2gmx"
pdb2gmx = gmxlib + os.sep + "pdb2gmx"
os.environ['GMXLIB'] = gmxlib
commandLine = "\"" + pdb2gmx + "\"" + " -ignh -merge -f " + pdbName + \
              " -ff " + ff + " -o " + pdbStem + "_temp.pdb -p " + pdbStem + \
              ".top"
os.system(commandLine)

# read both files storing the necessary information
grofileName = pdbStem + "_temp.pdb"
grofile = open(grofileName, "r")
data = grofile.readlines()
grofile.close()
grodata = []
for line in data: # skip the first two lines of the file and the last one
  if line[:4] == "ATOM":
    grodata.append(line)

topfileName = pdbStem + ".top"  
topfile = open(topfileName, "r")
data = topfile.readlines()
topfile.close()
topdata = []
flag = 0 # flag to signal the beginning of the atom section
for line in data:
  if flag and line == "\n": # end of the atom section
    break
  if flag: # atom section
    fields = line.split()
    atom_type, charge = fields[1], fields[6]
    topdata.append([atom_type, charge])
  if line.find("atom") != -1 and line.find("type") != -1 and \
     line.find("charge") != -1:
    flag = True
    continue

## merge the two files into a new one
numOfAtoms = len(grodata)
if numOfAtoms != len(topdata):
  print("GROFILE and TOPFILE contain a different number of atoms")
  sys.exit(1)
outfileName = pdbStem + ".easymifs"
outfile = open(outfileName, "w")
for i in range(numOfAtoms):
  # prepare the line according to the PDB specification
  line = grodata[i][:12] + topdata[i][0] + " " * (4 - len(topdata[i][0])) + \
         grodata[i][16:55] + topdata[i][1]
  outfile.write(line + "\n") 
outfile.close()

## remove temporary files
os.remove(pdbStem + "_temp.pdb")
os.remove("posre.itp")
os.remove(pdbStem + ".top")
if os.path.isfile(pdbStem + "_stripped.pdb"):
  os.remove(pdbStem + "_stripped.pdb")
