#! /usr/bin/env python
#  Filename: rename_pdb_for_chimera.py

######################################################################
#                                                                    #
# Created: Maximilian N. Andrews, July 18, 2005                      #
#  Ciamician, Bologna University                                     #
#                                                                    #
# Function:                                                          #
#                                                                    #
#  This script adds a chain ID to all the subunits (column 22) in    #
#  the PDB file.                                                     #
#  The chains are named 'A','B','C',etc. until all of the subunits   #
#  have been labeled.                                                #
#                                                                    #
#  For the subunits to be renamed properly the PDB must contain a    #
#  line beginning with TER after each subunit, for example:          #
#  ...                                                               #
#  ATOM   1570  OXT HIE   103      24.590  -9.249  30.054            #
#  TER                                                               #
#  ATOM   1571  N   SER   104      16.936   7.767  30.236            #
#  ...                                                               #
#                                                                    #
# Notes:                                                             #
#                                                                    #
#  The original file is first read by the script, then copied to a   #
#  backup file - with the same name, but with "~" added at the end   #
#  of it - and then it is renamed and saved with the original name.  #
#                                                                    #
# Usage:                                                             #
#                                                                    #
#  rename_pdb.py arg1 arg2                                           #
#  arg1 is the the number of chains that have to be renamed;         #
#  arg2 is the name of the input file;                               #
#  example: rename_pdb.py 4 MY_PDB.pdb                               #
#                                                                    #
######################################################################


import sys
import os
import shutil

def BAILING_OUT(N_ARG):
     """ The srcipt needs two arguments:
     [int] The number of subunits to be renamed;
     [str] The name of the PDB file to be renamed.

     The correct use of this script is:

     rename_pdb.py arg1 arg2

     arg1 is the number of chains in the protein;
     arg2 is the name of the input file;

     example: rename_pdb.py 4 MY_PDB.pdb
     """
     if N_ARG != 3:
             print BAILING_OUT.__doc__
             sys.exit()
     else:
         return

BAILING_OUT(len(sys.argv))

NUMBER_CHAINS = sys.argv[1]
PDB_IN = sys.argv[2]

print 'Number of chains is: ',sys.argv[1]
print 'The PDB file to rename is: ',sys.argv[2]
PDB_OUT = os.path.splitext(PDB_IN)[0]+'-collect.pdb'
PDB_BAK = PDB_IN+'~'
f = open(PDB_IN,'r')
lines = f.readlines()
f.close()

print 'Making backup file '+PDB_BAK+'...'
shutil.copyfile(PDB_IN,PDB_BAK)
print '...done!'

print 'Renaming '+PDB_IN+'...'
CHAIN = 65 # this is the integer that corrisponds to the character "A"
MAX_CHAIN = int(NUMBER_CHAINS) + CHAIN
LETTER = chr(CHAIN)
NEW_PDB = []

for i in range(len(lines)):
     line = lines[i]
     if line.startswith('TER'):
         CHAIN = CHAIN + 1

     if CHAIN < MAX_CHAIN:
             LETTER = chr(CHAIN)
     else:
             LETTER = ' '

     if line.startswith('ATOM'):
          newline=line[0:21]+LETTER+line[22:]
          NEW_PDB.append(newline)

     elif line.startswith('TER') and len(line) > 21: # just in case  there's a chain identifier in the "TER" line... (the condition len (line) is added just in case the TER line is a just few characters  long... i.e. it would add a blank space at the start of the following  line)
         newline=line[0:21]+' '+line[22:]
         NEW_PDB.append(newline)
     else: # if line doesn't start with either ATOM or TER, it will  print it as is
         NEW_PDB.append(line)

print '...done!'

print 'Writing renamed file '+PDB_OUT+'...'
f = open(PDB_OUT,'w')
for i in range(len(NEW_PDB)):
     f.write(NEW_PDB[i])
f.close()
print '...done!'