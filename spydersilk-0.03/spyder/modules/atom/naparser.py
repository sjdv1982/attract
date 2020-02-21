# Copyright 2008, 2009, 2011 Marc van Dijk, Sjoerd de Vries
# This file is part of the Spyder module: "atom" 
# For licensing information, see LICENSE.txt 

from __future__ import print_function

import os,sys,re
from naparser_params import *


class fileoutputclass:
  
        def __init__(self):
                self.content = ""
        def write(self, text):
                self.content += text
                
class naparser:
        
        def __init__(self, pdbcontents, warnings=True):
                
                self.pdbcontents = pdbcontents
                                
                self.atcounter = 0
                self.allatom_line = []
                self.title = []
                self.header = []
                self.footer = []
                self.end = []
                self.model = []
                self.label = []
                self.atnum = []
                self.atname = []
                self.atalt = []
                self.resname = []
                self.chain = []
                self.resnum = []
                self.coord = []
                self.occ = []
                self.b = []
                self.segid= []
                
                self.warnings = warnings
                
                self.__readPDB()
                self.__NAresid1to3()
                self.__IUPACtoCNS()
                
                self.pdbcontents = self.writePDB()
                                
        
        def __NAresid1to3(self):
                
                """
                Convert list of 1-letter nucleic-acid code sequence to 3-letter code and update resname
                """
                
                seq3 = []
                
                for resid1 in self.resname:
                        try:
                                resid3 = NAres3[NAres1.index(resid1.upper())]         # If NAresid is one-letter code, convert to three-letter code
                                seq3.append(resid3)
                        except ValueError as err:
                                if resid1.upper() in AAres3:                                        # If resid is amino-acid three letter code, just append
                                        seq3.append(resid1.upper())                                        # Amino-acid one letter code in PDB not accepted(expected)
                                elif resid1.upper() == 'HOH ':                                        # Waters are neglected, just append.
                                        seq3.append(resid1.upper())
                                elif resid1.upper() in NAres3:                                         # If NAresid allready in three letter code, just append
                                        seq3.append(resid1.upper())
                                else:
                                        if self.warnings:
                                          print("WARNING: no match for residue: %s" % (resid1),file=sys.stderr)         # If not of the above, raise exception.
                                        seq3.append(resid1.upper())
                
                if len(seq3) == len(self.resname):
                        self.resname = seq3
                else:
                        pass
        
        def __IUPACtoCNS(self):
                
                """
                Convert IUPAC atom type notation to CNS atom type notation. Get info from IUPAC and CNS lists from params.py.
                Currently only conversion of nucleic-acid atom types.
                """
                
                newatomseq = []
                  
                for atom in self.atname:
                        try:
                                newatom = CNS[IUPAC.index(atom)]
                                newatomseq.append(newatom)
                        except ValueError:
                                newatomseq.append(atom)
                
                if len(newatomseq) == len(self.atname):
                        self.atname = newatomseq
                else:
                        pass
        
        def __readPDB(self):
                
                
                lines = self.pdbcontents.splitlines(True)
                self.__readPDBlines(lines)
        
        def __readPDBlines(self, lines):
                
                i = 0
                atom_hetatm = re.compile('(ATOM  |HETATM)')
                head = re.compile('^(HEADER|COMPND|SOURCE|JRNL|HELIX|REMARK|SEQRES|CRYST1|SCALE|ORIG)')
                foot = re.compile('(CONECT |MASTER)')
                end = re.compile('(END)')
                model = re.compile('(MODEL)')
                for line in lines:
                    if atom_hetatm.match(line):
                      line = line[:-1]
                      self.allatom_line.append(line)
                      self.label.append(line[0:6])                         #atom label
                      self.atnum.append(int(line[6:12]))                 #atom number
                      self.atname.append(line[12:16])                        #atom type
                      self.atalt.append(line[16:17])
                      self.resname.append(line[17:21])                 #residu name
                      self.chain.append(line[21])                         #chain
                      self.resnum.append(int(line[22:26]))        #residu number

                      try:
                              self.coord.append((float(line[30:38]), float(line[38:46]), float(line[46:54]))) #X,Y,Z coordinates
                      except:
                              if not line[0:3] == 'TER' or line[0:5] == 'MODEL':
                                      print("ERROR: coordinate error in line:",file=sys.stderr)
                                      print("     ", line,file=sys.stderr)

                      try:
                              self.occ.append(float(line[54:60]))
                      except:
                              self.occ = ([  1.00]*len(self.atnum))

                      try:
                              self.b.append(float(line[60:66]))                         #B factor
                      except:
                              self.b = ([  0.00]*len(self.atnum))

                      try:
                              self.segid.append(line[72])                        #SEGID
                      except:
                              self.segid.append(blank)

                      self.atcounter += 1
                      i = i + 1

                    elif head.match(line):
                      self.header.append(line[:-1])
                    elif foot.match(line):
                      self.footer.append(line[:-1])
                    elif end.match(line):
                      self.end.append(line[:-1])
                    elif model.match(line[:-1]):
                      self.model.append(line)
          
        def writePDB(self,noheader=False,nohetatm=False):
                
                """
                Saves the Protein class object to a PDB-format file
                if noheader = True, no header (REMARK etc.) or footer lines are written
                if nohetatm = True, no hetero atoms are written
                """
                out = fileoutputclass()
                
                if noheader == False:
                        for i in range(len(self.title)):
                                out.write('%s\n' % self.title[i])
                        for i in range(len(self.header)):
                                out.write("%s\n" % self.header[i])
                
                if nohetatm == False:
                        for i in xrange(len(self.resnum)):
                                if self.label[i] == 'ATOM  ':
                                        self.__writePDBline(out,i)
                                elif self.label[i] == 'TER   ':
                                        out.write('TER   \n')
                                elif self.label[i] == 'HETATM':
                                        out.write("%s\n" % self.allatom_line[i])
                else:
                        for i in xrange(len(self.resnum)):
                                if self.label[i] == 'ATOM  ':
                                          self.__writePDBline(out,i)
                                elif self.label[i] == 'TER   ':
                                        out.write('TER   \n')
                
                if nofooter == False:
                        for i in range(len(self.footer)):
                                        out.write("%s\n" % self.footer[i])
                
                if len(self.end) == 3:
                        for i in range(len(self.end)):
                                out.write("%s\n" % self.end[i])
                else:
                        self.end = ['END']
                        for i in range(len(self.end)):
                                out.write("%s\n" % self.end[i])
                
                return out.content
        
        def __writePDBline(self,FD,i):
                
                FD.write('%-6s%5i %-4s%1s%-4s%1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f%12s\n' % 
                        (self.label[i],self.atnum[i],self.atname[i],self.atalt[i],self.resname[i],self.chain[i], 
                         self.resnum[i],self.coord[i][0],self.coord[i][1],self.coord[i][2],self.occ[i],self.b[i],self.segid[i]))


