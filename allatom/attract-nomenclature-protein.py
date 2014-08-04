#!/usr/bin/python2.7
###########################
###  L I P R A R I E S  ###
###########################
import sys
import os
#################################################
pdbout=open('attract-'+sys.argv[1],'w')
L=open(sys.argv[1],'r').readlines()
HIP=[]
HIE=[]

for line in L:
	if len(line.split())<4:continue
	if line[17]=='H':
		if line.split()[2]=='N':hd1=0
		if line.split()[2]=='HD1': hd1=1
		if line.split()[2]=='HE2':
			if hd1==1:HIP.append(line.split()[4])
			else:HIE.append(line.split()[4])

for line in L:
	if len(line.split())<4:
		pdbout.write(line)
		continue
	if line[17]=='H':
		hid=True
		for res in HIP:
			if line.split()[4]==res:
				pdbout.write(line[:17]+'HIP'+line[20:])
				hid=False
		for res in HIE:
			if line.split()[4]==res:
				pdbout.write(line[:17]+'HIE'+line[20:])
				hid=False
		if hid:	pdbout.write(line[:17]+'HID'+line[20:])
	elif line[13:17]=='HG1 ' and (line[17]=='C' or line[17]=='S'): pdbout.write(line[:13]+'HG '+line[16:])
	elif line[13:18]=='CD  I': pdbout.write(line[:13]+'CD1'+line[16:])
	elif line[13:16]=='OT1': pdbout.write(line[:13]+'O  '+line[16:])
	elif line[13:16]=='OT2': pdbout.write(line[:13]+'OXT'+line[16:])
	else:	pdbout.write(line)

pdbout.close()
