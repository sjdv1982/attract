import sys
def read_topology(atom, resn):
  bonds = []
  angles = []
  bbbonds = [('N','HN'),('N','CA'),('CA','C'),('C','O')]
  bbangles = [('HN','N','CA'),('N','CA','C'),('CA','C','O')]
  if resn == 'CYS':
    allbonds = bbbonds + [('CA','HA'),('CA','CB'),('CB','HB1'),('CB','HB2'),('CB','SG'),('SG','HG')]
    allangles = bbangles + [('N','CA','HA'),('N','CA','CB'),('CB','CA','HA'),('C','CA','HA'),('CA','CB','HB1'),('CA','CB','HB2'),('HB1','CB','HB2'),
			    ('HB1','CB','SG'),('HB2','CB','SG'),('CA','CB','SG'),('CB','SG','HG'),('CB','CA','C')]
    bonds = [x for x in allbonds if atom in x]
    angles = [x for x in allangles if atom in x]
    
  elif resn == 'ALA':
    allbonds = bbbonds + [('CA','HA'),('CA','CB'),('CB','HB1'),('CB','HB2'),('CB','HB3')]
    allangles = bbangles + [('N','CA','HA'),('N','CA','CB'),('CB','CA','HA'),('C','CA','HA'),('CB','CA','C'),('CA','CB','HB1'),
			    ('CA','CB','HB2'),('CA','CB','HB3'),('HB1','CB','HB2'),('HB1','CB','HB3'),('HB2','CB','HB3')]
    bonds = [x for x in allbonds if atom in x]
    angles = [x for x in allangles if atom in x]
    
  elif resn == 'ARG':
    allbonds = bbbonds + [('CA','HA'),('CA','CB'),('CB','HB1'),('CB','HB2'),('CB','CG'),('CG','HG1'),
			  ('CG','HG2'),('CG','CD'),('CD','HD1'),('CD','HD2'),('CD','NE'),
			  ('NE','HE'),('NE','CZ'),('CZ','NH1'),('NH1','HH11'),('NH1','HH12'),
			  ('CZ','NH2'),('NH2','HH21'),('NH2','HH22')]
    allangles = bbangles + [('N','CA','HA'),('N','CA','CB'),('CB','CA','HA'),('C','CA','HA'),('CB','CA','C'),('CA','CB','HB1'),
			    ('CA','CB','HB2'),('HB1','CB','HB2'),('CA','CB','CG'),('HB1','CB','CG'),('HB2','CB','CG'),('CB','CG','HG1'),
			    ('CB','CG','HG2'),('HG1','CG','HG2'),('HG1','CG','CD'),('HG2','CG','CD'),('CG','CD','HD1'),
			    ('CG','CD','HD2'),('HD1','CD','HD2'),('CB','CG','CD'),('HD1','CD','NE'),('HD2','CD','NE'),
			    ('CG','CD','NE'),('CD','NE','HE'),('CD','NE','CZ'),('NE','CZ','NH1'),('NE','CZ','NH2'),
			    ('HE','NE','CZ'),('CZ','NH1','HH11'),('CZ','NH1','HH12'),('HH11','NH1','HH12'),
			    ('CZ','NH2','HH21'),('CZ','NH2','HH22'),('HH21','NH2','HH22'), ('NH1','CZ','NH2')]
    bonds = [x for x in allbonds if atom in x]
    angles = [x for x in allangles if atom in x]    
    
  elif resn == 'ASN':
    allbonds = bbbonds + [('CA','HA'),('CA','CB'),('CB','HB1'),('CB','HB2'),('CB','CG'),
			  ('CG','OD1'),('CG','ND2'),('ND2','HD21'),('ND2','HD22')]
    allangles = bbangles + [('N','CA','HA'),('N','CA','CB'),('CB','CA','HA'),('C','CA','HA'),('CB','CA','C'),('CA','CB','HB1'),
			    ('CA','CB','HB2'),('HB1','CB','HB2'),('CA','CB','CG'),('HB1','CB','CG'),('HB2','CB','CG'),('CB','CG','OD1'),
			    ('CB','CG','ND2'),('CG','ND2','HD21'),('CG','ND2','HD22'),('HD21','ND2','HD22'),('OD1','CG','ND2')]
    bonds = [x for x in allbonds if atom in x]
    angles = [x for x in allangles if atom in x] 
    
  elif resn == 'ASP':
    allbonds = bbbonds + [('CA','HA'),('CA','CB'),('CB','HB1'),('CB','HB2'),('CB','CG'),
			  ('CG','OD1'),('CG','OD2'),('OD2','HD2')]
    allangles = bbangles + [('N','CA','HA'),('N','CA','CB'),('CB','CA','HA'),('C','CA','HA'),('CB','CA','C'),('CA','CB','HB1'),
			    ('CA','CB','HB2'),('HB1','CB','HB2'),('CA','CB','CG'),('HB1','CB','CG'),('HB2','CB','CG'),('CB','CG','OD1'),
			    ('CB','CG','OD2'),('OD1','CG','OD2'),('CG','OD2','HD2')]
    bonds = [x for x in allbonds if atom in x]
    angles = [x for x in allangles if atom in x]
    
  elif resn == 'GLN':
    allbonds = bbbonds + [('CA','HA'),('CA','CB'),('CB','HB1'),('CB','HB2'),('CB','CG'),('CG','HG1'),
			  ('CG','HG2'),('CG','CD'),('CD','OE1'),('CD','NE2'),('NE2','HE21'),('NE2','HE22')]
    allangles = bbangles + [('N','CA','HA'),('N','CA','CB'),('CB','CA','HA'),('C','CA','HA'),('CB','CA','C'),('CA','CB','HB1'),
			    ('CA','CB','HB2'),('HB1','CB','HB2'),('CA','CB','CG'),('HB1','CB','CG'),('HB2','CB','CG'),('CB','CG','HG1'),
			    ('CB','CG','HG2'),('HG1','CG','HG2'),('CB','CG','CD'),('HG1','CG','CD'),('HG2','CG','CD'),('CG','CD','OE1'),
			    ('CG','CD','NE2'),('OE1','CD','NE2'),('CD','NE2','HE21'),('CD','NE2','HE22'),('HE21','NE2','HE22')]
    bonds = [x for x in allbonds if atom in x]
    angles = [x for x in allangles if atom in x]
  
  elif resn == 'GLU':
    allbonds = bbbonds + [('CA','HA'),('CA','CB'),('CB','HB1'),('CB','HB2'),('CB','CG'),('CG','HG1'),
			  ('CG','HG2'),('CG','CD'),('CD','OE1'),('CD','OE2'),('OE2','HE2')]
    allangles = bbangles + [('N','CA','HA'),('N','CA','CB'),('CB','CA','HA'),('C','CA','HA'),('CB','CA','C'),('CA','CB','HB1'),
			    ('CA','CB','HB2'),('HB1','CB','HB2'),('CA','CB','CG'),('HB1','CB','CG'),('HB2','CB','CG'),('CB','CG','HG1'),
			    ('CB','CG','HG2'),('HG1','CG','HG2'),('CB','CG','CD'),('HG1','CG','CD'),('HG2','CG','CD'),('CG','CD','OE1'),
			    ('CG','CD','OE2'),('OE1','CD','OE2'),('CD','OE2','HE2')]
    bonds = [x for x in allbonds if atom in x]
    angles = [x for x in allangles if atom in x]
    
  elif resn == 'GLY':
    allbonds = bbbonds + [('CA','HA1'),('CA','HA2')]
    allangles = bbangles + [('N','CA','HA1'),('C','CA','HA1'),('N','CA','HA2'),('C','CA','HA2'),('HA1','CA','HA2')]
    bonds = [x for x in allbonds if atom in x]
    angles = [x for x in allangles if atom in x]
    
  elif resn == 'HIS':
    allbonds = bbbonds + [('CA','HA'),('CA','CB'),('CB','HB1'),('CB','HB2'),('CB','CG'),('CG','ND1'),('ND1','HD1'),
			  ('ND1','CE1'),('CE1','HE1'),('CG','CD2'),('CD2','HD2'),('CD2','NE2'),('NE2','HE2'),('CE1','NE2')]
    allangles = bbangles + [('N','CA','HA'),('N','CA','CB'),('CB','CA','HA'),('C','CA','HA'),('CB','CA','C'),('CA','CB','HB1'),
			    ('CA','CB','HB2'),('HB1','CB','HB2'),('CA','CB','CG'),('HB1','CB','CG'),('HB2','CB','CG'),('CB','CG','ND1'),
			    ('CB','CG','CD2'),('CD2','CG','ND1'),('CG','ND1','HD1'),('CG','ND1','CE1'),('CG','CD2','HD2'),('CG','CD2','NE2'),
			    ('HD1','ND1','CE1'),('ND1','CE1','HE1'),('ND1','CE1','NE2'),('CE1','NE2','HE2'),('HE1','CE1','NE2'),
			    ('HD2','CD2','NE2'),('CD2','NE2','HE2'),('CD2','NE2','CE1')]
    bonds = [x for x in allbonds if atom in x]
    angles = [x for x in allangles if atom in x]
    
  elif resn == 'ILE':
    allbonds = bbbonds + [('CA','HA'),('CA','CB'),('CB','HB'),('CB','CG1'),('CG1','HG11'),('CG1','HG12'),
			  ('CB','CG2'),('CG2','HG21'),('CG2','HG22'),('CG2','HG23'),('CG1','CD1'),('CD1','HD11'),
			  ('CD1','HD12'),('CD1','HD13')]
    allangles = bbangles + [('N','CA','HA'),('N','CA','CB'),('CB','CA','HA'),('C','CA','HA'),('CB','CA','C'),('CA','CB','HB'),('CA','CB','CG1'),
			    ('CA','CB','CG2'),('HB','CB','CG1'),('HB','CB','CG2'),('CG1','CB','CG2'),('CB','CG1','HG11'),('CB','CG1','HG12'),('HG11','CG','HG12'),
			    ('CB','CG2','HG21'),('CB','CG2','HG22'),('CB','CG2','HG23'),('HG21','CG2','HG22'),('HG21','CG2','HG23'),('HG22','CG2','HG23'),
			    ('CB','CG1','CD1'),('HG11','CG1','CD1'),('HG12','CG1','CD1'),('CG1','CD1','HD11'),('CG1','CD1','HD12'),
			    ('CG1','CD1','HD13'),('HD11','CD1','HD12'),('HD11','CD1','HD13'),('HD12','CD1','HD13')]
    bonds = [x for x in allbonds if atom in x]
    angles = [x for x in allangles if atom in x]
    
  elif resn == 'LEU':
    allbonds = bbbonds + [('CA','HA'),('CA','CB'),('CB','HB1'),('CB','HB2'),('CB','CG'),('CG','HG'),('CG','CD1'),
			  ('CD1','HD11'),('CD1','HD12'),('CD1','HD13'),('CG','CD2'),('CD2','HD21'),('CD2','HD22'),('CD2','HD23')]
    allangles = bbangles + [('N','CA','HA'),('N','CA','CB'),('CB','CA','HA'),('C','CA','HA'),('CB','CA','C'),('CA','CB','HB1'),('CA','CB','HB2'),('HB1','CB','HB2'),('CA','CB','CG'),
			    ('HB1','CB','CG'),('HB2','CB','CG'),('CB','CG','CD1'),('CB','CG','CD2'),('CB','CG','HG'),('HG','CG','CD1'),
			    ('HG','CG','CD2'),('CG','CD1','HD11'),('CG','CD1','HD12'),('CG','CD1','HD13'),('HD11','CD1','HD12'),
			    ('HD11','CD1','HD13'),('HD12','CD1','HD13'),('CG','CD2','HD21'),('CG','CD2','HD22'),('CG','CD2','HD23'),
			    ('HD21','CD2','HD22'),('HD21','CD2','HD23'),('HD22','CD2','HD23')]
    bonds = [x for x in allbonds if atom in x]
    angles = [x for x in allangles if atom in x]
    
  elif resn == 'LYS':
    allbonds = bbbonds + [('CA','HA'),('CA','CB'),('CB','HB1'),('CB','HB2'),('CB','CG'),('CG','HG1'),('CG','HG2'),('CG','CD'),
			  ('CD','HD1'),('CD','HD2'),('CD','CE'),('CE','HE1'),('CE','HE2'),('CE','NZ'),('NZ','HZ1'),('NZ','HZ2'),('NZ','HZ3')]
    allangles = bbangles + [('N','CA','HA'),('N','CA','CB'),('CB','CA','HA'),('C','CA','HA'),('CB','CA','C'),('CA','CB','HB1'),('CA','CB','HB2'),('CA','CB','CG'),
			    ('HB1','CB','CG'),('HB2','CB','CG'),('HB1','CB','HB2'),('CB','CG','HG1'),('CB','CG','HG2'),('HG1','CG','HG2'),('CB','CG','CD'),
			    ('HG1','CG','CD'),('HG2','CG','CD'),('CG','CD','CE'),('CG','CD','HD1'),('CG','CD','HD2'),('HD1','CD','HD2'),('HD1','CD','CE'),
			    ('HD2','CD','CE'),('CD','CE','HE1'),('CD','CE','HE2'),('HE1','CE','HE2'),('CD','CE','NZ'),('HE1','CE','NZ'),('HE2','CE','NZ'),
			    ('CE','NZ','HZ1'),('CE','NZ','HZ2'),('CE','NZ','HZ3'),('HZ1','NZ','HZ2'),('HZ1','NZ','HZ3'),('HZ2','NZ','HZ3')]
    bonds = [x for x in allbonds if atom in x]
    angles = [x for x in allangles if atom in x]
    
  elif resn == 'MET':
    allbonds = bbbonds + [('CA','HA'),('CA','CB'),('CB','HB1'),('CB','HB2'),('CB','CG'),('CG','HG1'),('CG','HG2'),('CG','SD'),('SD','CE'),
			  ('CE','HE1'),('CE','HE2'),('CE','HE3')]
    allangles = bbangles + [('N','CA','HA'),('N','CA','CB'),('CB','CA','HA'),('C','CA','HA'),('CB','CA','C'),('CA','CB','HB1'),('CA','CB','HB2'),('HB1','CB','HB2'),('CA','CB','CG'),
			    ('HB1','CB','CG'),('HB2','CB','CG'),('CB','CG','HG1'),('CB','CG','HG2'),('HG1','CG','HG2'),('CB','CG','SD'),('HG1','CG','SD'),('HG2','CG','SD'),
			    ('CG','SD','CE'),('SD','CE','HE1'),('SD','CE','HE2'),('SD','CE','HE3'),('HE1','CE','HE2'),('HE1','CE','HE3'),('HE2','CE','HE3')]
    bonds = [x for x in allbonds if atom in x]
    angles = [x for x in allangles if atom in x]
    
  elif resn == 'PHE':
    allbonds = bbbonds + [('CA','HA'),('CA','CB'),('CB','HB1'),('CB','HB2'),('CB','CG'),('CG','CD1'),('CD1','HD1'),('CG','CD2'), ('CD2','HD2'),
			 ('CD1','CE1'),('CE1','HE1'),('CD2','CE2'),('CE2','HE2'),('CE1','CZ'),('CZ','HZ'),('CE2','CZ')]
    allangles = bbangles + [('N','CA','HA'),('N','CA','CB'),('CB','CA','HA'),('C','CA','HA'),('CB','CA','C'),('CA','CB','HB1'),('CA','CB','HB2'),('HB1','CB','HB2'),('CA','CB','CG'),
			    ('HB1','CB','CG'),('HB2','CB','CG'),('CB','CG','CD1'),('CB','CG','CD2'),('CD1','CG','CD2'),('CG','CD1','HD1'),('CG','CD2','HD2'),('CG','CD1','CE1'),('CG','CD2','CE2'),
			    ('HD1','CD1','CE1'),('HD2','CD2','CE2'),('CD1','CE1','HE1'),('CD2','CE2','HE2'),('CD1','CE1','CZ'),('CD2','CE2','CZ'),('CE1','CZ','HZ'),
			    ('CE2','CZ','HZ'),('CE1','CZ','CE2'),('HE1','CE1','CZ'),('HE2','CE2','CZ')]
    bonds = [x for x in allbonds if atom in x]
    angles = [x for x in allangles if atom in x]
    
  elif resn == 'PRO':
    allbonds = [('N','CA'),('CA','HA'), ('CA','CB'),('CB','HB1'),('CB','HB2'),('CB','CG'),('CG','HG1'),('CG','HG2'),('CG','CD'),
		('CD','HD1'),('CD','HD2'),('CD','N'),('CA','C'),('C','O')]
    allangles = [('N','CA','C'),('CA','C','O'),('N','CA','HA'),('N','CA','CB'),('N','CD','HD1'),('N','CD','HD2'),('N','CD','CG'),
		 ('HA','CA','CB'),('CA','CB','CG'),('CA','CB','HB1'),('CA','CB','HB2'),('HB1','CB','HB2'),('HB1','CB','CG'),('HB2','CB','CG'),
		 ('CB','CG','CD'),('CB','CG','HG1'),('CB','CG','HG2'),('HG1','CG','HG2'),('HG1','CG','CD'),('HG2','CG','CD'),('CG','CD','HD1'),
		 ('CG','CD','HD2'),('HD1','CD','HD2')]
    bonds = [x for x in allbonds if atom in x]
    angles = [x for x in allangles if atom in x]
    
  elif resn == 'SER':
    allbonds = bbbonds + [('CA','HA'), ('CA','CB'),('CB','HB1'),('CB','HB2'),('CB','OG'),('OG','HG')]
    allangles = bbangles + [('N','CA','HA'),('N','CA','CB'),('CB','CA','HA'),('C','CA','HA'),('CB','CA','C'),('CA','CB','HB1'),('CA','CB','HB2'),('HB1','CB','HB2'),
			    ('HB1','CB','OG'),('HB2','CB','OG'),('CA','CB','OG'),('CB','OG','HG')]
    bonds = [x for x in allbonds if atom in x]
    angles = [x for x in allangles if atom in x]
    
  elif resn == 'THR':
    allbonds = bbbonds + [('CA','HA'),('CA','CB'),('CB','HB'),('CB','OG1'),('OG1','HG1'),('CB','CG2'),('CG2','HG21'),('CG2','HG22'),('CG2','HG23')]
    allangles = bbangles + [('N','CA','HA'),('N','CA','CB'),('CB','CA','HA'),('C','CA','HA'),('CB','CA','C'),('CA','CB','HB'),('CA','CB','OG1'),('CA','CB','CG2'),
			    ('HB','CB','OG1'),('HB','CB','CG2'),('CB','OG1','HG1'),('CB','CG2','HG21'),('CB','CG2','HG22'),('CB','CG2','HG23'),('HG21','CG2','HG22'),
			    ('HG21','CG2','HG23'),('HG22','CG2','HG23')]
    bonds = [x for x in allbonds if atom in x]
    angles = [x for x in allangles if atom in x]
    
  elif resn == 'TRP':
    allbonds = bbbonds + [('CA','HA'),('CA','CB'),('CB','HB1'),('CB','HB2'),('CB','CG'),('CG','CD1'),('CD1','HD1'),('CG','CD2'),('CD1','NE1'),('NE1','HE1'),
			  ('NE1','CE2'),('CD2','CE2'),('CD2','CE3'),('CE3','HE3'),('CE2','CZ2'),('CZ2','HZ2'),('CE3','CZ3'),('CZ3','HZ3'),('CZ2','CH2'),('CH2','HH2'),('CZ3','CH2')]
    allangles = bbangles + [('N','CA','HA'),('N','CA','CB'),('CB','CA','HA'),('C','CA','HA'),('CB','CA','C'),('CA','CB','HB1'),('CA','CB','HB2'),('CA','CB','CG'),
			    ('HB1','CB','CG'),('HB2','CB','CG'),('CB','CG','CD1'),('CB','CG','CD2'),('CD1','CG','CD2'),('CG','CD1','HD1'),('CG','CD1','NE1'),('CG','CD2','CE2'),('CG','CD2','CE3'),('CD1','NE1','HE1'),
			    ('HE1','NE1','CE2'),('NE1','CE2','CD2'),('NE1','CE2','CZ2'),('CD2','CE2','CZ2'),('CD2','CE3','HE3'),('CD2','CE3','CZ3'),('CE2','CZ2','HZ2'),('CE2','CD2','CE3'),
			    ('CE2','CZ2','CH2'),('CE2','NE1','CD1'),('HE3','CE3','CZ3'),('CE3','CZ3','HZ3'),('CE3','CZ3','CH2'),('HZ2','CZ2','CH2'),('CZ2','CH2','CZ3'),('HZ3','CZ3','CH2'),('CZ2','CH2','HH2'),
			    ('CZ3','CH2','HH2')]
    bonds = [x for x in allbonds if atom in x]
    angles = [x for x in allangles if atom in x]
    
  elif resn == 'TYR':
    allbonds = bbbonds + [('CA','HA'),('CA','CB'),('CB','HB1'),('CB','HB2'),('CB','CG'),('CG','CD1'),('CD1','HD1'),('CG','CD2'),('CD2','HD2'),('CD1','CE1'),('CE1','HE1'),
			  ('CD2','CE2'),('CE2','HE2'),('CE1','CZ'),('CE2','CZ'),('CZ','OH'),('OH','HH')]
    allangles = bbangles + [('N','CA','HA'),('N','CA','CB'),('CB','CA','HA'),('C','CA','HA'),('CB','CA','C'),('CA','CB','HB1'),('CA','CB','HB2'),('CA','CB','CG'),
			    ('HB1','CB','CG'),('HB2','CB','CG'),('CB','CG','CD1'),('CB','CG','CD2'),('CD1','CG','CD2'),('CG','CD1','HD1'),('CG','CD2','HD2'),('CG','CD1','CE1'),('CG','CD2','CE2'),
			    ('HD1','CD1','CE1'),('CD1','CE1','HE1'),('CD1','CE1','CZ'),('HD2','CD2','CE2'),('CD2','CE2','HE2'),('CD2','CE2','CZ'),('HE1','CE1','CZ'),
			    ('HE2','CE2','CZ'),('CE1','CZ','OH'),('CE1','CZ','CE2'),('CE2','CZ','OH'),('CZ','OH','HH')]
    bonds = [x for x in allbonds if atom in x]
    angles = [x for x in allangles if atom in x]
    
  elif resn == 'VAL':
    allbonds = bbbonds + [('CA','HA'),('CA','CB'),('CB','HB'),('CB','CG1'),('CG1','HG11'),('CG1','HG12'),('CG1','HG13'),('CB','CG2'),('CG2','HG21'),('CG2','HG22'),('CG2','HG23')]
    allangles = bbangles + [('N','CA','HA'),('N','CA','CB'),('CB','CA','HA'),('C','CA','HA'),('CB','CA','C'),('CA','CB','HB'),('HB','CB','CG1'),('HB','CB','CG2'),
			    ('CA','CB','CG1'),('CA','CB','CG2'),('CG1','CB','CG2'),('CB','CG1','HG11'),('CB','CG1','HG12'),('CB','CG1','HG13'),('HG11','CG1','HG12'),('HG11','CG1','HG13'),
			    ('HG12','CG1','HG13'),('CB','CG2','HG21'),('CB','CG2','HG22'),('CB','CG2','HG23'),('HG21','CG2','HG22'),('HG21','CG2','HG23'),('HG22','CG2','HG23')]
    bonds = [x for x in allbonds if atom in x]
    angles = [x for x in allangles if atom in x]
    
  else:
    print "ERROR: residue name not found!"
    sys.exit(1)
    
  return bonds, angles