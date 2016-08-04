import numpy as np
import matplotlib.pyplot as plt
import sys

def fgen(stepsize, start, rcutoff, steps, polfunc, potpars, power1, power2, indmax, change):

    x=np.arange(steps)*stepsize+start
    
    #normal attract potential without attractive shift
    if polfunc==0:
	def function(r,ai,bi,rmini,emini,ivori):
	  if r<=rmini:
	    f=(ai/(r**power1))-(bi/(r**power2))+(ivori-1)*emini  
	  elif r>=rmini:
	      f=ivori*((ai/(r**(power1)))-(bi/(r**power2)))
	  return f
	
	vecfunc=np.vectorize(function)

	if change==-1:
	    A=potpars[1]*potpars[0]**power1
	    B=potpars[1]*potpars[0]**power2
	    Emin=((power2/power1)**(power1/(power1-power2))-(power2/power1)**(power2/(power1-power2)))*potpars[1] #only for power=8!!!!!!!!!!!!!1
	    Rmin=((power1/power2)**(1./(power1-power2)))*potpars[0]
	    Ivor=potpars[2]
	    scores=np.zeros((steps,indmax))
	    for i in range(steps):
		scores[i,:]=vecfunc(x[i],A,B,Rmin,Emin,Ivor)
	    return scores
	else:
	    a=abs(potpars[1])*potpars[0]**power1
	    b=abs(potpars[1])*potpars[0]**power2
	    emin=((power2/power1)**(power1/(power1-power2))-(power2/power1)**(power2/(power1-power2)))*abs(potpars[1])
	    rmin=((power1/power2)**(1./(power1-power2)))*potpars[0]
	    ivor=np.sign(potpars[1])*potpars[2]
	    changedscores=vecfunc(x, a, b, rmin, emin, ivor)
	    return changedscores
    
    elif polfunc==1:
	
	def function(r, epsi, sigi, hi, rmini, emini, ivori):
	  if r<=rmini:
	    f=epsi*((sigi**power1/r**(power1))-(sigi**power2/r**power2))+(ivori-1)*emini+(ivori*hi*np.exp(-(r-((2.*rmini)-sigi))**2/(rmini-sigi)**2)) 
	  elif r>=rmini:
	      f=ivori*epsi*((sigi**power1/r**(power1))-(sigi**power2/r**power2))+(ivori*hi*np.exp(-(r-((2.*rmini)-sigi))**2/(rmini-sigi)**2))
	  return f
	
	vecfunc=np.vectorize(function)

	if change==-1:
	    Eps=potpars[1]
	    Sig=potpars[0]
	    H=potpars[2]
	    Emin=((power2/power1)**(power1/(power1-power2))-(power2/power1)**(power2/(power1-power2)))*potpars[1] #only for power=8!!!!!!!!!!!!!1
	    Rmin=((power1/power2)**(1./(power1-power2)))*potpars[0]
	    Ivor=potpars[3]
	    scores=np.zeros((steps,indmax))
	    for i in range(steps):
		scores[i,:]=vecfunc(x[i],Eps,Sig,H,Rmin,Emin,Ivor)
	    return scores
	else:
	    eps=abs(potpars[1])
	    sig=potpars[0]
	    h=potpars[2]
	    emin=((power2/power1)**(power1/(power1-power2))-(power2/power1)**(power2/(power1-power2)))*abs(potpars[1])
	    rmin=((power1/power2)**(1./(power1-power2)))*potpars[0]
	    ivor=np.sign(potpars[1])*potpars[3]
	    changedscores=vecfunc(x, eps, sig, h, rmin, emin, ivor)
	    return changedscores
    
    
    elif polfunc==2:
	#attract potential with possible shift of the attractive potentials
	def function(r,a,b,rmini,emini,ivori,rmax):
	  if r<=rmini:
	    f=(a/(r**power1))-(b/(r**power2))+(ivori-1)*emini
	  elif r<=(rmini+rmax):
	    if ivori==1:
	      f=emini
	    else:
	      f=ivori*((a/(r**power1))-(b/(r**power2)))
	  else:
	    if ivori==1:
	      f=(a/((r-rmax)**power1))-(b/((r-rmax)**power2))
	    else:
	      f=ivori*((a/(r**power1))-(b/(r**power2)))
	  return f
    
    
    	vecfunc=np.vectorize(function)

	if change==-1:
	    A=potpars[1]*potpars[0]**power1
	    B=potpars[1]*potpars[0]**power2
	    Rmax=potpars[2]
	    Emin=((power2/power1)**(power1/(power1-power2))-(power2/power1)**(power2/(power1-power2)))*potpars[1] #only for power=8!!!!!!!!!!!!!1
	    Rmin=((power1/power2)**(1./(power1-power2)))*potpars[0]
	    Ivor=potpars[3]
	    scores=np.zeros((steps,indmax))
	    for i in range(steps):
		scores[i,:]=vecfunc(x[i],A,B,Rmin,Emin,Ivor,Rmax)
	    return scores
	else:
	    a=abs(potpars[1])*potpars[0]**power1
	    b=abs(potpars[1])*potpars[0]**power2
	    rmax=potpars[2]
	    emin=((power2/power1)**(power1/(power1-power2))-(power2/power1)**(power2/(power1-power2)))*abs(potpars[1])
	    rmin=((power1/power2)**(1./(power1-power2)))*potpars[0]
	    ivor=np.sign(potpars[1])*potpars[3]
	    changedscores=vecfunc(x, a, b, rmin, emin, ivor, rmax)
	    return changedscores
    
      