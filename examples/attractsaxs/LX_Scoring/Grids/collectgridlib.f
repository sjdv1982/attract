      subroutine stepgrid(rcut, coorrec, rlen, coorlig, llen,
     1 atypsrec, atypslig, bins, gtypes, grid)
cf2py intent(out) :: grid
cf2py intent(hide) :: rlen
cf2py intent(hide) :: llen
cf2py intent(hide) :: bins
      integer bins
      real rcut(bins)
      integer rlen, gtypes, atyps
      integer llen
      integer indexmap(100,100)
      integer i, j, b
      real coorrec(rlen,3)
      real coorlig(llen,3)
      integer grid(bins, gtypes)
      integer atypsrec(rlen)
      integer atypslig(llen)
      real xr, yr, zr, xl, yl, zl, dist2, rend2, r2
      integer ar, al, s
c gtypes is the number of parameter: gtypes = atyps*(atpys-1)/2 +atyps
c modatyp is atyps calculated from gtypes
      atyps=int(sqrt(0.25+2.*gtypes)-0.5)
      ind=1
      do 10, i=1, atyps
       do 20, j=i, atyps
        indexmap(i,j)=ind
        indexmap(j,i)=ind
        ind=ind+1
   20  continue
   10 continue
      rend2=rcut(bins)
      do 100, i=1, rlen
       xr=coorrec(i,1)
       yr=coorrec(i,2)
       zr=coorrec(i,3)
       ar=atypsrec(i)
       do 200, j=1, llen
        xl=coorlig(j,1)
        yl=coorlig(j,2)
        zl=coorlig(j,3)
        al=atypslig(j)
        dist2=(xr-xl)**2+(yr-yl)**2+(zr-zl)**2
        if (dist2 .le. rend2) then
	 do 300, b=1, bins
         r2=rcut(b)
         if (dist2 .le. r2) then
          in=indexmap(ar,al)
          grid(b,in)= grid(b,in)+1
          exit
         endif
  300   continue
        endif
  200  continue
  100 continue
      return
      end
c      
      subroutine gridnone(stepsize,start,rcut, steps, coorrec, rlen, 
     1 coorlig, llen, atypsrec, atypslig, gtypes, grid)
cf2py intent(out) :: grid
cf2py intent(hide) :: rlen
cf2py intent(hide) :: llen
      real rcut, stepsize, start
      integer steps, gtypes
      integer rlen
      integer llen
      integer i, j, ind, modatyp, pos
      integer indexmap(100,100)
      real coorrec(rlen,3)
      real coorlig(llen,3)
      integer atypsrec(rlen)
      integer atypslig(llen)
      integer grid(steps,gtypes)
      real xr, yr, zr, xl, yl, zl, dist2
      integer ar, al, s
c gtypes is the number of parameter: gtypes = atyps*(atpys-1)/2 +atyps
c modatyp is atyps calculated from gtypes
      modatyp=int(sqrt(0.25+2.*gtypes)-0.5)
      ind=1
      do 10, i=1, modatyp
       do 20, j=i, modatyp
        indexmap(i,j)=ind
        indexmap(j,i)=ind
        ind=ind+1
   20  continue
   10 continue
      do 100, i=1, rlen
       xr=coorrec(i,1)
       yr=coorrec(i,2)
       zr=coorrec(i,3)
       ar=atypsrec(i)
       do 200, j=1, llen
        xl=coorlig(j,1)
        yl=coorlig(j,2)
        zl=coorlig(j,3)
        al=atypslig(j)
        dist2=(xr-xl)**2+(yr-yl)**2+(zr-zl)**2
        if (dist2 .le. rcut**2) then
         pos=indexmap(ar,al)         
         if (dist2 .lt. start**2) then
          dist2=start**2
         endif
         s=int((sqrt(dist2)-start)/stepsize)+1
         grid(s,pos)=grid(s,pos)+1
        endif
  200  continue
  100 continue
      return
      end
c      
c
      subroutine spline(ninterpol, stepsize, start, ending, potsteps, 
     1 coorrec, rlen, coorlig, llen, atypsrec, atypslig, gtypes, grid)
cf2py intent(out) :: grid
cf2py intent(hide) :: rlen
cf2py intent(hide) :: llen
      real*8 stepsize, start, ending
      integer potsteps, gtypes
      integer rlen
      integer llen
      integer i, j, ind, modatyp, pos
      integer indexmap(100,100)
      real*8 coorrec(rlen,3)
      real*8 coorlig(llen,3)
      integer atypsrec(rlen)
      integer atypslig(llen)
      real*8 grid(potsteps,gtypes), nodes(potsteps,ninterpol)
      real*8 xr, yr, zr, xl, yl, zl, dist2, c, dist
      integer ar, al, s, upedge, edge, sg
c gtypes is the number of parameter: gtypes = atyps*(atpys-1)/2 +atyps
c modatyp is atyps calculated from gtypes
      upedge = potsteps - ninterpol +1
      do 30, i=1, potsteps
       if (i .lt. upedge) then
        edge = i
       else
        edge = potsteps - ninterpol +1
       endif
       do 40, j=1, ninterpol
        nodes(i,j)=start+(edge-1)*stepsize+(j-1)*stepsize
c        write(*,*) i,j,nodes(i,j)
   40  continue
   30 continue  
      modatyp=int(sqrt(0.25+2.*gtypes)-0.5)
      ind=1
      do 10, i=1, modatyp
       do 20, j=i, modatyp
        indexmap(i,j)=ind
        indexmap(j,i)=ind
        ind=ind+1
   20  continue
   10 continue
      do 100, i=1, rlen
       xr=coorrec(i,1)
       yr=coorrec(i,2)
       zr=coorrec(i,3)
       ar=atypsrec(i)
       do 200, j=1, llen
        xl=coorlig(j,1)
        yl=coorlig(j,2)
        zl=coorlig(j,3)
        al=atypslig(j)
        dist2=(xr-xl)**2+(yr-yl)**2+(zr-zl)**2
        if (dist2 .le. ending**2) then
         pos=indexmap(ar,al)
         dist = sqrt(dist2)
         if (dist .lt. start) then
          dist=start
         endif
         s=int((dist-start)/stepsize)+1
         if (s .lt. upedge) then
	  sg=int((dist-start)/stepsize)+1
	 else
	  sg= upedge
	 endif
         do 300, k=1, ninterpol
          c = 1.0d0
          do 400, l=1, ninterpol
	   if (l .ne. k) then
	    c = c*((dist-nodes(s,l))/(nodes(s,k)-nodes(s,l)))
	   endif
  400     continue
         grid(sg+k-1,pos)=grid(sg+k-1,pos)+c
  300    continue
        endif  
  200  continue
  100 continue
      return
      end
c      
c      
c
      subroutine distances(power, npower,stepsize,start, rcut, steps, 
     1 coorrec, rlen, coorlig, llen, atypsrec, atypslig, gtypes, grid)
cf2py intent(out) :: grid
cf2py intent(hide) :: rlen
cf2py intent(hide) :: llen
cf2py intent(hide) :: numpower
      real*8 rcut, stepsize, start, power(npower)
      integer steps, gtypes
      integer rlen
      integer llen
      integer i, j, ind, modatyp, pos
      integer indexmap(100,100)
      real*8 coorrec(rlen,3)
      real*8 coorlig(llen,3)
      integer atypsrec(rlen)
      integer atypslig(llen)
      real*8 grid(steps, npower, gtypes)
      real*8 xr, yr, zr, xl, yl, zl, dist2, p
      integer ar, al, s, k
c gtypes is the number of parameter: gtypes = atyps*(atpys-1)/2 +atyps
c modatyp is atyps calculated from gtypes
      modatyp=int(sqrt(0.25+2.*gtypes)-0.5)
      ind=1
      do 10, i=1, modatyp
       do 20, j=i, modatyp
        indexmap(i,j)=ind
        indexmap(j,i)=ind
        ind=ind+1
   20  continue
   10 continue
      do 100, i=1, rlen
       xr=coorrec(i,1)
       yr=coorrec(i,2)
       zr=coorrec(i,3)
       ar=atypsrec(i)
       do 200, j=1, llen
        xl=coorlig(j,1)
        yl=coorlig(j,2)
        zl=coorlig(j,3)
        al=atypslig(j)
        dist2=(xr-xl)**2+(yr-yl)**2+(zr-zl)**2
        if (dist2 .le. rcut**2) then
         pos=indexmap(ar,al)
         if (dist2 .lt. start**2) then
          dist2=start**2
         endif
         s=int((sqrt(dist2)-start)/stepsize)+1
         do 300, k = 1, npower
         p = power(k)
         grid(s,k,pos)=grid(s,k,pos)+1.0d0/dist2**p
  300    continue
         endif
  200  continue
  100 continue
      return
      end