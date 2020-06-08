      subroutine nonbon_soft(iab,xl,xr,fl,fr,wel,wer,chair,chail,ac,rc,
     1 emin,rmin2,iacir,iacil,nonr,nonl,ipon,nonp,
     2 potshape, cdie, swi_on,swi_off, enon,epote, softcore)
      implicit none

c     Parameters
      include "max.fin"
      integer iab,nonp,potshape
      integer cdie
      real swi_on, swi_off
      real*8 xl,xr,fl,fr,wel,wer,chair,chail,ac,rc
      real*8 emin,rmin2,enon,epote      
      integer iacir,iacil,ipon,nonl,nonr
      dimension nonr(maxmolpair),nonl(maxmolpair),ipon(99,99),
     1  iacir(maxatom), iacil(maxatom)
      dimension chair(maxatom),chail(maxatom),wer(maxatom),
     1 wel(maxatom),fl(maxatom),fr(maxatom),ac(99,99),rc(99,99),
     2 rmin2(99,99)
      dimension xl(maxatom),xr(maxatom),emin(99,99)

c MANU: Softcore Variable einführen
      real*8 softcore
      

c     Local variables
      real*8 dx,xnull,r2,rr1,rr2,rrd,rr23,rep,vlj,et,charge,alen,
     1 rlen,fb,fdb, rr2a, att, fb1, fb2
      real*8 fswi, r, shapedelta
      integer k,ik,i,j,ii,jj,it,jt,ivor
      dimension dx(3)
      real*8 e_min
c MANU: Alpha, beta einführen, beta ist softcore parameter für elektrostatik
      real*8 alpha, beta
      
      xnull=0.0d0
      enon=xnull
      epote=xnull
      r2=xnull
      do 100 ik=1,nonp
      i=nonr(ik)
      j=nonl(ik)
      it=iacir(i)
      jt=iacil(j)
      ii=3*(i-1)
      jj=3*(j-1)
      alen=wel(j)*wer(i)*ac(it,jt)
      rlen=wel(j)*wer(i)*rc(it,jt)
      e_min=wel(j)*wer(i)*emin(it,jt)
      ivor=ipon(it,jt)
      charge=wel(j)*wer(i)*chair(i)*chail(j)
      r2=xnull
      do 120 k=1,3
      dx(k)=xl(jj+k)-xr(ii+k)
      r2=r2+dx(k)**2
  120 continue
      if(r2.lt.0.001d0) r2=0.001d0

c MANU Beta berechnen

      fswi = 1.0d0
      if(charge.eq.0) then
       beta=0.0d0
      else
       beta = (3*abs(charge))/softcore
      endif
      if (swi_on.gt.0.or.swi_off.gt.0) then        
	if (r2.ge.(swi_on*swi_on)) then
	  if (r2.ge.(swi_off*swi_off)) then
	    fswi = 0.0d0
            beta = 0.0d0
	  else
	   r = sqrt(r2) 
	   fswi = 1.0d0-(r - swi_on)/(swi_off-swi_on)
           if(charge.eq.0) then
            beta = 0.0d0
           else
            beta = (3*abs(charge))/softcore
           endif
	  endif
	endif
      endif

      rr2=1.0d0/r2
      do 125 k=1,3
      dx(k)=rr2*dx(k)
  125 continue
      et=xnull
      if(charge.gt.0.001.or.charge.lt.-0.001) then
c MANU: if(cdie...) rausschmeißen? Beta als Softcore-Parameter für Softcore-Elektrostatik einführen
      if (cdie.eq.1) then
       rr1 = 1.0d0/(sqrt(r2)+beta)-1.0/(50.0d0+beta) 
c      (cap all distances at 50 A)
       if (rr1.lt.0) then 
       rr1 = 0
       endif
c MANU FRAGE: Soll hier und bei rr2a 1/50.0 auch in "Radius" stehen, wird ergebnis nicht verfäslcht?
       et=charge*rr1
c       write(*,*), sqrt(r2), et, epote
      else
c Auch hier Beta als SOftcore-Parameter für Softcore-Elektrostatik einführen.
       rr2a = 1.0d0/(r2+beta) - (1.0d0/(50.0d0**2+beta))
       if (rr2a.lt.0) then 
       rr2a = 0
       endif
c      (cap all distances at 50 A)
       et=charge*rr2a
      endif
      epote=epote+fswi*et
      if(iab.eq.1) then
      do 130 k=1,3
      if (cdie.eq.1) then
      if (rr1.le.0) then
      fdb = fswi * et * dx(k)
      else
      fb = fswi*charge*(rr1+1.0d0/(50.0d0+beta))**2*sqrt(r2)
      fdb=fb*dx(k)

c      write(*,*) "fdbe1 = ", fdb

      endif
      else
      if (rr2a.le.0) then
      fdb=fswi*2.0d0*et*dx(k)
      else
      fb =  fswi*2.0d0*charge*r2*(rr2a+1.0d0/(50.0d0**2+beta))**2
      fdb = fb*dx(k)

c      write(*,*) "fdbe2 = ", fdb

      endif
      endif
      fl(jj+k)=fl(jj+k)+fdb
      fr(ii+k)=fr(ii+k)-fdb
  130 continue
      endif
      endif
      rr23=rr2**3
c MANU: if(potshape...) rausgeworfen
c MANU: softcore ausrechnen, NOCH ZU ÜBERPRÜFEN: + oder - bei mitternachtsformel
	alpha = (2*rlen)/(alen+sqrt(alen**2+4*rlen*softcore))

c	(softcore-parameter alpha)
c	softcore-Nenner rrd
       rrd = 1.0d0/(r2**3+alpha)
       shapedelta = 6.0D0

       rep=rlen*(rrd**2)
       att = alen*rrd
       vlj = rep - att

c MANU: Größen ausgeben zur Kontrolle:
  
c       write(*,*) 'softcore = ', softcore
c       write(*,*) 'alpha = ', alpha
c       write(*,*) 'alen = ', alen
c       write(*,*) 'rlen = ', rlen
c       write(*,*) '##vlj = ', vlj

c       write(*,*) "epote = ", epote
c       write(*,*) "fswi = ", fswi
c       write(*,*) "charge = ", charge
c       write(*,*) "beta = ", beta

      if(r2.lt.rmin2(it,jt)) then
      enon=enon+fswi*(vlj+(ivor-1)*e_min)
!       write(*,*)'pair',i,j,it,jt,r2,vlj+(ivor-1)*e_min,e_min
c      write(*,*)'pair',i,j,it,jt,r2,vlj+(ivor-1)*emin(it,jt),et,
c     1 emin(it,jt)
      if(iab.eq.1) then
c	geändert
c	original: shapedelta*(rep*rr23)
      fb=(6.0D0*vlj+shapedelta*rep)*rrd*r2**3
      do 135 k=1,3

      fdb=fswi*fb*dx(k)

c      write(*,*) "fdb1 = ", fdb
      
      fl(jj+k)=fl(jj+k)+fdb
      fr(ii+k)=fr(ii+k)-fdb
  135 continue
      endif
      else
      enon=enon+fswi*ivor*vlj
c      write(*,*)'pair',i,j,it,jt,r2,ivor*vlj,et,
c     1 emin(it,jt)
      if(iab.eq.1) then
c	geändert
c	original: shapedelta*(rep*rr23)
      fb=(6.0D0*vlj+shapedelta*rep)*rrd*r2**3

      do 145 k=1,3

      fdb=fswi*ivor*fb*dx(k)
      fl(jj+k)=fl(jj+k)+fdb
      fr(ii+k)=fr(ii+k)-fdb

c      write(*,*) "fdb2 = ", fdb

  145 continue      
      endif
      endif
  100 continue
c      write(ERROR_UNIT, *) "nonbon_soft ", enon, epote

      return
      end
