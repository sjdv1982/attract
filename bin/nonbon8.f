      subroutine nonbon8(maxatom,maxmolpair,
     1 iab,xl,xr,fl,fr,wel,wer,chair,chail,ac,rc,       
     2 emin,rmin2,iacir,iacil,nonr,nonl,ipon,nonp,
     3 potshape, swi_on,swi_off, enon,epote)
      implicit none

c     Parameters      
      integer maxatom, maxmolpair,iab,nonp,potshape
      real swi_on, swi_off
      real*8 xl,xr,fl,fr,wel,wer,chair,chail,ac,rc
      real*8 emin,rmin2,enon,epote      
      integer iacir,iacil,ipon,nonl,nonr
      dimension nonr(maxmolpair),nonl(maxmolpair),ipon(99,99),
     1  iacir(maxatom), iacil(maxatom)
      dimension chair(maxatom),chail(maxatom),wer(maxatom),wel(maxatom),
     1  fl(maxatom),fr(maxatom),ac(99,99),rc(99,99),rmin2(99,99)
      dimension xl(maxatom),xr(maxatom),emin(99,99)

c     Local variables
      real*8 dx,xnull,r2,rr2,rrd,rr23,rep,vlj,et,charge,alen,rlen,fb,fdb
      real*8 fswi, r
      integer k,ik,i,j,ii,jj,it,jt,ivor
      dimension dx(3)
      
      xnull=0.0d0
      enon=xnull
      epote=xnull
      do 100 ik=1,nonp
      i=nonr(ik)
      j=nonl(ik)
      it=iacir(i)
      jt=iacil(j)
      ii=3*(i-1)
      jj=3*(j-1)
      alen=wel(j)*wer(i)*ac(it,jt)
      rlen=wel(j)*wer(i)*rc(it,jt)
      ivor=ipon(it,jt)
      charge=wel(j)*wer(i)*chair(i)*chail(j)
      r2=xnull
      do 120 k=1,3
      dx(k)=xl(jj+k)-xr(ii+k)
      r2=r2+dx(k)**2
  120 continue
      if(r2.lt.0.001d0) r2=0.001d0

      fswi = 1.0d0
      if (swi_on.gt.0.or.swi_off.gt.0) then        
	if (r2.ge.(swi_on*swi_on)) then
	  if (r2.ge.(swi_off*swi_off)) then
	    fswi = 0.0d0
	  else
	    r = sqrt(r2) 
	    fswi = 1.0d0-(r - swi_on)/(swi_off-swi_on)
	  endif
	endif
      endif

      rr2=1.0d0/r2
      do 125 k=1,3
      dx(k)=rr2*dx(k)
  125 continue
      et=xnull
      if(charge.gt.0.001.or.charge.lt.-0.001) then
      et=charge*rr2
      epote=epote+fswi*et
      if(iab.eq.1) then
      do 130 k=1,3
      fdb=fswi*2.0d0*et*dx(k)
      fl(jj+k)=fl(jj+k)+fdb
      fr(ii+k)=fr(ii+k)-fdb
  130 continue
      endif
      endif
      rr23=rr2**3
      if (potshape.eq.8) then
      rrd = rr2
      else if (potshape.eq.12) then
      rrd = rr23
      endif
      rep=rlen*rrd
      vlj=(rep-alen)*rr23
      if(r2.lt.rmin2(it,jt)) then
      enon=enon+fswi*(vlj+(ivor-1)*emin(it,jt))
c      write(*,*)'pair',i,j,it,jt,r2,vlj+(ivor-1)*emin(it,jt),et
      if(iab.eq.1) then
      fb=6.0D0*vlj+2.0D0*(rep*rr23)
      do 135 k=1,3
      fdb=fswi*fb*dx(k)
      fl(jj+k)=fl(jj+k)+fdb
      fr(ii+k)=fr(ii+k)-fdb
  135 continue
      endif
      else
      enon=enon+fswi*ivor*vlj
c      write(*,*)'pair',i,j,it,jt,r2,ivor*vlj,et
      if(iab.eq.1) then
      fb=6.0D0*vlj+2.0D0*(rep*rr23)
      do 145 k=1,3
      fdb=fswi*ivor*fb*dx(k)
      fl(jj+k)=fl(jj+k)+fdb
      fr(ii+k)=fr(ii+k)-fdb
  145 continue      
      endif
      endif
  100 continue
      return
      end
