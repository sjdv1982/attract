      subroutine read_parameters(paramfile,rbc,rc,ac,emin,rmin2,ipon,
     1 haspar, potshape, swi_on, swi_off)
      implicit none
      
c     Parameters
      character*100 paramfile
      real*8 rbc,rc,ac,emin,rmin2
      dimension rbc(99,99),rc(99,99),ac(99,99),
     1 emin(99,99),rmin2(99,99)
      integer ipon
      dimension ipon(99,99)
      integer potshape
      real swi_on, swi_off
      integer haspar
      dimension haspar(99,99)

c     Local variables
      integer i,j
      real*8 abc, iflo      
      dimension abc(99,99), iflo(99,99)
      integer maxpar
      
      open(11,file=paramfile)
      read(11,*) potshape,maxpar, swi_on, swi_off
      do 443 i=1,maxpar
      read(11,*)(rbc(i,j),j=1,maxpar)
      do 445 j=1,maxpar
      rbc(i,j)=1.0d0*rbc(i,j)
  445 continue
c     write(*,*)(rbc(i,j),j=1,maxpar)
  443 continue
      do 446 i=1,maxpar
      read(11,*)(abc(i,j),j=1,maxpar)
c     write(*,*)(abc(i,j),j=1,maxpar)
  446 continue
      do 447 i=1,maxpar
      read(11,*)(iflo(i,j),j=1,maxpar)
  447 continue
      close(11)
      do 450 j=1,maxpar
      do 448 i=1,maxpar
      if (rbc(i,j).ne.0) then
      haspar(i,j) = 1
      endif
      rc(i,j)=abc(i,j)*rbc(i,j)**potshape
      ac(i,j)=abc(i,j)*rbc(i,j)**6
      emin(i,j)= 0.0d0
      rmin2(i,j)= 0.0d0
      ipon(i,j)= 1.0d0       
      if (potshape.eq.8) then
       if (ac(i,j).gt.0.and.rc(i,j).gt.0) then
        emin(i,j)=-27.0d0*ac(i,j)**4/(256.0d0*rc(i,j)**3)
        rmin2(i,j)=4.0d0*rc(i,j)/(3.0d0*ac(i,j))
       endif 
       ipon(i,j)=iflo(i,j)
      elseif (potshape.eq.12) then
      emin(i,j)=-0.25d0*abc(i,j)
      rmin2(i,j)=1.25992105*rbc(i,j)*rbc(i,j)
      ipon(i,j)=iflo(i,j)
      else
        write(*,*), 'Unknown potential shape', potshape
        stop 
      endif            
  448 continue
  450 continue
      return
      end
