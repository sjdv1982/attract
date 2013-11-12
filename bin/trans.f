      subroutine trans(f,delta,natom)
      implicit none
c     Parameters
      include 'max.fin'
      integer natom
      real*8 f,delta
      integer, parameter :: max3atom = 3*maxatom
      dimension f(max3atom),delta(maxdof)
      
c     Local variables
      real*8 flim,ftr1,ftr2,ftr3,fbetr
      integer i,ii
      
c     In this subroutine the translational force components are calculated
      flim=1.0d18
      ftr1=0.0d0
      ftr2=0.0d0
      ftr3=0.0d0
      do 12 i=1,natom
      ii=3*(i-1)
      ftr1=ftr1+f(ii+1)
      ftr2=ftr2+f(ii+2)
      ftr3=ftr3+f(ii+3)
c      write(*,*),'ftr',ii,f(ii+1),f(ii+2),f(ii+3)
   12 continue
c force reduction, some times helps in case of very "bad" start structure
      do 20 i=1,3
      fbetr=ftr1**2+ftr2**2+ftr3**2
      if(fbetr.gt.flim) then
      ftr1=.01d0*ftr1
      ftr2=.01d0*ftr2
      ftr3=.01d0*ftr3
      end if
   20 continue
      delta(4)=ftr1
      delta(5)=ftr2
      delta(6)=ftr3
c      write(*,*),'trans',ftr1,ftr2,ftr3
      return
      end
