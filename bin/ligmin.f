      subroutine ligmin(f,natom,ijk,eig,nhm,delta)
      implicit none
      
c     Parameters
      include "max.fin"
      integer natom,ijk
      real*8 f,eig,delta
      integer nhm
      dimension delta(maxdof),
     1          f(3*maxatom),eig(maxlig,maxmode,3*maxatom)


c     Local variables
      real*8 xnull, force
      integer i,j
      integer, parameter :: ERROR_UNIT = 0

c
c  delta(6+i): skalar product between force f and eigenvector I
c  dei(i): skalar product between difference(actual pos.-ori-pos.) and eig(i).
      do 20 i=1,nhm      
      do 30 j=1,3*natom
      force = f(j)*eig(ijk,i,j)
      if (force.gt.1.or.force.lt.-1) then
      endif
      
      delta(6+i)=delta(6+i)-force
   30 continue
c      write(ERROR_UNIT,*) "Mode delta", i, delta(6+i)
   20 continue

      return 
      end
