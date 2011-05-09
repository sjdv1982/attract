      subroutine ligmin(maxlig,maxdof,maxmode,maxatom,
     1                  f,natom,ijk,eig,nhm,delta)
      implicit none
      
c     Parameters      
      integer maxlig,maxdof,maxmode,maxatom
      integer natom,ijk
      real*8 f,eig,delta
      integer nhm
      dimension delta(maxdof),
     1          f(3*maxatom),eig(maxlig,maxmode,3*maxatom)


c     Local variables
      real*8 xnull, force
      integer i,j

c
c  delta(6+i): skalar product between force f and eigenvector I
c  dei(i): skalar product between difference(actual pos.-ori-pos.) and eig(i).
      xnull=0.0d0 
      do 15 i=1,nhm
      delta(6+i)=xnull
   15 continue
      do 20 i=1,nhm      
      do 30 j=1,natom
      force = f(j)*eig(ijk+1,i,j)
      if (force.gt.1.or.force.lt.-1) then
      endif
      
      delta(6+i)=delta(6+i)-force
   30 continue
   20 continue

      return 
      end
