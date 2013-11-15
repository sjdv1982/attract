      subroutine ligmin_index(f,natom,ijk,eig,eig_val,nhm,nihm,delta)
      implicit none

c     Parameters
      include "max.fin"
      integer natom,ijk
      integer eig
      dimension eig(maxlig,maxindexmode,maxlenindexmode)
      real*8 f,eig_val,delta
      integer nhm, nihm
      dimension delta(maxdof),
     1        f(3*maxatom),eig_val(maxlig,maxindexmode,maxlenindexmode)


c     Local variables
      real*8 xnull, force
      integer i,j, jj
      integer, parameter :: ERROR_UNIT = 0
c
c  delta(6+i): skalar product between force f and eigenvector I
c  dei(i): skalar product between difference(actual pos.-ori-pos.) and eig(i).
      do 20 i=1,nihm
      force = 0
      do 30 j=1,maxlenindexmode
c     Only use relevant mode entries
      jj = eig(ijk,i,j)
      if (jj.eq.-1) then
      exit
      else
      force = f(jj+1)*eig_val(ijk,i,j)
      delta(6+nhm+i)=delta(6+nhm+i)-force
      endif

   30 continue
c      write(ERROR_UNIT,*) "Indexdelta", i, delta(6+nhm+i)
   20 continue

      return
      end
