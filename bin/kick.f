      subroutine kick(max3atom,maxatom,maxmode,
     2 natom,xori0,iaci,dligp,kickwork,nhm0,eig,dligp2,ekick)
      implicit none

c     parameters
      integer max3atom,maxatom,maxmode
      integer natom

      integer nhm0
      real*8 ekick
      
      integer iaci
      dimension iaci(maxatom)
      real*8 xori0,dligp,kickwork,eig,dligp2
      dimension xori0(max3atom),
     1  dligp(maxmode),eig(maxmode,max3atom),
     2  kickwork(maxatom,maxatom),dligp2(maxmode)

!
!
!
!
!

      return
      end
