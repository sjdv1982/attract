      subroutine kick(max3atom,maxatom,maxmode,
     2 natom,npsmgroups,dligp,
     3 spring,nhm0,eig,nrw,rw,rrv,rr
     4 dligp2,ekick)
      implicit none

c     parameters: input
      integer max3atom,maxatom,maxmode
      integer natom,npsmgroups,nhm0,nrw
      real*8 dligp,eig,rw,rrv
      dimension dligp(maxmode)
      dimension eig(maxmode,max3atom)
      dimension rw(maxatom,max3atom)
      dimension rrv(maxmode,maxmode,3)

c     parameters: output
      real*8 dligp2
      dimension dligp2(maxmode)
      real*8 ekick

c     local variables
      integer i,j,iii
      integer isel,idir,imode
      real*8 ampl,ampl_max!,rd,displ,cl,rw,rr

c initial max ampl ~sqrt(#active atoms)
      ampl_max=5.0
c random kick center atom
      isel=int(rand(0)*natom)+1
c random kick direction
      idir=int(rand(0)*2)*2-1
c random mode number
      imode=int(rand(0)*nhm0)+1
c kick amplitude
      ampl=rand(0)*ampl_max
      
c do the actual kick
      do 100 i=1,npsmgroups
c     generic displacement:
      displ=ampl*idir*rw(isel,i)
      iii=3*(i-1)
c     deformation:
      dligp2(iii+1)=dligp(iii+1)+displ*eig(imode,iii+1)
      dligp2(iii+2)=dligp(iii+2)+displ*eig(imode,iii+2)
      dligp2(iii+3)=dligp(iii+3)+displ*eig(imode,iii+3)
  100 continue

      call spring_energy(maxmode,npsmgroups,dligp2,spring,
     1 rrv,rr,ekick)

      return
      end
     
      