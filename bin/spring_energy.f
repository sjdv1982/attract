      subroutine spring_energy(maxmode,npsmgroups,dligp,spring,
     1 rrv,rr,espring)
      implicit none      
c calculate elastic energy for the deformed structure
c TODO: calculate gradients as well
c
c     parameters
      integer maxmode
      integer npsmgroups
      real*8 dligp
      dimension dligp(maxmode)
      real*8 spring
      dimension spring(maxatom,maxatom)
      real*8 rrv
      dimension rrv(maxmode,maxmode,3)                  
      real*8 rr
      dimension rr(maxmode,maxmode)                        
c     output 
      real*8 espring

c     local variables
      integer i,j,iii,jjj
      real*8 rd, rdx, rdy, rdz, 
      real*8 d, dd,dx, dy, dz
      real*8 refe
      
      espring=0
      do 200 i=1,npsmgroups-1
      iii=3*(i-1)
      do  150 j=i+1,npsmgroups
      if (spring(i,j).gt.0.000005) then
      jjj=3*(j-1)
      rdx = dligp(jjj+1)-dligp(iii+1)
      rdy = dligp(jjj+2)-dligp(iii+2)
      rdz = dligp(jjj+3)-dligp(iii+3)
      dx = rdx+rrv(i,j,1)
      dy = rdy+rrv(i,j,2)
      dz = rdz+rrv(i,j,3)
      d = sqrt(dx**2 + dy**2 + dz**2)
c precomputed:  rr(i,j) = sqrt(rrv(i,j,1)**2+rrv(i,j,2)**2+rrv(i,j,3)**2)
      dd = d - rr(i,j)
      rd = dd * dd
c note: elastic energy scales here as r^4
      if (rd.gt.0.000001) espring=espring+0.5*spring(i,j)*rd*rd
      endif
  150 continue
  200 continue

      return
      end
