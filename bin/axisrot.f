      subroutine axisrot(rr,rotmat)
      implicit none

c     parameters
      real*8 rotmat, rr
      dimension rotmat(0:9),rr(4)

      real*8 cs,cm,ss,r1r2,r1r3,r2r3

      cs=dcos(rr(4))
      ss=dsin(rr(4))
      cm=1.d0-cs
      r1r2=rr(1)*rr(2)
      r1r3=rr(1)*rr(3)
      r2r3=rr(2)*rr(3)

      rotmat(0) =  cs+rr(1)*rr(1)*cm     
      rotmat(1) = r1r2*cm-rr(3)*ss
      rotmat(2) = r1r3*cm+rr(2)*ss

      rotmat(3) = r1r2*cm+rr(3)*ss
      rotmat(4) = cs+rr(2)*rr(2)*cm
      rotmat(5) =  r2r3*cm-rr(1)*ss

      rotmat(6) = r1r3*cm-rr(2)*ss
      rotmat(7) = r2r3*cm+rr(1)*ss
      rotmat(8) = cs+rr(3)*rr(3)*cm

      return
      end
