c     first rotate using rot
c       this is a rotation around the Z axis
c         rotating Y into X and X into -Y
c     then rotate using ssi
c       this is a rotation around the (new) Y axis
c         rotating X into -Z and Z into X
c     finally, rotate using phi
c       this is a rotation around the (new) Z axis
c        rotating X into Y and Y into -X
c
c     so in principle, the reverse of rot,ssi,phi
c      should be phi,-ssi,rot
c
      subroutine euler2rotmat(phi,ssi,rot,rotmat)
      implicit none

c     parameters
      real*8 phi, ssi, rot, rotmat
      dimension rotmat(0:9)

      real*8 cs,cp,ss,sp,cscp,cssp,sscp,sssp,crot,srot,
     1 xar,yar

      cs=dcos(ssi)
      cp=dcos(phi)
      ss=dsin(ssi)
      sp=dsin(phi)
      cscp=cs*cp
      cssp=cs*sp
      sscp=ss*cp
      sssp=ss*sp
      crot=dcos(rot)
      srot=dsin(rot)
      
c      xar=xb(ii+1)*crot+xb(ii+2)*srot
c      yar=-xb(ii+1)*srot+xb(ii+2)*crot       
c      x(ii+1)=xar*cscp-yar*sp+xb(ii+3)*sscp
c      x(ii+2)=xar*cssp+yar*cp+xb(ii+3)*sssp
c      x(ii+3)=-xar*ss+xb(ii+3)*cs
      
      rotmat(0) = crot * cscp + srot * sp  
      rotmat(1) = srot * cscp - crot * sp
      rotmat(2) = sscp
      
      rotmat(3) = crot * cssp - srot * cp
      rotmat(4) = srot * cssp + crot * cp
      rotmat(5) = sssp
      
      rotmat(6) = -crot * ss
      rotmat(7) = -srot * ss
      rotmat(8) = cs
c      write(*,*), rotmat(0), rotmat(1), rotmat(2)
c      write(*,*), rotmat(3), rotmat(4), rotmat(5)
c      write(*,*), rotmat(6), rotmat(7), rotmat(8)

      return
      end
