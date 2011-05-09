      subroutine euler2torquemat(phi,ssi,rot,pm2)
      implicit none

c     parameters
      real*8 phi, ssi, rot,pm2
      dimension pm2(3,3,3)

c     local variables
      real*8 cp,crot,cs,cscp,cssp,sp,srot,ss,sscp,sssp
      real*8 xar,yar,xxx,yyy,zzz
c
c     calculates orientational force contributions
c     component 1: phi-angle
c     component 2: ssi-angle
c     component 3: rot-angle
c computes the partial derivative of Cartesians with respect to the 
c Euler variables 
      cs=cos(ssi)
      cp=cos(phi)
      ss=sin(ssi)
      sp=sin(phi)
      cscp=cs*cp
      cssp=cs*sp
      sscp=ss*cp
      sssp=ss*sp
      crot=cos(rot)
      srot=sin(rot)
      pm2(1,1,1)=-crot*cssp+srot*cp
      pm2(1,1,2)=-srot*cssp-crot*cp
      pm2(1,1,3)=-sssp
      pm2(2,1,1)=crot*cscp+srot*sp 
      pm2(2,1,2)=srot*cscp-crot*sp 
      pm2(2,1,3)=sscp 
      pm2(3,1,1)=0.0d0 
      pm2(3,1,2)=0.0d0 
      pm2(3,1,3)=0.0d0
      pm2(1,2,1)=-crot*sscp
      pm2(1,2,2)=-srot*sscp
      pm2(1,2,3)=cscp
      pm2(2,2,1)=-crot*sssp 
      pm2(2,2,2)=-srot*sssp 
      pm2(2,2,3)=cssp
      pm2(3,2,1)=-crot*cs 
      pm2(3,2,2)=-srot*cs 
      pm2(3,2,3)=-ss
      pm2(1,3,1)=-srot*cscp+crot*sp 
      pm2(1,3,2)=crot*cscp+srot*sp 
      pm2(1,3,3)=0.0d0 
      pm2(2,3,1)=-srot*cssp-crot*cp 
      pm2(2,3,2)=crot*cssp-srot*cp 
      pm2(2,3,3)=0.0d0
      pm2(3,3,1)=srot*ss 
      pm2(3,3,2)=-crot*ss
      pm2(3,3,3)=0.0d0 
      return
      end
