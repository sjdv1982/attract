      subroutine rota(x,f,delta,pm2, natom)

      implicit none

c     Parameters
      include 'max.fin'
      integer natom
      integer, parameter :: max3atom = 3*maxatom
      real*8 x,f,delta
      dimension x(max3atom),f(max3atom),delta(maxdof)
      real*8 pm2
      dimension pm2(3,3,3)

c     Local variables    
      integer i,ii,j,k,l
      real*8 torque(3,3)
c
c     calculates orientational force contributions
c     component 1: phi-angle
c     component 2: ssi-angle
c     component 3: rot-angle
c this subroutine requires the partial derivative of Cartesians with respect to the 
c Euler variables ,in case of ligand flex it is necessary to use the
c deformed ligand (in the orginal orientation, for ssi=0 etc.) 
c It depends on both the current Euler angles and the Cartesians in the orginal 
c orientation  
      
      do 5 k=1,3
      delta(k)=0.0d0
      do 4 l=1,3
      torque(k,l)=0.0d0
    4 continue
    5 continue
      
      do 10 i=1,natom
      ii=3*(i-1)
      do 20 k=1,3
      do 30 l=1,3
      torque(k,l) = torque(k,l) + x(ii+l)*f(ii+k)
   30 continue   
   20 continue
   10 continue

      do 120 j=1,3
      do 130 k=1,3
      do 140 l=1,3
      delta(j) = delta(j) + pm2(k,j,l) * torque(k,l)
  140 continue   
  130 continue
  120 continue
   
   
c     write(*,*)'rota finished'
c     write(*,*)'rotational force',(delta(j),j=1,3)
      return
      end
