c     applies an in-place force scale and rotation to 
c      the forces f of a selected ligand

      subroutine forcerotscale(max3atom,
     1 forcescale,rotmat,natom,f)

      implicit none

c     parameters
      integer max3atom, natom
      real*8 forcescale,rotmat,f
      dimension rotmat(0:8), f(max3atom)

c     local variables
      integer i,ii,j
      real*8 xx,yy,zz
      real*8 minv
      dimension minv(0:8)
      real*8 scalemin
      scalemin = 0.000001

      if (forcescale.le.scalemin) then 
      return
      endif
      
      minv(:) = rotmat(:)
      call matinv(minv)
      
      do 60 i=1,natom
      ii=3*(i-1)
      
      xx = f(ii+1) / forcescale
      yy = f(ii+2) / forcescale
      zz = f(ii+3) / forcescale
      f(ii+1) = xx * minv(0) + yy * minv(1) + zz * minv(2)
      f(ii+2) = xx * minv(3) + yy * minv(4) + zz * minv(5)
      f(ii+3) = xx * minv(6) + yy * minv(7) + zz * minv(8)
            
   60 continue
      
      return
      end
