c     applies an in-place rotation to x 
c      x can be coordinates or forces

      subroutine rotate(maxlig, totmax3atom,
     1 rotmat,xa,ya,za, pivot,
     2 ijk, ieins, x)
      implicit none

c     parameters
      integer maxlig, totmax3atom, ijk
      real*8 xa,ya,za
      integer ieins
      dimension ieins(0:maxlig-1)
      real*8 rotmat,x,pivot
      dimension rotmat(0:8), x(totmax3atom),pivot(0:maxlig-1,3)

c     local variables
      integer istart,i,ii,j
      real*8 xx,yy,zz


      istart = 0
      if (ijk.gt.0) istart = ieins(ijk-1)
      
c makes a step in orientational coordinates 

      do 60 i=istart+1,ieins(ijk)
      ii=3*(i-1)
      
      xx = x(ii+1)
      yy = x(ii+2)
      zz = x(ii+3)
      x(ii+1) = xx * rotmat(0) + yy * rotmat(1) + zz * rotmat(2)
      x(ii+2) = xx * rotmat(3) + yy * rotmat(4) + zz * rotmat(5)
      x(ii+3) = xx * rotmat(6) + yy * rotmat(7) + zz * rotmat(8)
      
      x(ii+1) = x(ii+1) + xa + pivot(ijk,1)
      x(ii+2) = x(ii+2) + ya + pivot(ijk,2)
      x(ii+3) = x(ii+3) + za + pivot(ijk,3)
   60 continue
      return
      end

c     applies an in-place rotation to selected ligand x
c      x can be coordinates or forces

      subroutine rotate1(max3atom,
     1 rotmat,xa,ya,za,pivot,natom,x)

      implicit none

c     parameters
      integer max3atom, natom
      real*8 xa,ya,za
      real*8 rotmat,x,pivot
      dimension rotmat(0:8), x(max3atom),pivot(3)

c     local variables
      integer i,ii,j
      real*8 xx,yy,zz

c makes a step in orientational coordinates 

      do 60 i=1,natom
      ii=3*(i-1)
      
      xx = x(ii+1)
      yy = x(ii+2)
      zz = x(ii+3)
      x(ii+1) = xx * rotmat(0) + yy * rotmat(1) + zz * rotmat(2)
      x(ii+2) = xx * rotmat(3) + yy * rotmat(4) + zz * rotmat(5)
      x(ii+3) = xx * rotmat(6) + yy * rotmat(7) + zz * rotmat(8)
      
      x(ii+1) = x(ii+1) + xa + pivot(1)
      x(ii+2) = x(ii+2) + ya + pivot(2)
      x(ii+3) = x(ii+3) + za + pivot(3)
   60 continue
      return
      end
