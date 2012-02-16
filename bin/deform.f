c     dligp is the part of dlig that applies to the deformations
c      for this ligand
c     ensdp is the part of ensd that applies to the ensemble differences
c      for this ligand


      subroutine deform(maxlig, max3atom,totmax3atom,
     1 maxatom, maxmode, ens, ensdp, cmorph, cmorphdp, dligp,
     2 nhm, ijk, ieins, eig,xb,x,xori,xori0, do_morph)
      implicit none

c     parameters
      integer maxlig, max3atom, totmax3atom,maxatom, maxmode
      integer ijk, ens
      real*8 cmorph
      integer nhm
      dimension nhm(0:maxlig-1)      
      integer ieins
      dimension ieins(0:maxlig-1)
      integer do_morph
      real*8 dligp, xb, eig,x,xori,xori0
      dimension dligp(maxmode),xb(totmax3atom),
     1  eig(0:maxlig-1,maxmode,max3atom), x(totmax3atom),
     2  xori(totmax3atom), xori0(totmax3atom)
      real*8 ensdp, cmorphdp
      dimension ensdp(max3atom)
      dimension cmorphdp(max3atom)

c     local variables
      integer istart,isize,i,ii,j, n, nn
      real*8 xx,yy,zz,xnull
      real*8 hmdeform
      integer hm
      dimension hmdeform(max3atom)
            
      xnull=0.0d0
      hm = 0

      istart = 0
      if (ijk.gt.0) istart = ieins(ijk-1)
      isize = ieins(ijk)-istart
c calculate deformations
      do 45 i=1,isize
      ii = 3 * (i-1)
      hmdeform(ii+1) = xnull
      hmdeform(ii+2) = xnull
      hmdeform(ii+3) = xnull
   45 continue 
      if (ens.gt.0.and.cmorph.lt.0) then
      hm = 1
      do 47 i=1,isize
      ii = 3 * (i-1)      
      hmdeform(ii+1) = hmdeform(ii+1) + ensdp(ii+1)
      hmdeform(ii+2) = hmdeform(ii+2) + ensdp(ii+2)
      hmdeform(ii+3) = hmdeform(ii+3) + ensdp(ii+3)
   47 continue
      endif
      if (do_morph.gt.0.and.cmorph.ge.0) then
      hm = 1
      do i=1,isize
      ii = 3 * (i-1)      
      hmdeform(ii+1) = hmdeform(ii+1) + ensdp(ii+1)
      hmdeform(ii+1) = hmdeform(ii+1) + cmorph*cmorphdp(ii+1)
      hmdeform(ii+2) = hmdeform(ii+2) + ensdp(ii+2)
      hmdeform(ii+2) = hmdeform(ii+2) + cmorph*cmorphdp(ii+2)
      hmdeform(ii+3) = hmdeform(ii+3) + ensdp(ii+3)
      hmdeform(ii+3) = hmdeform(ii+3) + cmorph*cmorphdp(ii+3)
      enddo
      endif
      
      do 55 n=1,nhm(ijk)
      hm = 1
      do 50 i=1,isize
      ii = 3 * (i-1)
      hmdeform(ii+1) = hmdeform(ii+1) + dligp(n) * eig(ijk,n,ii+1)
      hmdeform(ii+2) = hmdeform(ii+2) + dligp(n) * eig(ijk,n,ii+2)
      hmdeform(ii+3) = hmdeform(ii+3) + dligp(n) * eig(ijk,n,ii+3)      
   50 continue 
   55 continue

      if (hm.eq.0) then
c      x(3*istart+1:3*ieins(ijk)) = xb(3*istart+1:3*ieins(ijk))
      call memcpy(x(3*istart+1),xb(3*istart+1),
     1 (3*ieins(ijk)-3*istart)*8)
      
c      xori(3*istart+1:3*ieins(ijk)) = xori0(3*istart+1:3*ieins(ijk))
      call memcpy(xori(3*istart+1),xori0(3*istart+1),
     1 (3*ieins(ijk)-3*istart)*8)
      else
      do 60 i=istart+1,ieins(ijk)
      ii=3*(i-1)
      n = i-istart
      nn = 3 *(n-1)
      
      x(ii+1) = xb(ii+1) + hmdeform(nn+1)
      x(ii+2) = xb(ii+2) + hmdeform(nn+2)
      x(ii+3) = xb(ii+3) + hmdeform(nn+3)
      

      xori(ii+1) = xori0(ii+1) + hmdeform(nn+1)
      xori(ii+2) = xori0(ii+2) + hmdeform(nn+2)
      xori(ii+3) = xori0(ii+3) + hmdeform(nn+3)
      
   60 continue
      endif
      
      return
      end
