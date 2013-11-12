c     dligp is the part of dlig that applies to the deformations
c      for this ligand
c     ensdp is the part of ensd that applies to the ensemble differences
c      for this ligand


      subroutine deform(ens, ensdp, cmorph, cmorphdp, dligp,
     1 nhm,nihm, ijk, ieins, eig, index_eig, index_val,
     2 xb,x,xori,xori0, do_morph)
      implicit none

c     parameters
      include 'max.fin'
      integer,parameter :: max3atom = 3*maxatom
      integer,parameter :: totmax3atom = 3*totmaxatom
      integer ijk, ens
      real*8 cmorph
      integer nhm, nihm
      dimension nhm(0:maxlig-1)
      dimension nihm(0:maxlig-1)
      integer ieins
      dimension ieins(0:maxlig-1)
      integer do_morph
      integer index_eig
      dimension index_eig(0:maxlig-1,maxindexmode,maxlenindexmode)
      real*8 dligp, xb, eig,index_val,x,xori,xori0
      dimension dligp(maxmode+maxindexmode),xb(totmax3atom),
     1  eig(0:maxlig-1,maxmode,max3atom), x(totmax3atom),
     2  xori(totmax3atom), xori0(totmax3atom)
      dimension index_val(0:maxlig-1,maxindexmode,maxlenindexmode)
      real*8 ensdp, cmorphdp
      dimension ensdp(max3atom)
      dimension cmorphdp(max3atom)

c     local variables
      integer istart,isize,i,ii,j, n, nn
      real*8 xx,yy,zz,xnull
      real*8 hmdeform
      integer hm
      dimension hmdeform(max3atom)
      integer, parameter :: ERROR_UNIT = 0
            
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
c      if(eig(ijk,n,ii+1).ne.0) then
c      write(ERROR_UNIT,*) "Contribution",ii+1, eig(ijk,n,ii+1),dligp(n)
c      else if (eig(ijk,n,ii+2).ne.0) then
c      write(ERROR_UNIT,*) "Contribution",ii+2, eig(ijk,n,ii+2),dligp(n)
c      else if (eig(ijk,n,ii+3).ne.0) then
c      write(ERROR_UNIT,*) "Contribution",ii+3, eig(ijk,n,ii+3),dligp(n)
c      endif
      hmdeform(ii+1) = hmdeform(ii+1) + dligp(n) * eig(ijk,n,ii+1)
      hmdeform(ii+2) = hmdeform(ii+2) + dligp(n) * eig(ijk,n,ii+2)
      hmdeform(ii+3) = hmdeform(ii+3) + dligp(n) * eig(ijk,n,ii+3)
   50 continue 
   55 continue
c      write(ERROR_UNIT,*) "Index modes:", nihm(ijk)
      do 66 n=1,nihm(ijk)
      hm = 1
      do 65 i=1,maxlenindexmode
      ii = index_eig(ijk,n,i)
      if (ii.ne.-1) then
c      write(ERROR_UNIT,*) ii+1,index_val(ijk,n,i),dligp(6+nhm(ijk)+n)
      hmdeform(ii+1)=hmdeform(ii+1)+
     1 dligp(nhm(ijk)+n)*index_val(ijk,n,i)
      else
      exit
      endif
   65 continue
   66 continue

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
