      subroutine write_pdb(totmaxatom, maxlig, nlig,
     1 kai, tyi, rgi, iei, x, iaci, xlai, icop, we, ieins)
      implicit none

      integer totmaxatom, maxlig, nlig, i, ii, ijk, istart
      character*4 at, tyi,rgi
      integer iei, iaci, icop, kai, ieins
      real*8 x, xlai, we
      
      dimension kai(totmaxatom),
     1 iei(totmaxatom), x(totmaxatom), iaci(totmaxatom),
     2 xlai(totmaxatom), icop(totmaxatom), we(totmaxatom),
     3 tyi(totmaxatom), rgi(totmaxatom), ieins(0:maxlig-1)
       
      at = 'ATOM'
      do 80 ijk=0,nlig-1
      istart = 0
      if (ijk.gt.0) istart=ieins(ijk-1)
      do 60 i=istart, ieins(ijk)-1
      ii = 3 * i
      write(*,26) at,kai(i+1),tyi(i+1),rgi(i+1),iei(i+1),x(ii+1),
     1  x(ii+2),x(ii+3),iaci(i+1),xlai(i+1),icop(i+1),we(i+1)
   60 continue  
   80 write(*,'(a3)') 'TER'
      continue
      write(*,'(a3)') 'END'
   26 format(a4,i7,2x,a4,a4,2x,i3,4x,f8.3,f8.3,f8.3,i5,f8.3,i2,f5.2)      
      end
