      subroutine deformscore(totmaxatom, maxlig, nlig,
     1 kai, tyi, rgi, iei, x, iaci, xlai, icop, we, ieins,
     2 score)
      implicit none

      integer totmaxatom, maxlig, nlig, i, ii, ijk, istart
      character*4 at, tyi,rgi
      integer iei, iaci, icop, kai, ieins
      real*8 x, xlai, we

c     return value      
      real*8 score
      
      dimension kai(totmaxatom),
     1 iei(totmaxatom), x(totmaxatom), iaci(totmaxatom),
     2 xlai(totmaxatom), icop(totmaxatom), we(totmaxatom),
     3 tyi(totmaxatom), rgi(totmaxatom), ieins(0:maxlig-1)
      
      
      score = 0
       
      end
