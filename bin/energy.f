      subroutine energy(maxdof,maxmolpair,
     1  maxlig, maxatom,totmaxatom,maxmode,maxres,
     2 cartstatehandle, ministatehandle,
     3 iab,iori,itra,ieig,fixre,gridmode,
     4 phi,ssi,rot,xa,ya,za,dlig,seed,
     5 e,energies,delta)
     
      implicit none
      
c     Parameters
      integer cartstatehandle,ministatehandle
      integer maxlig,maxdof,maxatom,totmaxatom,maxmode,maxmolpair,
     1 maxres
      integer iori,itra,ieig,fixre,gridmode,seed
      real*8 energies, delta
      dimension energies(6), delta(maxdof)
      real*8 e, epair

      real*8 phi, ssi, rot, dlig, xa, ya, za
      dimension phi(maxlig), ssi(maxlig), rot(maxlig)
      dimension dlig(maxmode, maxlig)
      dimension xa(maxlig), ya(maxlig), za(maxlig)

c     Local variables      
      integer i, k, n, idr, idl,ii
      integer gridptr_dmmy
      pointer(gridptr, gridptr_dmmy)      
      integer nlig, molpairs, molpairhandle,molpairhandles
      dimension molpairhandles(maxmolpair)
      integer nhm
      dimension nhm(maxlig)
      pointer(ptr_nhm,nhm)
      real*8 pairenergies,deltar, deltal
      dimension pairenergies(6),deltar(maxdof),deltal(maxdof)
      integer jb,ju,jl,nmodes,iab,fixre2

c  return values:
c  energies is an array of double with a value for every energy type
c  (vdw, elec, ...)
c  delta contains all the DOF gradients

      call ministate_get_molpairhandles(
     1 ministatehandle, molpairhandles, molpairs)

      call cartstate_get_nlig_nhm(cartstatehandle,
     1 nlig,ptr_nhm)

c
c  all variables without lig-hm
c
      jb=3*iori*(nlig-fixre)+3*itra*(nlig-fixre)
c  all variables including lig-hm
      nmodes = 0
      do 5 i=1, nlig
      nmodes = nmodes + nhm(i)
5     continue
      ju=jb+ieig*nmodes
c  only trans or ori
      jl=3*iori*(nlig-fixre)

      do 10 i=1,6
      energies(i) = 0
  10  continue  
      
      do 11 i=1,maxdof
      delta(i) = 0
  11  continue  
             
c  iterate over all pairs: call pairenergy...
      do 30 k=1,molpairs
      
      molpairhandle = molpairhandles(k)
      call molpair_get_rl(molpairhandle,idr,idl,gridptr)

      if (gridmode.eq.1.and.gridptr.ne.0) then
      fixre2 = 1
      else
      fixre2 = 0
      if (idr.eq.0) fixre2 = fixre
      endif
      
      call pairenergy(maxlig,maxatom,totmaxatom,maxmode,
     1 maxdof,maxmolpair,maxres,
     2 cartstatehandle,molpairhandle,
     3 iab,fixre2,gridptr,
     4 phi(idr+1),ssi(idr+1),rot(idr+1),
     5 xa(idr+1),ya(idr+1),za(idr+1),dlig(:,idr+1),
     6 phi(idl+1),ssi(idl+1),rot(idl+1),
     7 xa(idl+1),ya(idl+1),za(idl+1),dlig(:,idl+1),     
     8 pairenergies,deltar,deltal)
     
c  ...and sum up the energies and deltas            
      if ((iori.eq.1).AND.(fixre2.eq.0)) then
      ii = 3 * (idr-fixre)
      delta(ii+1) = delta(ii+1) + deltar(1)
      delta(ii+2) = delta(ii+2) + deltar(2)
      delta(ii+3) = delta(ii+3) + deltar(3)
      endif
      if ((itra.eq.1).AND.(fixre2.eq.0)) then
      ii = jl + 3 * (idr-fixre)
      delta(ii+1) = delta(ii+1) + deltar(4)
      delta(ii+2) = delta(ii+2) + deltar(5)
      delta(ii+3) = delta(ii+3) + deltar(6)
      endif
      
      if ((iori.eq.1).AND.((idl+1).gt.fixre)) then
      ii = 3 * (idl-fixre)
      delta(ii+1) = delta(ii+1) + deltal(1)
      delta(ii+2) = delta(ii+2) + deltal(2)
      delta(ii+3) = delta(ii+3) + deltal(3)
      endif
      if ((itra.eq.1).AND.((idl+1).gt.fixre)) then
      ii = jl + 3 * (idl-fixre)
      delta(ii+1) = delta(ii+1) + deltal(4)
      delta(ii+2) = delta(ii+2) + deltal(5)
      delta(ii+3) = delta(ii+3) + deltal(6)
      endif

      
      if (ieig.eq.1) then
      ii = jb
      do 13 n=1,idr
      ii = ii + nhm(n)
13    continue  
      do 14 n=1,nhm(idr+1)  
      delta(ii+n) = delta(ii+n) + deltar(6+n)
14    continue

      ii = jb
      do 15 n=1,idl
      ii = ii + nhm(n)
15    continue  
      do 16 n=1,nhm(idl+1)  
      delta(ii+n) = delta(ii+n) + deltal(6+n)
16    continue

      endif
      
      epair = 0
      if (idr.lt.idl.OR.gridmode.eq.2.OR.(fixre.eq.1.and.idl.eq.0)) then
c     In case of grids, the receptor forces must be calculated using
c      an additional molpair, with receptor and ligands swapped
c      (unless gridmode is 2)
c      in that case, we are only interested in the forces, not the energies
      do 20 i=1,6
      energies(i) = energies(i) + pairenergies(i)
      epair = epair + pairenergies(i)
20    continue  
      endif
      
30    continue
   
      call globalenergy(
     1  maxlig,maxatom,totmaxatom,maxmode,maxdof,
     2	cartstatehandle, ministatehandle,
     3  phi,ssi,rot,xa,ya,za,dlig,seed,
     4  iab,iori,itra,ieig,fixre, energies, delta)
      
      e = 0
      do 990,i=1,6
      e = e + energies(i)
990   continue  
c      write(*,*), 'ENERGY', e,
c     1  delta(1),delta(2),delta(3),delta(4),delta(5),delta(6),
c     2  delta(7),delta(8),delta(9),delta(10),delta(11),delta(12),
c     3  delta(13)

c     e = nonbonded+potential+erest+esolv+enlig+eem 
      return
      end


