      subroutine energy(cartstatehandle, ministatehandle,
     1 iab,iori,itra,ieig,iindex,fixre,gridmode,
     2 ens,phi,ssi,rot,xa,ya,za,morph,dlig,
     3 locrests, has_locrests, seed,
     4 e,energies,delta,deltamorph)
     
      implicit none
      
c     Parameters
      integer cartstatehandle,ministatehandle
      include 'max.fin'
      integer iori,itra,ieig,iindex,fixre,gridmode,seed
      real*8 energies, delta, deltamorph
      dimension energies(6), delta(maxdof),deltamorph(maxlig)
      real*8 deltamorphr, deltamorphl      
      real*8 e
      
      real*8 locrests
      dimension locrests(3,maxlig)
      integer has_locrests      
      dimension has_locrests(maxlig)      

      integer ens
      dimension ens(maxlig)
      real*8 phi, ssi, rot, dlig, xa, ya, za, morph
      dimension phi(maxlig), ssi(maxlig), rot(maxlig)
      dimension dlig(maxmode+maxindexmode, maxlig)
      dimension xa(maxlig), ya(maxlig), za(maxlig)
      dimension morph(maxlig)

c     Local variables      
      real*8 epair, xnull
      integer i, k, n, idr, idl,ii
      integer gridptr_dmmy
      pointer(gridptr, gridptr_dmmy)      
      integer nlig, molpairs, molpairhandle,molpairhandles
      integer use_energy
      dimension molpairhandles(maxmolpair)
      integer nhm
      dimension nhm(maxlig)
      integer nihm
      dimension nihm(maxlig)
      pointer(ptr_nhm,nhm)
      pointer(ptr_nihm,nihm)
      real*8 pairenergies,deltar, deltal
      dimension pairenergies(6),deltar(6+maxmode+maxindexmode)
      dimension deltal(6+maxmode+maxindexmode)
      integer jb,jl,ju,iab,fixre2, ghost

      real*8 ff, frot
      dimension ff(maxlig)
      dimension frot(9,maxlig)      
      pointer(ptr_ff, ff)
      pointer(ptr_frot, frot)

      real*8 delta0
      dimension delta0(maxdof)
      integer nmodes,nimodes
      integer jn
      integer, parameter :: ERROR_UNIT = 0
c  return values:
c  energies is an array of double with a value for every energy type
c  (vdw, elec, ...)
c  delta0 contains all the preliminary DOF gradients
c  it will be folded into delta by axsym_fold_grads

      xnull = 0.0d0
      e = xnull

      call apply_axsym(cartstatehandle, morph, ens, 
     1 phi, ssi, rot, xa, ya, za, dlig, locrests)
           
      call ministate_ghost(ministatehandle, ghost)
      
      call ministate_get_molpairhandles(
     1 ministatehandle, molpairhandles, molpairs)

      call cartstate_get_nlig_nhm(cartstatehandle,
     1 nlig,ptr_nhm,ptr_nihm)
 
      call cartstate_get_forcerot(cartstatehandle,ptr_ff,ptr_frot)

c
c  all variables without lig-hm
c
      jb=3*iori*(nlig-fixre)+3*itra*(nlig-fixre)
c  only trans or ori
      jl=3*iori*(nlig-fixre)
      nmodes = 0
      nimodes = 0
      do i=1,nlig
      nmodes = nmodes + nhm(i)
      nimodes = nimodes + nihm(i)
      enddo
      ju = jb + ieig*nmodes
      jn = ju + iindex+nimodes
      do i=1,6
      energies(i) = xnull
      enddo 
      
      do i=1,maxdof
      delta0(i) = xnull
      delta(i) = xnull
   1  enddo
             
      do i=1,maxlig
      deltamorph(i)= xnull
      enddo       
                  
c  iterate over all pairs: call pairenergy...
      do 30 k=1,molpairs
      
      molpairhandle = molpairhandles(k)
      call molpair_get_rl(molpairhandle,idr,idl,gridptr,use_energy)

      if (gridmode.eq.1.and.gridptr.ne.0) then
      fixre2 = 1
      else
      fixre2 = 0
      if (idr.eq.0) fixre2 = fixre
      endif
      
      if (ghost.eq.0) then
      call pairenergy(cartstatehandle,molpairhandle,
     1 iab,fixre2,gridptr,
     2 ff(idr+1),frot(:,idr+1),ff(idl+1),frot(:,idl+1),
     3 ens(idr+1),phi(idr+1),ssi(idr+1),rot(idr+1),
     4 xa(idr+1),ya(idr+1),za(idr+1),morph(idr+1),dlig(:,idr+1),
     5 ens(idl+1),phi(idl+1),ssi(idl+1),rot(idl+1),
     6 xa(idl+1),ya(idl+1),za(idl+1),morph(idl+1),dlig(:,idl+1),
     7 pairenergies, deltar, deltal, deltamorphr, deltamorphl)
      write(ERROR_UNIT,*)'deltar',(deltar(i),i=1,jn)
      write(ERROR_UNIT,*)'deltal',(deltal(i),i=1,jn)
c  ...and sum up the energies and deltas            
      if ((iori.eq.1).AND.(fixre2.eq.0)) then
      ii = 3 * (idr-fixre)
      delta0(ii+1) = delta0(ii+1) + deltar(1)
      delta0(ii+2) = delta0(ii+2) + deltar(2)
      delta0(ii+3) = delta0(ii+3) + deltar(3)
      endif
      if ((itra.eq.1).AND.(fixre2.eq.0)) then
      ii = jl + 3 * (idr-fixre)
      delta0(ii+1) = delta0(ii+1) + deltar(4)
      delta0(ii+2) = delta0(ii+2) + deltar(5)
      delta0(ii+3) = delta0(ii+3) + deltar(6)
      endif
      
      if ((iori.eq.1).AND.((idl+1).gt.fixre)) then
      ii = 3 * (idl-fixre)
      delta0(ii+1) = delta0(ii+1) + deltal(1)
      delta0(ii+2) = delta0(ii+2) + deltal(2)
      delta0(ii+3) = delta0(ii+3) + deltal(3)
      endif
      
      if ((itra.eq.1).AND.((idl+1).gt.fixre)) then
      ii = jl + 3 * (idl-fixre)
      delta0(ii+1) = delta0(ii+1) + deltal(4)
      delta0(ii+2) = delta0(ii+2) + deltal(5)
      delta0(ii+3) = delta0(ii+3) + deltal(6)
      endif
      
      if (ieig.eq.1.OR.iindex.eq.1) then
      ii = jb
      do 13 n=1,idr
      ii = ii + nhm(n)
13    continue  
      do 14 n=1,nhm(idr+1)  
      delta0(ii+n) = delta0(ii+n) + deltar(6+n)
14    continue

      ii = ju
      do 37 n=1,idr
      ii = ii + nihm(n)
37    continue
      do 17 n=1,nihm(idr+1)
      delta0(ii+n) = delta0(ii+n) + deltar(6+nhm(idr+1)+n)
17    continue

      ii = jb
      do 15 n=1,idl
      ii = ii + nhm(n)
15    continue  
      do 16 n=1,nhm(idl+1)  
      delta0(ii+n) = delta0(ii+n) + deltal(6+n)
16    continue

      ii = ju
      do 38 n=1,idl
      ii = ii + nihm(n)
38    continue
      do 19 n=1,nihm(idl+1)
      delta0(ii+n) = delta0(ii+n) + deltal(6+nhm(idl+1)+n)
19    continue
      endif

      deltamorph(idl+1) = deltamorph(idl+1) + deltamorphl
      deltamorph(idr+1) = deltamorph(idr+1) + deltamorphr
      
      epair = 0
      if (use_energy.gt.0) then
      do 20 i=1,6
c     write(ERROR_UNIT, *) "Pair energy:", pairenergies(i)
      energies(i) = energies(i) + pairenergies(i)
      epair = epair + pairenergies(i)
20    continue  
      endif

      endif 
c     endif ghost.eq.0
       write(ERROR_UNIT,*)'delta0',(delta0(i),i=1,jn)
30    continue
      call globalenergy(
     1	cartstatehandle, ministatehandle,
     2  ens,phi,ssi,rot,xa,ya,za,morph,dlig,
     3  locrests, has_locrests, seed,
     4  iab,iori,itra,ieig,iindex,fixre,
     5  energies, delta0, deltamorph)
      write(ERROR_UNIT,*)'delta0',(delta0(i),i=1,jn)
      e = 0
      do 990,i=1,6
c      write(ERROR_UNIT, *) "Global energy:", energies(i)
      e = e + energies(i)
990   continue  
      call axsym_fold_grads(ministatehandle, cartstatehandle, 
     1 delta0, delta, deltamorph)
      write(ERROR_UNIT,*)'delta',(delta(i),i=1,jn)
c      write(ERROR_UNIT,*), 'Final ENERGY', e
c     1  delta(1),delta(2),delta(3),delta(4),delta(5),delta(6),
c     2  delta(7),delta(8),delta(9),delta(10),delta(11),delta(12),
c     3  delta(13)

c     e = nonbonded+potential+erest+esolv+enlig+eem 

      return
      end


