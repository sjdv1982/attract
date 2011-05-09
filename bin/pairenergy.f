       subroutine pairenergy(maxlig,maxatom,totmaxatom,maxmode,
     1 maxdof,maxmolpair,maxres,
     2 cartstatehandle,molpairhandle,iab,fixre,gridptr,
     3 phir, ssir, rotr, xar, yar, zar, dligr,
     4 phil, ssil, rotl, xal, yal, zal, dligl,
     5 energies, deltar, deltal)

       implicit none

c      Parameters
       integer maxlig,maxatom,totmaxatom,maxmode,maxmolpair,maxdof,
     1  maxres
       integer cartstatehandle,molpairhandle
       integer iab, fixre
       integer gridptr_dmmy
       pointer(gridptr, gridptr_dmmy)
       real*8 phir, ssir, rotr, xar, yar, zar, dligr
       real*8 phil, ssil, rotl, xal, yal, zal, dligl
       dimension dligr(maxdof), dligl(maxdof)
       real*8 energies
       dimension energies(6)
       real*8 deltar(6+maxmode),deltal(6+maxmode)
      
c      Handle variables: molpair
       integer idr, idl
       integer iactr,iactl,nonr,nonl
       dimension iactr(maxatom),iactl(maxatom)
       dimension nonr(maxmolpair),nonl(maxmolpair)
       integer nonp
       pointer(ptr_iactr,iactr)
       pointer(ptr_iactl,iactl)
       pointer(ptr_nonr,nonr)
       pointer(ptr_nonl,nonl)

c      Handle variables: cartstate       
       real*8 xr,wer,chair,fr,pivotr
       dimension xr(maxatom),wer(maxatom),chair(maxatom)
       dimension fr(maxatom), pivotr(3)
       integer nmaxcor,ieir,iacir,natcor,icopr
       dimension nmaxcor(maxres),ieir(maxatom),
     1 iacir(maxatom),natcor(maxres),icopr(maxatom)       
       integer nresr, natomr
       pointer(ptr_ieir,ieir)
       pointer(ptr_xr,xr)
       pointer(ptr_fr,fr)
       pointer(ptr_iacir,iacir)
       pointer(ptr_icopr,icopr)
       pointer(ptr_wer,wer)
       pointer(ptr_chair,chair)
       pointer(ptr_nmaxcor,nmaxcor)
       pointer(ptr_natcor,natcor)

       real*8 xl,wel,chail,fl,pivotl
       dimension xl(maxatom),wel(maxatom),chail(maxatom)
       dimension fl(maxatom), pivotl(3)
       integer nmaxcol,ieil,iacil,natcol,icopl
       dimension nmaxcol(maxres),ieil(maxatom),
     1  iacil(maxatom),natcol(maxres),icopl(maxatom)
       integer nresl, natoml
       pointer(ptr_ieil,ieil)
       pointer(ptr_xl,xl)
       pointer(ptr_fl,fl)
       pointer(ptr_iacil,iacil)
       pointer(ptr_icopl,icopl)
       pointer(ptr_wel,wel)
       pointer(ptr_chail,chail)
       pointer(ptr_nmaxcol,nmaxcol)
       pointer(ptr_natcol,natcol)
       
       integer ncopr, ncopl
       dimension ncopr(0:10,0:20,maxres), ncopl(0:10,0:20,maxres)
       pointer(ptr_ncopr,ncopr)
       pointer(ptr_ncopl,ncopl)

c      Handle variables: forcefield parameters
       integer potshape
       real*8 rbc,rc,ac,emin,rmin2,vlj
       dimension rbc(99,99),rc(99,99),ac(99,99),emin(99,99),
     1  rmin2(99,99)
       integer ipon
       dimension ipon(99,99)
       pointer(ptr_rbc,rbc)
       pointer(ptr_rc,rc)
       pointer(ptr_ac,ac)
       pointer(ptr_emin,emin)
       pointer(ptr_rmin2,rmin2)
       pointer(ptr_ipon,ipon)

c      Handle variables: full coordinates and modes
       real*8 xb,x,xori,xori0,eig
       integer ieins,nhm
       dimension xb(3*totmaxatom),x(3*totmaxatom)
       dimension xori(3*totmaxatom),xori0(3*totmaxatom)
       dimension eig(maxlig,maxmode,3*maxatom)
       dimension ieins(maxlig),nhm(maxlig)
       pointer(ptr_xb,xb)
       pointer(ptr_x,x)
       pointer(ptr_xori,xori)
       pointer(ptr_xori0,xori0)       
       pointer(ptr_eig,eig)
       pointer(ptr_ieins,ieins)
       pointer(ptr_nhm,nhm)

c      Local variables: matrices
       real*8 rotmatr,rotmatl,rotmatd,rotmatd2,rotmatrinv,rotmatlinv
       dimension rotmatr(9),rotmatl(9),rotmatd(9),rotmatrinv(9),
     1  rotmatlinv(9),rotmatd2(9)
       real*8 tr,tl,td,td0,xnull
       dimension tr(3),tl(3),td(3),td0(3)
       real*8 pivotd,pivotd0,pivotnull
       dimension pivotd(3),pivotd0(3),pivotnull(3)
       real*8 pm2(3,3,3), pm20(3,3,3), mat(9)

c  energies is an array of double with a value for every energy type
c  (vdw, elec, ...) for this pair
c
c phir, ssir, ssil etc. are scalars
c dligr / dligl is an array of nhm[ligand]/nhm[receptor]
c deltar/deltar is an array of 6 + nhm[ligand]/nhm[receptor]

c     Local variables: other      
      integer i, j,l,n, use_grid, rigid, true, false
      real*8 enon,epote,zero,e
      real*8 flcopy,xl0
      dimension flcopy(maxatom),xl0(3*maxatom)
      real*8 deltar0(maxdof)
      
      zero=0.0d0
      true = 1
      false = 0
      pivotnull(1) = zero
      pivotnull(2) = zero
      pivotnull(3) = zero

c reset forces

      call reset_forces(cartstatehandle)

c Are we using a grid?
       use_grid = 0
       if (gridptr.ne.0) then
c TODO?: use dligr to determine RMSD
c  if RMSD is zero: rigid mode? or keep non-rigid?
c  if RMSD is larger: non-rigid mode
       rigid = 0
       use_grid = 1
       endif

c get molpair 

       call molpair_get_values(molpairhandle,idr,idl,
     1  ptr_iactr,ptr_iactl,nonp,ptr_nonr,ptr_nonl)

c get full coordinates and normal modes

       call cartstate_f_pairenergy(cartstatehandle,ptr_nhm,ptr_ieins,
     1  ptr_eig, ptr_xb, ptr_x, ptr_xori, ptr_xori0)

c apply normal mode deformations to receptor and ligand
       call deform(maxlig,3*maxatom,3*totmaxatom,maxatom,maxmode,
     1  dligr,nhm,idr,ieins,eig,xb,x,xori,xori0)
       call deform(maxlig,3*maxatom,3*totmaxatom,maxatom,maxmode,
     1  dligl,nhm,idl,ieins,eig,xb,x,xori,xori0)

c select deformed coordinates: select non-pivotized coordinates for receptor

       call cartstate_select_ligand(cartstatehandle,idr,
     1  natomr,nresr,ptr_ieir,ptr_xr,ptr_fr,pivotr,ptr_iacir,
     2  ptr_icopr,ptr_wer,ptr_chair,ptr_ncopr,
     3  ptr_nmaxcor,ptr_natcor, true)

       call cartstate_select_ligand(cartstatehandle,idl,
     1  natoml,nresl,ptr_ieil,ptr_xl,ptr_fl,pivotl,ptr_iacil,
     2  ptr_icopl,ptr_wel,ptr_chail,ptr_ncopl,
     3  ptr_nmaxcol,ptr_natcol, false)
       xl0(:) = xl(:)
     
       call cartstate_get_parameters(cartstatehandle,
     1  ptr_rbc,ptr_rc,ptr_ac,ptr_emin,ptr_rmin2,ptr_ipon,potshape)
            
       xnull = 0.0d0   
       do 10 i=1,6+maxmode
       deltal(i) = xnull
       deltar(i) = xnull 
       deltar0(i)= xnull 
10     continue  

c rotate ligand into receptor frame 
       call euler2rotmat(phir,ssir,rotr,rotmatr)
       call euler2rotmat(phil,ssil,rotl,rotmatl)
       	      
       call matcopy(rotmatr,rotmatrinv)

       tr(1) = xar
       tr(2) = yar
       tr(3) = zar
       call matinv(rotmatrinv)
       
       call matcopy(rotmatl,rotmatlinv)
       tl(1) = xal
       tl(2) = yal
       tl(3) = zal
       call matinv(rotmatlinv)
       
       call matmult(rotmatl, rotmatrinv,rotmatd)
       td0(1) = tl(1)-tr(1)
       td0(2) = tl(2)-tr(2)
       td0(3) = tl(3)-tr(3)       
       call vecmatmult(td0,rotmatrinv,td)
       
       pivotd0(1) = pivotl(1)-pivotr(1)
       pivotd0(2) = pivotl(2)-pivotr(2)
       pivotd0(3) = pivotl(3)-pivotr(3)
       call vecmatmult(pivotd0,rotmatrinv,pivotd)
       pivotd(1) = pivotd(1) + pivotr(1)
       pivotd(2) = pivotd(2) + pivotr(2)
       pivotd(3) = pivotd(3) + pivotr(3)

c calculate nonbonded energies and forces
       
       call rotate1(3*maxatom,
     1  rotmatd,td(1),td(2),td(3),
     2	 pivotd, natoml,xl)

c       write(*,'(a2,3f8.3)'), 'P',xr(3*902+1),xr(3*902+2),xr(3*902+3)
c       write(*,'(a2,3f8.3)'), 'P',xl(3*317+1),xl(3*317+2),xl(3*317+3)
c       write(*,*),(xr(3*902+1)-xl(3*317+1))*(xr(3*902+1)-xl(3*317+1))
c     1  +(xr(3*902+2)-xl(3*317+2))*(xr(3*902+2)-xl(3*317+2))+
c     2  (xr(3*902+3)-xl(3*317+3))*(xr(3*902+3)-xl(3*317+3))
       
       if (use_grid.eq.1) then
       if (fixre.eq.0) then
       call euler2torquemat(phir,ssir,rotr,pm2)
c      our pm2 matrix is correct for global frame forces; however, our forces
c       will be in the receptor frame
c      therefore, we must rotate the pm2 matrix
       pm20 = pm2
       do 110 j=1,3
       do 120 l=1,3
       mat = rotmatrinv
       pm2(1,j,l) = mat(1)*pm20(1,j,l)+mat(2)*pm20(2,j,l)
     1  + mat(3)*pm20(3,j,l)
       pm2(2,j,l) = mat(4)*pm20(1,j,l)+mat(5)*pm20(2,j,l)
     1  + mat(6)*pm20(3,j,l)
       pm2(3,j,l) = mat(7)*pm20(1,j,l)+mat(8)*pm20(2,j,l)
     1  + mat(9)*pm20(3,j,l)
120    continue
110    continue       
       endif
       
       call nonbon_grid(gridptr,rigid,iab,fixre,xl,xr,pivotr,tr,
     1  wel,wer,chail,chair,iacil,iacir,natoml,natomr,
     2  rc,ac,emin,rmin2,ipon,potshape,fl,enon,epote,
     3  fr,pm2,deltar0) 
c      rotate delta-translate back into global frame
       call rotate1(3*maxatom,
     1  rotmatr,zero,zero,zero,
     2	 pivotnull, 1,deltar0(4))     
       else  
       call molpair_pairgen(molpairhandle,cartstatehandle,nonp)
       call nonbon8(maxatom,maxmolpair,
     1  iab,xl,xr,fl,fr,wel,wer,chair,chail,ac,rc,
     2  emin,rmin2,iacir,iacil,nonr,nonl,ipon,nonp,
     3  potshape,enon,epote)
       endif
       
       energies(1) = enon
       energies(2) = epote
       
c calculate receptor mode deltas
       call ligmin(maxlig,maxdof,maxmode,maxatom,
     1  fr,natomr,idr,eig,nhm(idr+1),deltar) 
       
c rotate ligand forces into ligand frame
       call matmult(rotmatrinv,rotmatl,rotmatd2)
       flcopy(:) = fl(:)
       call rotate1(3*maxatom,
     1  rotmatd2,zero,zero,zero,
     2	 pivotnull, natoml,flcopy)

c calculate ligand mode deltas
       call ligmin(maxlig,maxdof,maxmode,maxatom,
     1  flcopy,natoml,idl,eig,nhm(idl+1),deltal) 

c rotate all forces into global frame
       call rotate1(3*maxatom,
     1  rotmatr,zero,zero,zero,
     2	 pivotnull, natomr,fr)

       call rotate1(3*maxatom,
     1  rotmatr,zero,zero,zero,
     2	 pivotnull, natoml,fl)

c calculate rot/transl delta for receptor/ligand

       call cartstate_select_ligand(cartstatehandle,idr,
     1  natomr,nresr,ptr_ieir,ptr_xr,ptr_fr,pivotr,ptr_iacir,
     2  ptr_icopr,ptr_wer,ptr_chair,ptr_ncopr,
     3  ptr_nmaxcor,ptr_natcor, false)

       if (fixre.eq.0) then
       call euler2torquemat(phir,ssir,rotr,pm2)
       call trans(3*maxatom,maxdof,fr,deltar,natomr)
       call rota(3*maxatom,maxdof,
     1  xr,fr,deltar,pm2,natomr)
       do 20 i=1,maxdof
       deltar(i) = deltar(i)+deltar0(i)
20     continue  
       endif
       
       call euler2torquemat(phil,ssil,rotl,pm2)
       call rota(3*maxatom,maxdof,
     1  xl0,fl,deltal,pm2,natoml)       
       call trans(3*maxatom,maxdof,fl,deltal,natoml)       

      e = zero
      do 990,i=1,6
      e = e + energies(i)
990   continue  
c     e = nonbonded+potential+erest+esolv+enlig+eem 
c      write(*,*),'deltar',deltar(1),deltar(2),deltar(3),
c     1 deltar(4),deltar(5),deltar(6)
c      write(*,*),'deltal',deltal(1),deltal(2),deltal(3),
c     1 deltal(4),deltal(5),deltal(6)
      
      return
      end

      subroutine printmat(m) 
       real*8 m      
       dimension m(9)
       write (*,*), m(1),m(2),m(3)
       write (*,*), m(4),m(5),m(6)
       write (*,*), m(7),m(8),m(9)
      end
