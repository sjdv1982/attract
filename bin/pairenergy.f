       subroutine pairenergy(
     1 cartstatehandle,molpairhandle,iab,fixre,gridptr,
     2 ffr, frotr, ffl, frotl,
     3 ensr, phir, ssir, rotr, xar, yar, zar, morphr, dligr,
     4 ensl, phil, ssil, rotl, xal, yal, zal, morphl, dligl,
     5 energies, deltar, deltal, deltamorphr, deltamorphl)

       implicit none

c      Parameters
       include "max.fin"
       integer cartstatehandle,molpairhandle
       integer iab, fixre
       integer gridptr_dmmy
       pointer(gridptr, gridptr_dmmy)
       integer ensl, ensr
       real*8 ffr, frotr, ffl, frotl
       dimension frotr(9), frotl(9)
       real*8 phir, ssir, rotr, xar, yar, zar, morphr, dligr
       real*8 phil, ssil, rotl, xal, yal, zal, morphl, dligl
       dimension dligr(maxmode+maxindexmode)
       dimension dligl(maxmode+maxindexmode)
       real*8 energies
       dimension energies(6)
       real*8 deltar(6+maxmode+maxindexmode)
       real*8 deltal(6+maxmode+maxindexmode)
       real*8 deltamorphr, deltamorphl
       
      
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
       dimension xr(3*maxatom),wer(maxatom),chair(maxatom)
       dimension fr(3*maxatom), pivotr(3)
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
       dimension xl(3*maxatom),wel(maxatom),chail(maxatom)
       dimension fl(3*maxatom), pivotl(3)
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

       real*8 ensdr
       dimension ensdr(3*maxatom)
       pointer(ptr_ensdr, ensdr)
       real*8 ensdl
       dimension ensdl(3*maxatom)
       pointer(ptr_ensdl, ensdl)
       real*8 cmorphr, cmorphdr
       pointer(ptr_cmorphdr,cmorphdr)
       dimension cmorphdr(3*maxatom)
       real*8 cmorphl, cmorphdl
       pointer(ptr_cmorphdl,cmorphdl)
       dimension cmorphdl(3*maxatom)
c      variables to integrate softcore potential
       integer use_softcore
       pointer(ptr_use_softcore,use_softcore)
       real*8 softcore
       pointer(ptr_softcore,softcore)

c      Handle variables: forcefield parameters
       integer potshape, cdie
       real swi_on, swi_off
       real*8 rbc,rc,ac,emin,rmin2,vlj,epsilon
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
       integer index_eig
       real*8 index_val
       integer ieins,nhm,nihm
       dimension xb(3*totmaxatom),x(3*totmaxatom)
       dimension xori(3*totmaxatom),xori0(3*totmaxatom)
       dimension eig(maxmode,3*maxatom,maxlig)
       dimension index_eig(maxlenindexmode,maxindexmode,maxlig)
       dimension index_val(maxlenindexmode,maxindexmode,maxlig)
       dimension ieins(maxlig),nhm(maxlig), nihm(maxlig)
       pointer(ptr_xb,xb)
       pointer(ptr_x,x)
       pointer(ptr_xori,xori)
       pointer(ptr_xori0,xori0)       
       pointer(ptr_eig,eig)
       pointer(ptr_index_eig,index_eig)
       pointer(ptr_index_val,index_val)
       pointer(ptr_ieins,ieins)
       pointer(ptr_nhm,nhm)
       pointer(ptr_nihm,nihm)
       integer dmmy
       pointer(ptr_dmmy, dmmy)

c      Local variables: matrices
       real*8 rotmatr,rotmatl,rotmatd,rotmatd2,rotmatrinv,rotmatlinv
       dimension rotmatr(9),rotmatl(9),rotmatd(9),rotmatrinv(9),
     1  rotmatlinv(9),rotmatd2(9)
       real*8 tr,tl,td,td0,xnull
       dimension tr(3),tl(3),td(3),td0(3)
       real*8 pivotd,pivotd0,pivotnull
       dimension pivotd(3),pivotd0(3),pivotnull(3)
       real*8 pm2(3,3,3), pm20(3,3,3), pm2a(3,3,3), mat(9)
       real*8 frotrinv, frotlinv
       dimension frotrinv(9), frotlinv(9)

c  energies is an array of double with a value for every energy type
c  (vdw, elec, ...) for this pair
c
c phir, ssir, ssil etc. are scalars
c dligr / dligl is an array of nhm[ligand]/nhm[receptor]
c deltar/deltar is an array of 6 + nhm[ligand]/nhm[receptor]

c     Local variables: other      
      integer i, j,l,n, use_grid, torquegrid, rigid, true, false
      integer dummy2
      real*8 enon,epote,zero,e
      real*8 flcopy,xl0
      dimension flcopy(3*maxatom),xl0(3*maxatom)
      real*8 deltar0(6+maxmode+maxindexmode)
      integer, parameter :: ERROR_UNIT = 0
      
      deltamorphr = 0
      deltamorphl = 0
      
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
       rigid = 0
       use_grid = 1
       endif

c get molpair 

       call molpair_get_rl(molpairhandle,idr,idl,ptr_dmmy,
     1  torquegrid,dummy2)

c get full coordinates and normal modes

       call cartstate_f_pairenergy(cartstatehandle,ptr_nhm,ptr_nihm,
     1  ptr_ieins, ptr_eig, ptr_index_eig, ptr_index_val,
     2  ptr_xb, ptr_x, ptr_xori, ptr_xori0)

c apply ensemble/normal mode deformations to receptor and ligand
       call cartstate_get_ensd(cartstatehandle,idr,ensr,ptr_ensdr,
     1  morphr,cmorphr,ptr_cmorphdr)
       call deform(ensr,ensdr,cmorphr,cmorphdr,dligr,nhm,nihm,idr,ieins,
     1  eig,index_eig,index_val,xb,x,xori,xori0,1)
       call cartstate_get_ensd(cartstatehandle,idl,ensl,ptr_ensdl,
     1  morphl,cmorphl,ptr_cmorphdl)
       call deform(ensl,ensdl,cmorphl,cmorphdl,dligl,nhm,nihm,idl,ieins,
     1  eig,index_eig,index_val,xb,x,xori,xori0,1)

c select deformed coordinates: select non-pivotized coordinates for receptor

       call cartstate_select_ligand(cartstatehandle,idr,
     1  natomr,nresr,ptr_ieir,ptr_xr,ptr_fr,pivotr,ptr_iacir,
     2  ptr_icopr,ptr_wer,ptr_chair,ptr_ncopr,
     3  ptr_nmaxcor,ptr_natcor, true)

       call cartstate_select_ligand(cartstatehandle,idl,
     1  natoml,nresl,ptr_ieil,ptr_xl,ptr_fl,pivotl,ptr_iacil,
     2  ptr_icopl,ptr_wel,ptr_chail,ptr_ncopl,
     3  ptr_nmaxcol,ptr_natcol, false)
       xl0(1:3*natoml) = xl(1:3*natoml)
     
       call cartstate_get_parameters(cartstatehandle,
     1  ptr_rbc,ptr_rc,ptr_ac,ptr_emin,ptr_rmin2,ptr_ipon,potshape,
     2  cdie,epsilon,swi_on, swi_off, ptr_use_softcore, ptr_softcore)
            
       xnull = 0.0d0   
       do 10 i=1,6+maxmode+maxindexmode
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
       do j=1,3
       do l=1,3
       mat = rotmatrinv
       pm2(1,j,l) = mat(1)*pm20(1,j,l)+mat(2)*pm20(2,j,l)
     1  + mat(3)*pm20(3,j,l)
       pm2(2,j,l) = mat(4)*pm20(1,j,l)+mat(5)*pm20(2,j,l)
     1  + mat(6)*pm20(3,j,l)
       pm2(3,j,l) = mat(7)*pm20(1,j,l)+mat(8)*pm20(2,j,l)
     1  + mat(9)*pm20(3,j,l)
       enddo
       enddo
       endif

       call nonbon_grid(gridptr,torquegrid,rigid,iab,fixre,
     1  xl,xr,pivotr,tr,
     2  wel,wer,chail,chair,iacil,iacir,natoml,natomr,
     3  rc,ac,emin,rmin2,ipon,potshape,cdie,epsilon,
     4  swi_on, swi_off,
     5  fl,enon,epote,fr,pm2,deltar0) 
c      rotate delta-translate back into global frame

       call rotate1(3*maxatom,
     1  rotmatr,zero,zero,zero,
     2	 pivotnull, 1,deltar0(4))     
       else  
       call molpair_pairgen(molpairhandle,cartstatehandle,nonp)
       call molpair_get_values(molpairhandle,idr,idl,
     1  ptr_iactr,ptr_iactl,nonp,ptr_nonr,ptr_nonl) 
c       write(ERROR_UNIT,*), xl(3*132+1:3*132+3), xr(3*1704+1:3*1704+3)
       if (use_softcore.eq.0) then
      call nonbon8(iab,xl,xr,fl,fr,wel,wer,chair,chail,ac,rc,
     1     emin,rmin2,iacir,iacil,nonr,nonl,ipon,nonp,
     2     potshape,cdie,swi_on,swi_off,enon,epote,natomr,natoml)
       else
       call nonbon_soft(iab,xl,xr,fl,fr,wel,wer,chair,chail,ac,rc,
     1     emin,rmin2,iacir,iacil,nonr,nonl,ipon,nonp,
     2     potshape, cdie, swi_on,swi_off, enon,epote, softcore)
       endif
       endif
       
       energies(1) = enon
       energies(2) = epote
c       write(ERROR_UNIT,*) "Pair", energies(1), energies(2)
c calculate receptor harmonic mode deltas
       call ligmin(fr,natomr,idr+1,eig,nhm(idr+1),deltar)
c       write(ERROR_UNIT,*)'deltar',(deltar(i),i=1,12)

c calculate receptor index mode deltas
       call ligmin_index(fr,natomr,idr+1,index_eig,index_val,
     1  nhm(idr+1),nihm(idr+1),deltar)
c       write(ERROR_UNIT,*)'deltar',(deltar(i),i=1,12)
       call grad_morph(fr,natomr,morphr,cmorphdr,
     1  deltamorphr)       
       
c rotate ligand forces into ligand frame
       call matmult(rotmatr,rotmatlinv,rotmatd2)
       flcopy(:) = fl(:)
       call rotate1(3*maxatom,
     1  rotmatd2,zero,zero,zero,
     2	 pivotnull, natoml,flcopy)
       
c calculate ligand harmonic mode deltas
       call ligmin(flcopy,natoml,idl+1,eig,nhm(idl+1),deltal)
c       write(ERROR_UNIT,*)'deltal',(deltal(i),i=1,12)
c calculate ligand index mode deltas
       call ligmin_index(flcopy,natoml,idl+1,index_eig,index_val,
     1  nhm(idl+1),nihm(idl+1),deltal)
c       write(ERROR_UNIT,*)'deltal',(deltal(i),i=1,12)
       call grad_morph(flcopy,natoml,morphl,cmorphdl,
     1  deltamorphl)
       
c rotate all forces into global frame
       call rotate1(3*maxatom,
     1  rotmatr,zero,zero,zero,
     2	 pivotnull, natomr,fr)
       call forcerotscale(3*maxatom,
     1  ffr,frotr,natomr,fr)

       call rotate1(3*maxatom,
     1  rotmatr,zero,zero,zero,
     2	 pivotnull, natoml,fl)
       call forcerotscale(3*maxatom,
     1  ffl,frotl,natoml,fl)
       
c calculate rot/transl delta for receptor/ligand

       call cartstate_select_ligand(cartstatehandle,idr,
     1  natomr,nresr,ptr_ieir,ptr_xr,ptr_fr,pivotr,ptr_iacir,
     2  ptr_icopr,ptr_wer,ptr_chair,ptr_ncopr,
     3  ptr_nmaxcor,ptr_natcor, false)
       
       if (fixre.eq.0) then
       call euler2torquemat(phir,ssir,rotr,pm2)
c      multiply the torque matrix with the inverse axsym force rotation matrix
       frotrinv(:) = frotr(:)
       call matinv(frotrinv)
       mat = frotrinv
       do j=1,3
       do l=1,3
       pm2a(1,j,l) = mat(1)*pm2(1,j,l)+mat(2)*pm2(2,j,l)
     1  + mat(3)*pm2(3,j,l)
       pm2a(2,j,l) = mat(4)*pm2(1,j,l)+mat(5)*pm2(2,j,l)
     1  + mat(6)*pm2(3,j,l)
       pm2a(3,j,l) = mat(7)*pm2(1,j,l)+mat(8)*pm2(2,j,l)
     1  + mat(9)*pm2(3,j,l)
       enddo
       enddo
       
       call trans(fr,deltar,natomr)
       call rota(xr,fr,deltar,pm2a,natomr)
       do 20 i=1,6+maxmode+maxindexmode
       deltar(i) = deltar(i)+deltar0(i)
20     continue  
       endif
              
       call euler2torquemat(phil,ssil,rotl,pm2)      
c      multiply the torque matrix with the inverse axsym force rotation matrix
       frotlinv(:) = frotl(:)
       call matinv(frotlinv)
       mat = frotlinv
       do j=1,3
       do l=1,3
       pm2a(1,j,l) = mat(1)*pm2(1,j,l)+mat(2)*pm2(2,j,l)
     1  + mat(3)*pm2(3,j,l)
       pm2a(2,j,l) = mat(4)*pm2(1,j,l)+mat(5)*pm2(2,j,l)
     1  + mat(6)*pm2(3,j,l)
       pm2a(3,j,l) = mat(7)*pm2(1,j,l)+mat(8)*pm2(2,j,l)
     1  + mat(9)*pm2(3,j,l)
       enddo
       enddo

       call rota(xl0,fl,deltal,pm2a,natoml)
       call trans(fl,deltal,natoml)
c      write(ERROR_UNIT,*)'deltar',(deltar(i),i=1,12)
c          write(ERROR_UNIT,*)'deltal',(deltal(i),i=1,12)
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
