       subroutine pairgen(maxatom,maxres,maxmolpair,
     1  molpairhandle,cartstatehandle,rcut)
      
       implicit none

c      Parameters
       integer maxatom,maxres,maxmolpair,
     1  molpairhandle,cartstatehandle
       real*8 rcut

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
       dimension fl(maxatom),pivotl(3)       
       integer nmaxcol,ieil,iacil,natcol,icopl
       dimension nmaxcol(maxres),ieil(maxatom),
     1 iacil(maxatom),natcol(maxres),icopl(maxatom)
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
c      paramters for softcore potential
       integer use_softcore
       real*8 softcore
c      Handle variables: forcefield parameters
       integer potshape
       real swi_on, swi_off
       real*8 rbc,rc,ac,emin,rmin2,xnull,vlj
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

c     Local variables       
       integer dmmy1
       real*8 dmmy2
       integer i,j,jj,jjj,k,ii,iii,it,jt,false,true
       real*8 rcut0, rcut1,rd
              
       true = 1
       false = 0
	      
       call molpair_get_values(molpairhandle,idr,idl,
     1  ptr_iactr,ptr_iactl,nonp,ptr_nonr,ptr_nonl)
     
       call cartstate_select_ligand(cartstatehandle,idr,
     1  natomr,nresr,ptr_ieir,ptr_xr,ptr_fr,pivotr,ptr_iacir,
     2  ptr_icopr,ptr_wer,ptr_chair,ptr_ncopr,
     3  ptr_nmaxcor,ptr_natcor,true)
       call cartstate_select_ligand(cartstatehandle,idl,
     1  natoml,nresl,ptr_ieil,ptr_xl,ptr_fl,pivotl,ptr_iacil,
     2  ptr_icopl,ptr_wel,ptr_chail,ptr_ncopl,
     3  ptr_nmaxcol,ptr_natcol,false)
       call cartstate_get_parameters(cartstatehandle,
     1  ptr_rbc,ptr_rc,ptr_ac,ptr_emin,ptr_rmin2,ptr_ipon,potshape,
     2  dmmy1,dmmy2,swi_on, swi_off, use_softcore, softcore)


c This subroutine generates a ligand-receptor pairlist
c based on a cutoff (rcut) and on the active atoms (in case of multiple copies)
c active atoms are in iactr(i):receptor and iactl(j):ligand.
c Atom code 99 is a universal code for 'dummy'
c Atom code 32 can also be a dummy, if no parameters above 31 are specified

       jj=0
       k=1
       rcut0=dsqrt(rcut)

       do 250 i=1,natomr
       if(iactr(i).eq.1.and.iacir(i).ne.0.and.iacir(i).ne.99) then
       iii=3*(i-1)
       it=iacir(i)
       do 255 j=1,natoml
       if(iactl(j).eq.1.and.iacil(j).ne.0.and.iacil(j).ne.99) then
       jjj=3*(j-1)
       jt=iacil(j)
       
c      TODO: use rcut or rcut1? use rcut for now....
c       rcut1=(rcut0+rbc(it,jt))**2 

       rd=(xr(iii+1)-xl(jjj+1))**2+
     1   (xr(iii+2)-xl(jjj+2))**2+
     2   (xr(iii+3)-xl(jjj+3))**2 

c       if(rd.le.rcut1) then
c      write(*,*)'rcut', iii,jjj,rd, rcut1
       if(rd.le.rcut) then
       jj=jj+1
       if (jj.ge.maxmolpair) then
       write(*,*) 'ERROR: maximum number of pairs reached', maxmolpair
       stop
       endif
       nonr(jj)=i
       nonl(jj)=j
       endif
       endif
  255  continue
       endif
  250  continue
       nonp=jj
c       write(*,*)'nonp',nonp
       
       call molpair_set_nonp(molpairhandle,nonp)
       
       return
       end
