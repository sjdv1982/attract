       subroutine select(maxatom, maxres, maxmolpair,
     1  molpairhandle, cartstatehandle,rcut)
      
       implicit none

c      Parameters
       integer maxatom,maxres,maxmolpair
       integer molpairhandle,cartstatehandle
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
       dimension fl(maxatom), pivotl(3)       
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
c define paramters for use of softcore potential
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

c      Local variables       
       integer dmmy1
       real*8 dmmy2
       real*8 enr, enl, false
       dimension enr(maxres,maxres),enl(maxres,maxres)
       integer ilowr, ilowl
       dimension ilowr(maxres), ilowl(maxres)
       
       integer i,ii,iii,it,j,jj,jjj,jt,k,kk,ipt,ivor,irid
       real*8 r2,rr2,rrd,rr23,rlen,rep,et,enref,charge,alen
       real*8 fswi, r
                     
       false = 0
       
       call molpair_get_values(molpairhandle,idr,idl,
     1  ptr_iactr,ptr_iactl,nonp,ptr_nonr,ptr_nonl)
     
       call cartstate_select_ligand(cartstatehandle,idr,
     1  natomr,nresr,ptr_ieir,ptr_xr,ptr_fr,pivotr,ptr_iacir,
     2  ptr_icopr,ptr_wer,ptr_chair,ptr_ncopr,
     3  ptr_nmaxcor,ptr_natcor, false)
       call cartstate_select_ligand(cartstatehandle,idl,
     1  natoml,nresl,ptr_ieil,ptr_xl,ptr_fl,pivotl,ptr_iacil,
     2  ptr_icopl,ptr_wel,ptr_chail,ptr_ncopl,
     3  ptr_nmaxcol,ptr_natcol, false)
       call cartstate_get_parameters(cartstatehandle,
     1  ptr_rbc,ptr_rc,ptr_ac,ptr_emin,ptr_rmin2,ptr_ipon,potshape,
     2  dmmy1,dmmy2,swi_on, swi_off, use_softcore, softcore)
          
c This subroutine goes through all residues and all side chain copies
c The atoms belonging to the energetically most favorable side chain copy
c are selected and labeled active: iactl(i)=1 (iactr(i)=1:receptor) , atoms
c belonging to high energy copies are labeled iactl(i)=0 
c go through all residue pairs with side chain copies
c and select the best copies for building an atom-atom pairlist
c set energy arrays for receptor and ligand (enr(i,j),enl(j,i)) to zero
c
c TODO: this selection is pair-specific, so in principle, a different
c  active copy can be chosen for interaction pair A-B than for A-C
c
      xnull=0.0d0
c      write(*,*)'nresr,nresl,natomr,natoml',nresr,nresl,natomr,natoml  
      do 50 i=1,nresr
      do 60 j=1,nmaxcor(i)
      enr(i,j)=xnull
   60 continue
   50 continue
      do 70 i=1,nresl
      do 80 j=1,nmaxcol(i)
      enl(i,j)=xnull
   80 continue
   70 continue 
c go through all receptor residues
      do 100 i=1,nresr
c go through side chain copies
      do 105 ii=1,nmaxcor(i)
c go copy atoms
      do 107 jj=1,natcor(i)
      irid=ncopr(ii,jj,i)
      it=iacir(irid)
      iii=3*(irid-1)
      et=xnull
c go over all interactions with the ligand
      do 110 j=1,natoml
      jjj=3*(j-1)
      r2=xnull
      do 120 k=1,3
      r2=r2+(xl(jjj+k)-xr(iii+k))**2
  120 continue
      if(r2.lt.rcut) then
      jt=iacil(j)
      alen=wel(j)*wer(irid)*ac(it,jt)
      rlen=wel(j)*wer(irid)*rc(it,jt)
      ivor=ipon(it,jt)
      charge=wel(j)*wer(irid)*chair(irid)*chail(j)
      r2=r2+0.001d0
      fswi = 1.0d0
      if (swi_on.gt.0.or.swi_off.gt.0) then
        if (r2.ge.(swi_on*swi_on)) then
	  if (r2.ge.(swi_off*swi_off)) then
	    fswi = 0.0d0
	  else
	    r = sqrt(r2) 
	    fswi = 1.0d0 - (r - swi_on)/(swi_off-swi_on)
	  endif
	endif
      endif      
      rr2=1.0d0/r2
      if(abs(charge).gt.xnull) et=charge*rr2
      if (potshape.eq.8) then
      rrd = rr2
      else if (potshape.eq.12) then
      rrd = rr23
      endif
      rep=rlen*rrd
      rr23=rr2**3
      if(r2.lt.rmin2(it,jt)) then
      vlj=(rep-alen)*rr23+(ivor-1)*emin(it,jt)
      else
      vlj=ivor*(rep-alen)*rr23
      endif
      enr(i,ii)=enr(i,ii)+fswi*(vlj+et)
c     write(*,*)'i,j,vlj',i,ii,j,ncopr(ii,jj,i),vlj,et
      endif
  110 continue
  107 continue
  105 continue
  100 continue
c now the same for the ligand protein
      do 200 i=1,nresl
c     write(*,*)'nmaxcol(i),natcol(i)',nmaxcol(i),natcol(i)
      do 205 ii=1,nmaxcol(i)
      do 207 jj=1,natcol(i)
      irid=ncopl(ii,jj,i)
      it=iacil(irid)
c     write(*,*)'irid',irid,it
      iii=3*(irid-1)
      et=xnull
c go over all interactions with the receptor
      do 210 j=1,natomr
      jjj=3*(j-1)
      r2=xnull
      do 220 k=1,3
      r2=r2+(xl(iii+k)-xr(jjj+k))**2
  220 continue
      if(r2.lt.rcut) then
      jt=iacir(j)
      alen=wer(j)*wel(irid)*ac(jt,it)
      rlen=wer(j)*wel(irid)*rc(jt,it)
      ivor=ipon(jt,it)
      charge=wer(j)*wel(irid)*chail(irid)*chair(j)
      r2=r2+0.001d0
      fswi = 1.0d0
      if (swi_on.gt.0.or.swi_off.gt.0) then
        if (r2.ge.(swi_on*swi_on)) then
	  if (r2.ge.(swi_off*swi_off)) then
	    fswi = 0.0d0
	  else
	    r = sqrt(r2) 
	    fswi = 1.0d0 - (r - swi_on)/(swi_off-swi_on)
	  endif
	endif
      endif            
      rr2=1.0d0/r2
      if(abs(charge).gt.xnull) et=charge*rr2

      rr23=rr2**3
      if (potshape.eq.8) then
      rrd = rr2
      else if (potshape.eq.12) then
      rrd = rr23
      endif      
      
      if(r2.lt.rmin2(jt,it)) then      
      rep=rlen*rrd
      vlj=(rep-alen)*rr23+(ivor-1)*emin(jt,it)      
      else     
      rep=rlen*rrd
      vlj=ivor*(rep-alen)*rr23     
      endif
      
      enl(i,ii)=enl(i,ii)+fswi*(vlj+et)
c     write(*,*)'i,j,vlj',i,j,ii,jj,ncopl(ii,jj,i),vlj,et,r2
      endif
  210 continue
c     write(*,*)'i,ii,jj,vlj',i,ii,jj,nmaxcol(i),natcol(i),
c    1 ncopl(ii,jj,i),enl(i,ii),et
  207 continue
  205 continue
  200 continue
c selection of best copy
      do 250 i=1,nresr
      if(nmaxcor(i).ne.0) then
      enref=enr(i,1)
      ilowr(i)=1
      endif
      do 260 j=2,nmaxcor(i)
      if(enr(i,j).lt.enref) then
      enref=enr(i,j)
      ilowr(i)=j
      endif
  260 continue
  250 continue
      do 350 i=1,nresl
      if(nmaxcol(i).ne.0) then 
      enref=enl(i,1)
c     write(*,*)'enl(i,1)',i,enl(i,1)
      ilowl(i)=1
      endif
      do 360 j=2,nmaxcol(i)
c     write(*,*)'enl(i,j)',i,j,enl(i,j)
      if(enl(i,j).lt.enref) then
      enref=enl(i,j)
      ilowl(i)=j
      endif
  360 continue
  350 continue
c determine active atoms
c the active atom list is used in subroutine pairgen to
c generate an active pairlist which is used in nonbon8 to calculate
c energy and forces
      do 400 i=1,natomr
      ipt=ieir(i)
      iactr(i)=0
      if(nmaxcor(ipt).ne.0.and.icopr(i).eq.ilowr(ipt)) iactr(i)=1
      if(nmaxcor(ipt).ne.0.and.icopr(i).eq.0) iactr(i)=1
      if(nmaxcor(ipt).eq.0) iactr(i)=1
  400 continue
      do 450 i=1,natoml
      ipt=ieil(i)
      iactl(i)=0
      if(nmaxcol(ipt).ne.0.and.icopl(i).eq.ilowl(ipt)) iactl(i)=1
      if(nmaxcol(ipt).ne.0.and.icopl(i).eq.0) iactl(i)=1
      if(nmaxcol(ipt).eq.0) iactl(i)=1
c      write(*,*)'selection',i,iactl(i)     
  450 continue
      return
      end
