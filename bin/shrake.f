      program shrake
      implicit none
c
      character name1*80,name2*80,at*4,b*80
      character*4 ty(30000),rg(30000)
c
      integer i,j,im,natoms,nt,nflex,nstru,nphi,ncosth,iemax
      integer ongrid(30000),neigh(5000),numneh,ie(30000),ire(5000)
      real xcrd(30000), ycrd(30000), zcrd(30000),srad,cha
      real radi(30000),csth(80),snth(80),phgh(80)
      real sasa(30000),sres(6000)
c
c
      real xp,yp,zp,ccdist,ges,rr,gesd,gesp
      real phigh, costh, sinth
      logical coverd
      integer ii,jj,kk,ll,ij,ik,idna,iold
      real pi
      pi= 3.141592654
c
c
c   precalculate costh,sinth and phigh for different ncosth...
c
      nphi =  30
      ncosth = 15
      do 350,ik=1,ncosth
          csth(ik) = -1.0 + (float(ik)-0.5)*(2.0/float(ncosth))
          snth(ik) = sqrt(1.0-csth(ik)*csth(ik))
  350 continue
c
      do 360,ik=1,nphi
          phgh(ik) = 2.0*pi*float(ik-1)/float(nphi)
  360 continue 
c
c read in data
c
      call getarg(1,name1)
      call getarg(2,name2)
      read(name2,*) srad
      srad = 1.4
      do 400 i=1,30000
      radi(i)=0.0
  400 continue
      open(42,file=name1)
      i=0
      iold=-10
      idna=0
   40 read(42,20,end=30) b
   20 format(a80)
      if(b(:4).eq.'ATOM') then
      i=i+1
      read(b,26) at,kk,ty(i),rg(i),ie(i),xcrd(i),ycrd(i),zcrd(i)
      if(ty(i)(:1).eq.'C'.or.ty(i)(2:2).eq.'C') radi(i)=1.7 
      if(ty(i)(:1).eq.'H'.or.ty(i)(2:2).eq.'H') radi(i)=0.0 
      if(ty(i)(:1).eq.'P'.or.ty(i)(2:2).eq.'P') radi(i)=2.0 
      if(ty(i)(:1).eq.'S'.or.ty(i)(2:2).eq.'S') radi(i)=2.0 
      if(ty(i)(:1).eq.'N'.or.ty(i)(2:2).eq.'N') radi(i)=1.6
      if(ty(i)(:1).eq.'O'.or.ty(i)(2:2).eq.'O') radi(i)=1.5
      if(ie(i).ne.iold) then
      iold=ie(i)
      idna=idna+1
      endif 
      endif 
      ire(i)=idna 
      goto 40
   30 close(42)
      natoms=i
      iemax=idna
      write(*,*) 'Zahl',natoms,iemax
      do 50 i=1,natoms
      ongrid(i)=1
      if(radi(i).lt.0.001) write(*,*)'atom ',i,' has no radius!'
   50 continue
c    
c     go over all atoms in list
      do 100 ii = 1, natoms
        sasa(ii) = 0.0
        if (ongrid(ii).eq.1) then
c       make a neighbor list 
        numneh = 1
        do 200 jj = 1, natoms
        if (jj.ne.ii) then
            ccdist =   ( xcrd(ii) - xcrd(jj) )**2
     .                     +( ycrd(ii) - ycrd(jj) )**2
     .                     +( zcrd(ii) - zcrd(jj) )**2 
        rr=(radi(ii)+radi(jj)+2.0*srad)*(radi(ii)+radi(jj)+2.0*srad) 
        if (ccdist.le.rr) then
              neigh(numneh) = jj
              numneh = numneh + 1
c             write(*,*)'numneh', numneh
              if (numneh.gt.200) then
                  write(*,*) 'ERROR: Atom has too many neighbors'
              endif
            endif
          endif
200     continue
        numneh = numneh - 1
c       now go over cosine(thetas) and do phi circle for each.
c
        do 500 jj = 1, ncosth
          costh = csth(jj)
          sinth = snth(jj)
          do 600 kk = 1, nphi
            phigh = phgh(kk)
            xp = (radi(ii)+srad)*sinth*cos(phigh)
            yp = (radi(ii)+srad)*sinth*sin(phigh)
            zp = (radi(ii)+srad)*costh
c           go over neighbors
            ll = 1
            coverd = .false.
700         if( (.not.coverd).and.(ll.le.numneh) ) then
              if (  ( xcrd(ii) + xp - xcrd(neigh(ll)) ) **2
     .                  +( ycrd(ii) + yp - ycrd(neigh(ll)) ) **2
     .                  +( zcrd(ii) + zp - zcrd(neigh(ll)) ) **2 .lt.
     .               ( radi(neigh(ll)) + srad )**2 )  coverd = .true.
c
              ll = ll + 1
              goto 700
            endif
            if (.not.coverd) sasa(ii) = sasa(ii) + 1.0
600       continue
500     continue
        sasa(ii) = (radi(ii)+srad)**2*(2.0/float(ncosth))*
     .             (2.0*pi/nphi)*sasa(ii)
      endif
100   continue
      ges=0.0
      do 78 i=1,iemax
      sres(i)=0.0
   78 continue
c     write(*,*)'sas for each atom'
      do 80 i=1,natoms
      sres(ire(i))=sres(ire(i))+sasa(i)
      ges=ges+sasa(i)
   80 continue
      write(*,*)'for each residue'
      do 85 i=1,iemax
      write(*,*) i,sres(i)
   85 continue       
      write(*,*) 'SAS:',ges
c 
   26 format(a4,i7,1x,a4,1x,a4,2x,i4,4x,3f8.3)
      end
