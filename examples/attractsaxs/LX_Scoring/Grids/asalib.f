      subroutine asa(xcrd, ycrd, zcrd, radi, srad, natoms, ges)
cf2py intent(out) :: ges
cf2py intent(hide) :: natoms
      real radi(natoms)
      real modatyp(natoms)
c
      integer i, j, natoms, nphi, ncosth
      integer neigh(5000),numneh,ongrid(30000)
      real xcrd(natoms), ycrd(natoms), zcrd(natoms),srad
      real csth(80),snth(80),phgh(80)
      real sasa(30000)
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
c    
c     go over all atoms in list
      do 50 i=1,natoms
      ongrid(i)=1
   50 continue
      do 100 ii = 1, natoms
        sasa(ii) = 0.0
        if (ongrid(ii).eq.1) then
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
      do 80 i=1,natoms
      ges=ges+sasa(i)
   80 continue     
c      write(*,*) 'SAS:',ges
c 
   26 format(a4,i7,1x,a4,1x,a4,2x,i4,4x,3f8.3)
      end
c
c
c     calculate surfacearea of atom-type
c
      subroutine asaatom(xcrd,ycrd,zcrd,radi,atyp,natyp,srad,natoms,ges)
cf2py intent(out) :: ges
cf2py intent(hide) :: natoms
      real radi(natoms)
      real modatyp(natoms)
      integer atyp(natoms)
c
      integer i, j, u, natoms, nphi, ncosth, natyp
      integer neigh(5000),numneh,ongrid(30000)
      real xcrd(natoms), ycrd(natoms), zcrd(natoms),srad
      real csth(80),snth(80),phgh(80)
      real sasa(30000)
c
c
      real xp,yp,zp,ccdist,ges(natyp),rr,gesd,gesp
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
c    
c     go over all atoms in list
      do 50 i=1,natoms
      ongrid(i)=1
   50 continue
      do 100 ii = 1, natoms
        sasa(ii) = 0.0
        if (ongrid(ii).eq.1) then
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
      do 800 u=1, natyp
      ges(u)=0.0
800   continue      
      do 80 i=1,natoms
      ges(atyp(i))=ges(atyp(i))+sasa(i)
   80 continue     
c      write(*,*) 'SAS:',ges
c 
   26 format(a4,i7,1x,a4,1x,a4,2x,i4,4x,3f8.3)
      end
