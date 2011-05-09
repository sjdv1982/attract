      program translate
c
c  program to generate start points on the surface of a receptor structure
c  $path/translate receptor.pdb ligand.pdb > translate.dat
c  author: Martin Zacharias, Jacobs University Bremen
c  employs Shrake and Rupley method to define start positions at receptor surface
c
      character name1*80,at*4,b*80,name2*80
      character*4 ty(8000),rg(8000)
      integer i,j,im,natoms,nt,nflex,nstru,nphi,ncosth,iemax,kkk
      integer ongrid(8000),neigh(8000),numneh,ie(8000)
      real xcrd(8000), ycrd(8000), zcrd(8000),srad,cha,ddd
      real radi(8000),csth(400),snth(400),phgh(400),cx,cy,cz
      real xp(99000),yp(99000),zp(99000),sasa(8000)
c
      real xq,yq,zq,ccdist,ges,rr
      real phigh, costh, sinth
      logical coverd
      integer ii,jj,kk,ll,ij,ik
      real pi
      pi= 3.141592654
c
c begin code
c
      nphi =   90
      ncosth = 60
      do 400 i=1,8000
      radi(i)=0.0
  400 continue
c
c   precalculate costh,sinth and phigh for different ncosth...
c
      do 350,ik=1,ncosth
          csth(ik) = -1.0 + (float(ik)-0.5)*(2.0/float(ncosth))
          snth(ik) = sqrt(1.0-csth(ik)*csth(ik))
  350 continue
c
      do 360,ik=1,nphi
      phgh(ik) = 2.0*pi*float(ik-1)/float(nphi)
  360 continue 
c get name of receptor
      call getarg(1,name1)
c get name of ligand
      call getarg(2,name2)
      open(41,file=name2)
      i=0
      cx=0.0
      cy=0.0
      cz=0.0
   30 read(41,20,end=35) b
   20 format(a80)
      if(b(:4).eq.'ATOM') then
      i=i+1
      read(b,26) at,kk,ty(i),rg(i),ie(i),xcrd(i),ycrd(i),zcrd(i)
c     write(*,26) at,kk,ty(i),rg(i),ie(i),xcrd(i),ycrd(i),zcrd(i)
      cx=cx+xcrd(i) 
      cy=cy+ycrd(i) 
      cz=cz+zcrd(i) 
      endif 
      goto 30
   35 close(42)
      natoms=i
      cx=cx/float(natoms)
      cy=cy/float(natoms)
      cz=cz/float(natoms)
c determine largest distance from center
      dold=0.0
      do 36 i=1,natoms
      ddd=(xcrd(i)-cx)**2+(ycrd(i)-cy)**2+(zcrd(i)-cz)**2
      if (dold.lt.ddd) dold=ddd
   36 continue  
c add 1.0 to mas distance to obtain srad 
      srad=sqrt(dold)+1.0
c     write(*,*)'srad=', natoms,srad
      open(unit=42,file=name1)
      i=0
   40 read(42,20,end=45) b
      if(b(:4).eq.'ATOM') then
      i=i+1
      read(b,26) at,kk,ty(i),rg(i),ie(i),xcrd(i),ycrd(i),zcrd(i)
c     write(*,26) at,kk,ty(i),rg(i),ie(i),xcrd(i),ycrd(i),zcrd(i)
      radi(i)=2.0
      if(ty(i)(:1).eq.'C') radi(i)=2.5 
      endif 
      goto 40
   45 close(42)
      natoms=i
      iemax=ie(i)
      do 50 i=1,natoms
      ongrid(i)=1
c     write(*,*)'radi',natoms,radi(i)
      if(radi(i).le.0.0) write(*,*)'atom ',i,' has no radius!'
   50 continue
c    
c     go over all atoms in list
      ijk=0
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
              if (numneh.gt.8000) then
                  write(*,*) 'ERROR: Atom has too many neighbors'
              endif
            endif
          endif
200     continue
        numneh = numneh - 1
        do 500 jj = 1, ncosth 
          costh = csth(jj)
          sinth = snth(jj)
          do 600 kk = 1, nphi
            phigh = phgh(kk)
            xq = (radi(ii)+srad)*sinth*cos(phigh)
            yq = (radi(ii)+srad)*sinth*sin(phigh)
            zq = (radi(ii)+srad)*costh
c           go over neighbors
            ll = 1
            kkk= 1
            coverd = .false.
700         if( (.not.coverd).and.(ll.le.numneh) ) then
            ddd= ( xcrd(ii) + xq - xcrd(neigh(ll)) ) **2 
     .          +( ycrd(ii) + yq - ycrd(neigh(ll)) ) **2
     .          +( zcrd(ii) + zq - zcrd(neigh(ll)) ) **2
             if(ddd .lt. (radi(neigh(ll)) + srad )**2 ) kkk=kkk+1
             if(ddd .lt. (radi(neigh(ll)) + srad )**2 )  coverd = .true.
              ll = ll + 1
              goto 700
            endif
c      if(.not.coverd.and.kkk.gt.1) then
       if(.not.coverd) then
       ijk=ijk+1 
       xp(ijk)=xcrd(ii) + xq
       yp(ijk)=ycrd(ii) + yq
       zp(ijk)=zcrd(ii) + zq
       endif
600     continue
500     continue
      endif
100   continue
      ipmax=ijk
      call dreieck(xp,yp,zp,ipmax)
      do 900 i=1,natoms
      write(*,26) at,kk,ty(i),rg(i),ie(i),xcrd(i),ycrd(i),zcrd(i)
  900 continue
      write(*,'(a3)') 'TER'
   26 format(a4,i7,2x,a4,a4,1x,i4,4x,3f8.3)
   27 format(a4,i7,2x,a4,a4,1x,i4,4x,3f8.3,i5)
      end
c
      subroutine dreieck (xp,yp,zp,ipmax)
      integer ipmax,neigh(ipmax,30),neimax(ipmax),ielim(ipmax)
      real xp(ipmax),yp(ipmax),zp(ipmax),dis(ipmax,30)
c check distances and neighbors
      do i=1,ipmax
      ielim(i)=0
      enddo
      do i=1,ipmax-1
      k=1
      if(ielim(i).eq.0) then
      do j=i+1,ipmax
      dd=(xp(i)-xp(j))**2+(yp(i)-yp(j))**2+(zp(i)-zp(j))**2 
      if(dd.lt.200.0) then
      ielim(j)=1 
      k=k+1
      endif
      enddo
      endif
      enddo
      k=0
      do i=1,ipmax
      if(ielim(i).eq.0) k=k+1
      enddo
      write(*,'(3i5)') k
      k=0
      do i=1,ipmax
      if(ielim(i).eq.0) then
      write(*,28)"ATOM",k+1,"POSI","PROB",k+1,xp(i),yp(i),zp(i) 
      k=k+1
      endif
      enddo
      write(*,'(a3)')'TER'
   28 format(a4,i7,2x,a4,1x,a3,1x,i4,4x,3f8.3)
      return
      end
