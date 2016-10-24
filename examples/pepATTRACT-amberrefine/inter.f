      program interface
      character*80 name1,name2,name3,name4,filnam,b
      character*4 tt(150000),rg(150000),at
      integer re(150000),ka(150000),iai(150000),im1(150000),im2(150000)
      integer io1(51000),io2(51000)
      real x(500000),ceca(3)
c  Einlesen des pdb-files und Generierung des pdb-Output files
c   mit dem integer-type code und charges
      call getarg(1,name1)
      call getarg(2,name2)
      open(43,file=name1)
      read(name2,*) iflag 
      i=1
      k=1
      iatom=0
   40 read(43,20,end=35) b
   20 format(a80)
      if(b(:3).eq.'TER'.and.iatom.eq.1) then
      isplit=k-1
      iatom=2
      endif
      if(b(:4).eq.'ATOM') then
      if(iatom < 1) iatom=1
      ii=3*(i-1)
      read(b,26) at,ka(i),tt(i),rg(i),re(i),
     1           x(ii+1),x(ii+2),x(ii+3)
      if(tt(i).eq."CA  ") then
c     if(tt(i).eq."CA  ".or.tt(i)(:2).eq."ND".or.tt(i)(:2).eq."NE") then
      iai(k)=i
      k=k+1
      endif
      i=i+1
      endif 
      goto 40
   35 close(43)
      nca=k-1
      natom=i-1
      nat3=3*natom
c     write(*,*)'natom',nca,natom,isplit
c  determine equi distances for first protein
      do i=3,isplit-3,2
      ii=3*(iai(i)-1)
      do j=1,isplit-2
      jj=3*(iai(j)-1)
      dis=sqrt((x(ii+1)-x(jj+1))**2+(x(ii+2)-x(jj+2))**2
     1+(x(ii+3)-x(jj+3))**2)
      if(dis.lt.12.0.and.dis.gt.7.0) then
      write(*,174) ka(iai(i)),ka(iai(j)),dis-3.5,dis,dis,
     1             dis+3.5,2.0,2.0
      endif
      end do
      end do
c   determine equi distances for second protein
      if(iflag.eq.1) then
      do i=isplit+2,nca-1,2
      ii=3*(iai(i)-1)
      do j=isplit+1,nca-2,2
      jj=3*(iai(j)-1)
      dis=sqrt((x(ii+1)-x(jj+1))**2+(x(ii+2)-x(jj+2))**2
     1+(x(ii+3)-x(jj+3))**2)
      if(dis.lt.12.0.and.dis.gt.7.0) then
      write(*,174) ka(iai(i)),ka(iai(j)),dis-3.5,dis,dis,
     1             dis+3.5,2.0,2.0
      endif
      end do
      end do
      endif
c this is for limiting the mobility of the ligand relative to receptor 
      do i=1,nca
      im1(i)=0
      im2(i)=0
      end do
      do i=1,isplit-2,4
      ii=3*(iai(i)-1)
      do j=isplit+2,nca-1,4
      jj=3*(iai(j)-1)
      dis=sqrt((x(ii+1)-x(jj+1))**2+(x(ii+2)-x(jj+2))**2
     1+(x(ii+3)-x(jj+3))**2)
      if(dis.lt.20.0.and.dis.gt.15.0) then
      write(*,174) ka(iai(i)),ka(iai(j)),dis-13.5,dis-10.0,dis+10.0,
     1             dis+13.5,0.25,0.25
      endif
      end do
      end do
c this is for center to center pushing
      xcen1=0.0
      ycen1=0.0
      zcen1=0.0
      xcen2=0.0
      ycen2=0.0
      zcen2=0.0
      icen=0
      do i=1,nca
      im1(i)=0
      im2(i)=0
      end do
      do i=1,isplit-2,1
      ii=3*(iai(i)-1)
      do j=isplit+2,nca-1,1
      jj=3*(iai(j)-1)
      dis=sqrt((x(ii+1)-x(jj+1))**2+(x(ii+2)-x(jj+2))**2
     1+(x(ii+3)-x(jj+3))**2)
      if(dis.lt.20.0.and.dis.gt.10.0) then
      im1(i)=1
      im2(j)=1
      xcen1=xcen1+x(ii+1)
      xcen2=xcen2+x(jj+1)
      ycen1=ycen1+x(ii+2)
      ycen2=ycen2+x(jj+2)
      zcen1=zcen1+x(ii+3)
      zcen2=zcen2+x(jj+3)
      icen=icen+1
      endif
      end do
      end do
      xcen1=xcen1/float(icen)
      ycen1=ycen1/float(icen)
      zcen1=zcen1/float(icen)
      xcen2=xcen2/float(icen)
      ycen2=ycen2/float(icen)
      zcen2=zcen2/float(icen)
      dclold=100000000.0
      do i=1,isplit-2,1
      ii=3*(iai(i)-1)
      dcl=(x(ii+1)-xcen1)**2+(x(ii+2)-ycen1)**2+(x(ii+3)-zcen1)**2
      if(dcl.lt.dclold) then
      iold=i
      dclold=dcl
      endif
      enddo 
      dclold=100000000.0
      do j=isplit+2,nca-1,1
      jj=3*(iai(j)-1)
      dcl=(x(jj+1)-xcen2)**2+(x(jj+2)-ycen2)**2+(x(jj+3)-zcen2)**2
      if(dcl.lt.dclold) then
      jold=j
      dclold=dcl
      endif
      enddo 
      write(*,174) ka(iai(iold)),ka(iai(jold)),0.0,0.0,0.0,
     1             90.0,0.50,0.50
  174 format('&rst iat=',i5,',',i5,',r1=',f4.1,',r2=',f5.2,',r3=',f5.2,
     1',r4=',f4.1,',rk2=',f5.2,',rk3=',f5.2,',ir6=0',',ialtd=0, /')
  178 format('r1=',f6.3,',','r2=',f6.3,',','r3=',f6.3,',',
     1       'r4=',f6.3,',','rk2=',f6.3,',','rk3=',f6.3,',')
  180 format(' igr1=',50(i5,','))
  181 format(' igr2=',50(i5,','))
   25 format(a4,i7,1x,a4,1x,a4,i5,4x,3f8.3)
   26 format(a4,i7,2x,a4,a4,1x,i4,4x,3f8.3)
 2000 end
