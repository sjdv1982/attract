      program compare
c
c program to generate a structure deformed in a subset of normal modes
c to optimally fit to a target structure.
c usage: $path/compare start.pdb target.pdb eigen.dat  number of modes 
c (eigen.dat contains the normal modes calculated for structure start.pdb)
c author: Martin Zacharias, Jacobs University Bremen
c
      character*80 name1,name2,name3,name4,filnam,b
      character*4 tt(2000),rg(2000),at
      integer re(2000),ka(2000),iai(2000)
      real val(2000),x(2000),z(2000),y(6000,6000),cha(2000),
     1     zn(2000),dxz(2000),skalar(50)
c  Einlesen des pdb-files und Generierung des pdb-Output files
c   mit dem integer-type code und charges
      call getarg(1,name1)
      call getarg(2,name2)
      call getarg(3,name3)
      call getarg(4,name4)
      read(name4,*) ne
      open(43,file=name1)
      i=1
   40 read(43,20,end=35) b
      if(b(:4).eq.'ATOM') then
      ii=3*(i-1)
      read(b,26) at,ka(i),tt(i),rg(i),re(i),
     1           x(ii+1),x(ii+2),x(ii+3)
c     write(*,26) at,ka(i),tt(i),rg(i),re(i),
c    1           x(ii+1),x(ii+2),x(ii+3)
      i=i+1
      endif 
      goto 40
   35 close(43)
      open(44,file=name2)
      i=1
   42 read(44,20,end=45) b
      if(b(:4).eq.'ATOM') then
      ii=3*(i-1)
      read(b,26) at,ka(i),tt(i),rg(i),re(i),
     1           z(ii+1),z(ii+2),z(ii+3)
c     write(*,26) at,ka(i),tt(i),rg(i),re(i),
c    1           z(ii+1),z(ii+2),z(ii+3)
      i=i+1
      endif 
      goto 42
   45 close(44)
      natom=i-1
      nat3=3*natom
c determine difference vector between the two structures
      do i=1,nat3
      dxz(i)=x(i)-z(i)
      end do
c start reading eigenvectors
      itest=int(nat3/6)
      irest=nat3-6*itest
      write(*,*) nat3,itest,irest
      open(45,file=name3)
      do l=1,nat3
      read(45,*) kk,val(l) 
      do i=1,itest
      read(45,*)(y(l,6*(i-1)+j),j=1,6)
      end do
      if(irest.gt.0) read(45,*)(y(l,6*itest+j),j=1,irest)
c     write(*,*) 'eigenvector read ',l,kk,val(l)
      end do
      close(45)
      write(*,*)'eigenvalues',(val(i),i=nat3-6,nat3-25,-1)
c do the scalar product between ne softest modes and structural difference vector
      ikk=0
      do i=1,nat3
      zn(i)=0.0
      end do
      do in=nat3-6,nat3-5-ne,-1
      ikk=ikk+1
      skalar(ikk)=0.0 
      do i=1,nat3
      skalar(ikk)=skalar(ikk)+dxz(i)*y(in,i) 
      end do
c     skalar(ikk)=skalar(ikk)/nat3
      write(*,*)' contribution of: ',ikk,skalar(ikk)
c generate best possible approximation to second structure
      do i=1,nat3
      zn(i)=zn(i)-skalar(ikk)*y(in,i)
      end do
      end do  
      open(55,file='final.pdb')
      do i=1,natom
      ii=3*(i-1)
      write(55,26) at,i,tt(i),rg(i),re(i),
     1    x(ii+1)+zn(ii+1),x(ii+2)+zn(ii+2),x(ii+3)+zn(ii+3)
      end do
      write(55,'(a3)')'TER'
      write(55,'(a3)')'END'
      close(55)
   20 format(a80)
   25 format(a4,i7,1x,a4,1x,a4,i5,4x,3f8.3)
   26 format(a4,i7,2x,a4,a4,2x,i3,4x,3f8.3,i5,f8.3)
      goto 2000
 1000 write(*,*)'muell: nat3.ne.j',nat3,j
 2000 end
