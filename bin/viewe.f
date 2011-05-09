      program eigenr
      character*80 name1,name2,name3,name4,filnam,b
      character*4 tt(1500),rg(1500),at
      integer re(1500),ka(1500),iai(1500)
      real val(4500),x(4500),y(4500,4500),cha(4500)
c  Einlesen des pdb-files und Generierung des pdb-Output files
c   mit dem integer-type code und charges
      call getarg(1,name2)
      read(name2,*) delta
      call getarg(2,name3)
      call getarg(3,name4)
      open(43,file=name3)
      i=1
   40 read(43,20,end=35) b
   20 format(a80)
      if(b(:4).eq.'ATOM') then
      ii=3*(i-1)
      read(b,26) at,ka(i),tt(i),rg(i),re(i),
     1           x(ii+1),x(ii+2),x(ii+3)
      write(*,26) at,ka(i),tt(i),rg(i),re(i),
     1           x(ii+1),x(ii+2),x(ii+3)
      i=i+1
      endif 
      goto 40
   35 close(43)
      natom=i-1
      nat3=3*natom
c     if(381.ne.nat3) goto 1000 
      itest=int(nat3/6)
      irest=nat3-6*itest
      write(*,*) nat3,itest,irest
      open(45,file=name4)
      do l=1,nat3
      read(45,*) kk,val(l) 
      do i=1,itest
      read(45,*)(y(l,6*(i-1)+j),j=1,6)
      end do
      if(irest.gt.0) read(45,*)(y(l,6*itest+j),j=1,irest)
      write(*,*) 'eigenvector read ',l,kk,val(l)
      end do
      close(45)
      write(*,*)'eigenvalues',(val(i),i=nat3-6,nat3-15,-1)
      ikk=0
      open(45,file='eignew.out')
      do in=nat3-6,nat3-25,-1
      write(45,*) in,val(in) 
      del=(1.0/val(in))*delta
      ikk=ikk+1
      vor=1.0
      do 1020 k=1,itest
      kk=6*(k-1)
      write(45,'(6f15.10)') (y(in,kk+j),j=1,6)
 1020 continue
      if(irest.gt.0) write(45,'(3f15.10)') (y(in,3*(nat3-1)+j),j=1,3)
      write(filnam,177) ikk
  177 format('mode',i3.3,'.pdb')
      open(55,file=filnam)
      write(55,'(a11)')'Header from'
c     do ijl=1,15
      del=-1.0*del
      do m=-3,3
      do i=1,natom
      ii=3*(i-1)
      xx=x(ii+1)+m*vor*del*y(in,ii+1)
      xy=x(ii+2)+m*vor*del*y(in,ii+2)
      xz=x(ii+3)+m*vor*del*y(in,ii+3)
      write(55,26) at,i,tt(i),rg(i),re(i),xx,xy,xz
      end do
      write(55,'(a3)')'TER'
      write(55,'(a3)')'END'
      end do
c     end do
      close(55)
      end do
   25 format(a4,i7,1x,a4,1x,a4,i5,4x,3f8.3)
   26 format(a4,i7,2x,a4,a4,2x,i3,4x,3f8.3,i5,f8.3)
      goto 2000
 1000 write(*,*)'muell: nat3.ne.j',nat3,j
 2000 end
