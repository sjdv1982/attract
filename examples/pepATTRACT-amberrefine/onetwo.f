      program onetwo 
      character*80 name1,name2,name3,name4,filnam,b
      integer re(25000),ka(25000),iai(25000)
      real*8 x(75000),ceca(3)
c  Einlesen des pdb-files und Generierung des pdb-Output files
c   mit dem integer-type code und charges
      call getarg(1,name1)
      call getarg(2,name2)
      call getarg(3,name3)
      call getarg(4,name4)
      open(43,file=name1)
      read(name2,*) isplit
      i=1
      k=1
      read(43,*) 
      read(43,'(i5)') natom
      iup=natom - isplit
      natom=3*natom 
      iup=3*iup
      isplit=3*isplit
      isix=int(natom/6)
      irest=int(natom-6*isix)
c     write(*,*),isix,irest
      do i=1,isix
      ii=6*(i-1)
      read(43,'(6f12.7)') (x(ii+j),j=1,6)
c     write(*,'(6f12.7)') (x(ii+j),j=1,6)
      end do
      if(irest.ne.0) read(43,'(3f12.7)') (x(6*isix+j),j=1,irest)
      close(43)
      isix=int(isplit/6)
      irest=isplit-6*isix
      open(44,file=name3)
      write(44,*)
      write(44,'(i5)') isplit/3
      do i=1,isix
      ii=6*(i-1)
      write(44,'(6f12.7)') (x(ii+j),j=1,6)
      end do
      if(irest.gt.0) write(44,'(6f12.7)') (x(6*isix+j),j=1,irest)
      close(44)
      open(44,file=name4)
      isix=int(iup/6)
      irest=iup-6*isix
      write(44,*) 
      write(44,'(i5)') iup/3
      do i=1,isix
      ii=6*(i-1)+isplit
      write(44,'(6f12.7)') (x(ii+j),j=1,6)
      end do
      if(irest.gt.0)write(44,'(6f12.7)')(x(6*isix+j+isplit),j=1,irest)
      close(44)
      end
