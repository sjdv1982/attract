      program redatom
c  run by: ./redatom input.pdb output.pdb output.seq
c  needs build_atom.dat in the same directory
      character*80 b,eingabe,ausgabe,ausgabe2
      character*4 names(30),pbtyp(30,30),atyp(37),rg(19000),ty(19000),
     1            rgn(19000),tyn(19000),rgold,anam,at
      character*1 ch(0:19000),chn(0:19000),short(30),sequ(15000)
      integer iatyp(30,37),re(19000),ires(19000),imax(30),iac(19000),
     1        ilam(30,30),ilp(19000),ice(19000),ren(19000),iren(19000),
     2        reold,iacn(19000),icys(30,2)
      real xc(19000),yc(19000),zc(19000),charge(30,30),cha(19000),
     1     xn(19000),yn(19000),zn(19000),chan(19000)
      data(names(j),j=1,30)/'ALA ','ARG ','ASN ','ASP ','CYS ','GLN ',
     1                      'GLU ','GLY ','HIS ','ILE ','LEU ','LYS ',
     2               'LYSH','MET ','PHE ','PRO ','SER ','THR ','TRP ',
     3                      'TYR ','VAL ','ADE ','CYT ','GUA ','THY ',
     4                      'DADE','DCYT','DGUA','DTHY','NET '/
      data(short(j),j=1,30)/'A','R','N','D','C','Q','E','G','H','I','L',
     1                     'K','K','M','F','P','S','T','W','Y','V','Y',
     2                     'C','Z','B','X','X','X','X','X'/
c
      data (atyp(j),j=1,37)/'O   ','OM  ','OA  ','OW  ','N   ','NT  ',
     1      'NL  ','NR5 ','NR5*','NX  ','C   ','CH1 ','CH2 ','CH3 ',
     2      'CR51','CR61','CB  ','H   ','HO  ','HW  ','HT  ','S   ',
     3      'FE  ','ZN  ','NZ  ','NE  ','P   ','OS  ','CS1 ','NR6 ',
     4      'NR6*','CS2 ','SI  ','NA  ','CL  ','CA  ','MG  '/
      call getarg(1,eingabe) 
      call getarg(2,ausgabe)
      call getarg(3,ausgabe2)
      open(42,file=eingabe)
      i=1
      iflater=0
   40 read(42,20,end=30) b
   20 format(a80)
      if(b(:4).eq.'ATOM'.and.b(14:16).ne.'OXT'.and.b(14:14).ne.'H') then
      read(b,25) at,k,ty(i),rg(i),ch(i),re(i),xc(i),yc(i),zc(i)
c  25 format(a4,i7,2x,a4,a4,1x,i4,4x,3f8.3)
   25 format(a4,i7,2x,a4,a4,a1,i4,4x,3f8.3)
      i=i+1
      endif 
      if(b(:3).eq.'TER'.and.iflater.eq.0) then
      iter=i
      iflater=1
      endif
      goto 40
   30 close(42)
      natom=i-1
c automatic addition of hydrogens
c treat first N separately
      reold=100010
      ireold=0
      rgold='    '
      do 150 i=1,natom
      if(re(i).ne.reold.or.rg(i).ne.rgold) then
      reold=re(i)
      rgold=rg(i)
      ireold=ireold+1
      endif
      re(i)=ireold
      if(ty(i).eq.atyp(11)) then
      ires=re(i)
      icarb=i
      endif
  150 continue
      open(45,file='build_atom.dat')
c
c     generierung der pdbcade tabelle und der integercode tabellen
c     fuer jede Aminosaeure;
c
   43 read(45,50,end=80) b
   50 format(a80)
      do 70 i=1,30
      if(b(:4).eq.names(i)) then
      read(45,*) imax(i)
      write(*,*) names(i),imax(i)
      do 55 j=1,imax(i)
      read(45,60) pbtyp(i,j),iatyp(i,j),charge(i,j),ilam(i,j)
      write(*,60) pbtyp(i,j),iatyp(i,j),charge(i,j),ilam(i,j)
   60 format(5x,a4,1x,i4,f10.4,i5)
   55 continue
      endif
   70 continue
      goto 43 
   80 close(45)
      do 200 i=1,natom
       do 210 jj=1,30
         if(names(jj).eq.rg(i)) then
              write(*,*)'name found',i,jj,re(i),rg(i),ty(i),names(jj) 
              ires(i)=jj
          goto 200
        endif
  210 continue
  215 format(1x,'Illegal 3 letter code -')
  200 continue
c    rearrange order in pdb-file if necessary
      j=1
      reold=0
      iseq=0
      do 225 i=1,natom
      if(re(i).ne.reold) then 
      ik=ires(i)
      reold=re(i)
      if(i.eq.iter) nseq=iseq 
      iseq=iseq+1
      do 227 jj=1,30
       if(names(jj).eq.rg(i)) then
       sequ(iseq)=short(jj)
       endif
  227 continue
      write(*,*)'i,re(i),ik,imx',i,re(i),ik,imax(ik)
      do 230 k=1,imax(ik)
      kk=i-1+k
      do 235 n=1,imax(ik)
      if(ty(kk).eq.pbtyp(ik,n).and.re(kk).eq.reold) then 
      j=i-1+ilam(ik,n)
      tyn(j)=ty(kk)
      rgn(j)=rg(kk)
      chn(j)=ch(kk)
      ren(j)=re(kk)
      xn(j)=xc(kk)
      yn(j)=yc(kk)
      zn(j)=zc(kk)
      iacn(j)=iatyp(ik,n)
      chan(j)=charge(ik,n)
      write(*,*)'detected',j,tyn(j),rgn(j),ik,imax(ik)
      endif
  235 continue
  230 continue
      endif
  225 continue
      do 237 i=1,natom
  237 write(*,25)at,i,tyn(i),rgn(i),chn(i),ren(i),xn(i),yn(i),zn(i)
      num=1
      chn(0)='n'
      reold=100010
      ireold=0
      rgold='    '
      do 260 i=1,natom
      if(ren(i).ne.reold) then
      reold=ren(i)
      ireold=ireold+1
      endif
      ren(i)=ireold
      if(tyn(i).eq.atyp(11)) then
c     ires=re(i)
      icarb=i
      endif
      if(chn(i).ne.chn(i-1).and.tyn(i).eq.'N   '
     1 .and.tyn(i+1).eq.'CA  ') then
      xc(num)=xn(i)
      yc(num)=yn(i)
      zc(num)=zn(i) 
      rg(num)=rgn(i) 
      ty(num)=tyn(i) 
      ch(num)=chn(i) 
      re(num)=ren(i)
      iac(num)=iacn(i)
      cha(num)=chan(i)
      num=num+1
      xc(num)=1.66*xn(i)-0.66*xn(i+1)
      yc(num)=1.66*yn(i)-0.66*yn(i+1)
      zc(num)=1.66*zn(i)-0.66*zn(i+1)
      rg(num)=rgn(i) 
      ty(num)=atyp(18) 
      ch(num)=chn(i) 
      re(num)=ren(i)
      iac(num)=32
      cha(num)=0.28
      num=num+1
c     write(*,*)'i,num',i,num,natom
      elseif(tyn(i).eq.atyp(5).and.tyn(i+1).eq.atyp(36).and.rgn(i)
     1 .ne.'PRO ') then
      xc(num)=xn(i)
      yc(num)=yn(i)
      zc(num)=zn(i)
      rg(num)=rgn(i)
      ty(num)=tyn(i)
      ch(num)=chn(i)
      re(num)=ren(i)
      iac(num)=iacn(i)
      cha(num)=chan(i)
      num=num+1
      xc(num)=0.75*(-xn(i+1)-xn(icarb)+2*xn(i))+xn(i)
      yc(num)=0.75*(-yn(i+1)-yn(icarb)+2*yn(i))+yn(i)
      zc(num)=0.75*(-zn(i+1)-zn(icarb)+2*zn(i))+zn(i)
      rg(num)=rgn(i) 
      ty(num)=atyp(18) 
      ch(num)=chn(i) 
      re(num)=ren(i)
      iac(num)=32
      cha(num)=0.280
      num=num+1
      else
      xc(num)=xn(i)
      yc(num)=yn(i)
      zc(num)=zn(i)
      rg(num)=rgn(i) 
      ty(num)=tyn(i) 
      ch(num)=chn(i) 
      re(num)=ren(i)
      iac(num)=iacn(i)
      cha(num)=chan(i)
      num=num+1 
      endif
  260 continue
      write(*,*)'natom,num',natom,num 
      natom=num-1
      reold=0
      open(48,file=ausgabe)
      do 100 i=1,natom
      ilp(i)=1
      ice(i)=1
c     if(ty(i).eq.atyp(18)) then
c     cha(i)=0.280
c     iac(i)=32
c     endif
      if(ch(i).ne.ch(i-1).and.i-1.gt.0) write(48,'(a4)') "TER "
c     write(48,27) at,i,tyn(i),rgn(i),chn(i),ren(i),xn(i),yn(i),zn(i),
c    1             iac(i),cha(i),ilp(i),ice(i) 
c     write(48,27) at,i,ty(i),rg(i),ch(i),re(i),xc(i),yc(i),zc(i),
      write(48,27) at,i,ty(i),rg(i),re(i),xc(i),yc(i),zc(i),
c    1             iac(i),cha(i),ilp(i),ice(i) 
     1             iac(i),cha(i),0,1.00 
c  27 format(a4,i7,2x,a4,a4,a1,i4,4x,3f8.3,i5,f8.3,2i2)
   27 format(a4,i7,2x,a4,a4,i5,4x,3f8.3,i5,f8.3,i2,f5.2)
  100 continue
      write(48,'(A3)') 'TER'
      write(48,'(A3)') 'END'
      close(48)
c check for disulfid bonds
      write(*,*)'possible cys-cys SG-SG bonds'
      ic=0
      do 300 i=1, natom-1
      if(ty(i).eq.'SG  ') then
       do 350 j=i+1,natom
        if(ty(j).eq.'SG  ') then
         dis=(xc(i)-xc(j))**2+(yc(i)-yc(j))**2+(zc(i)-zc(j))**2
         if(dis.lt.13.0) then
           ic=ic+1
           icys(ic,1)=i
           icys(ic,2)=j
           write(*,'(2i5)') i,j
         endif
        endif
  350  continue
      endif
  300 continue 
      open(44,file=ausgabe2)
      if(nseq.lt.iseq) then
      write(44,'(2i5)') nseq,ic
      write(44,'(600a)')(sequ(j),j=1,nseq)
      write(44,'(i5)') iseq-nseq
      write(44,'(600a)')(sequ(j),j=nseq+1,iseq)
      else
      write(44,'(2i5)') iseq,ic
      write(44,'(600a)')(sequ(j),j=1,iseq)
      endif
      do i=1,ic
      write(44,'(2i7)') icys(i,1), icys(i,2)
      enddo
      close(44)
      end
