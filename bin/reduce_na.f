      program reduce
c
c  program to produce a reduced model of a protein structure
c  usage: $path/reduce structure.pdb
c  generates structurer.pdb as output pdb-file
c  author: Martin Zacharias, Jacobs University Bremen
c
      character*4 at,ty(15000),rg(15000),rgo,base(4,300)
      character*1 mark(15000)
      character*80 b
      character*40 name
      integer ilab(15000),iexi(15000),iexj(15000),ie(15000),
     1        ianf(15000)
      real xp(15000,3),xs(15000,3),xb(15000,3),sps(15000),
     1     psp(15000),psb(15000),bsp(900),x(15000,3),xi(3),
     2     xj(3),xk(3),xl(3),wpsp(20),wsps(20),wbsp(20),
     3     wpsb(20),wt1(20),wt2(20),wt3(20),wt4(20),
     4     wim(20)
c  get name of pdb-file
      call getarg(1,name)
      i=1
      k=2
      write(*,*) name
c     write(*,*) b
      kfi=index(name,' ')-5
      write(*,*) name(:kfi)//'.pdb'
      write(*,*) name(:kfi)//'r.pdb'
      open(42,file=name(:kfi)//'.pdb')
c  17 format(a,'.pdb')
      istart=0
      ianf(1)=1
c     open(42,file=b)
   40 read(42,20,end=30) b
      if(b(:4).eq.'ATOM') then
      read(b,26) at,ka,ty(i),rg(i),ie(i),mark(i),x(i,1),x(i,2),x(i,3)
      if(istart.eq.0) then
      iold=ie(i)
      istart=1
      endif
      if(ie(i).ne.iold) then
      ianf(k)=i
      k=k+1
      iold=ie(i)
      endif 
      i=i+1
      endif
      goto 40
   30 close(42)
      natom=i-1
      nres=k-1
      ianf(nres+1)=natom+1
      write(*,*) 'Number of atoms: ',natom
      write(*,*) 'Number of Residues: ', nres
      write(*,*) 'Start atom of each residue:'
      write(*,*) (ianf(i),i=1,nres+1)
      m=1
c     write(b,18) name
c  18 format(a4,'r.pdb')
      open(12,file=name(:kfi)//'r.pdb')
      do 100 i=1,nres
      xsum1=0.0
      ysum1=0.0
      zsum1=0.0
      xsum2=0.0
      ysum2=0.0
      zsum2=0.0
      xsum3=0.0
      ysum3=0.0
      zsum3=0.0
      xsum4=0.0
      ysum4=0.0
      zsum4=0.0
      xsum5=0.0
      ysum5=0.0
      zsum5=0.0
      kkk=1
      ima=0
      if(rg(ianf(i)).eq."  A ") then
      do 110 j=ianf(i),ianf(i+1)-1
      if(ty(j).eq."P   ") then 
      write(12,28)at,m,"GP1 ",rg(j),i,x(j,1),x(j,2),x(j,3),32,0.0,0,1.
      m=m+1
      o2counter=0
      endif
      if(ty(j).eq."C3' ".or.ty(j).eq."C4' ".or.ty(j).eq."C5' ") then
      xsum1=xsum1+x(j,1)
      ysum1=ysum1+x(j,2)
      zsum1=zsum1+x(j,3)
      endif 
      if(ty(j).eq."C3' ") then
      xsum1=0.3333*xsum1
      ysum1=0.3333*ysum1
      zsum1=0.3333*zsum1
      write(12,28)
     1 at,m,"GS1 ",rg(ianf(i)),i,xsum1,ysum1,zsum1,33,0.0,ima,1.
      m=m+1
      endif
      if(ty(j).eq."C2' ".or.ty(j).eq."C1' ".or.ty(j).eq."O4' ".or.ty(j)
     1.eq."O2' ") then
      xsum2=xsum2+x(j,1)
      ysum2=ysum2+x(j,2)
      zsum2=zsum2+x(j,3)
      o2counter=o2counter+1
      endif 
      if(ty(j).eq."C1' ") then
      xsum2=xsum2/o2counter
      ysum2=ysum2/o2counter
      zsum2=zsum2/o2counter
      write(12,28)
     1 at,m,"GS2 ",rg(ianf(i)),i,xsum2,ysum2,zsum2,34,0.0,ima,1.
      m=m+1
      o2counter=0
      endif
      if(ty(j).eq."N7  ") then 
      write(12,28)at,m,"GA1 ",rg(j),i,x(j,1),x(j,2),x(j,3),35,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."N3  ") then 
      write(12,28)at,m,"GA2 ",rg(j),i,x(j,1),x(j,2),x(j,3),36,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."N1  ") then 
      write(12,28)at,m,"GA3 ",rg(j),i,x(j,1),x(j,2),x(j,3),37,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."N6  ") then 
      write(12,28)at,m,"GA4 ",rg(j),i,x(j,1),x(j,2),x(j,3),38,0.0,0,1.
      m=m+1
      endif
  110 continue
      else if(rg(ianf(i)).eq."  G ") then
      do 112 j=ianf(i),ianf(i+1)-1
      if(ty(j).eq."P   ") then 
      write(12,28)at,m,"GP1 ",rg(j),i,x(j,1),x(j,2),x(j,3),32,0.0,0,1.
      m=m+1
      o2counter=0
      endif
      if(ty(j).eq."C3' ".or.ty(j).eq."C4' ".or.ty(j).eq."C5' ") then
      xsum1=xsum1+x(j,1)
      ysum1=ysum1+x(j,2)
      zsum1=zsum1+x(j,3)
      endif 
      if(ty(j).eq."C3' ") then
      xsum1=0.3333*xsum1
      ysum1=0.3333*ysum1
      zsum1=0.3333*zsum1
      write(12,28)
     1 at,m,"GS1 ",rg(ianf(i)),i,xsum1,ysum1,zsum1,33,0.0,ima,1.
      m=m+1
      endif
      if(ty(j).eq."C2' ".or.ty(j).eq."C1' ".or.ty(j).eq."O4' ".or.ty(j)
     1.eq."O2' ") then
      xsum2=xsum2+x(j,1)
      ysum2=ysum2+x(j,2)
      zsum2=zsum2+x(j,3)
      o2counter=o2counter+1
      endif 
      if(ty(j).eq."C1' ") then
      xsum2=xsum2/o2counter
      ysum2=ysum2/o2counter
      zsum2=zsum2/o2counter
      write(12,28)
     1 at,m,"GS2 ",rg(ianf(i)),i,xsum2,ysum2,zsum2,34,0.0,ima,1.
      m=m+1
      o2counter=0
      endif
      if(ty(j).eq."N7  ") then 
      write(12,28)at,m,"GG1 ",rg(j),i,x(j,1),x(j,2),x(j,3),39,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."N3  ") then 
      write(12,28)at,m,"GG2 ",rg(j),i,x(j,1),x(j,2),x(j,3),40,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."O6  ") then 
      write(12,28)at,m,"GG4 ",rg(j),i,x(j,1),x(j,2),x(j,3),42,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."N1  ".or.ty(j).eq."N2  ") then
      xsum3=xsum3+x(j,1)
      ysum3=ysum3+x(j,2)
      zsum3=zsum3+x(j,3)
      endif 
      if(ty(j).eq."N2 ") then
      xsum3=0.5*xsum3
      ysum3=0.5*ysum3
      zsum3=0.5*zsum3
      write(12,28)
     1 at,m,"GG3 ",rg(ianf(i)),i,xsum3,ysum3,zsum3,41,0.0,ima,1.
      m=m+1
      endif
  112 continue
      else if(rg(ianf(i)).eq."  C ") then
      do 114 j=ianf(i),ianf(i+1)-1
      if(ty(j).eq."P   ") then 
      write(12,28)at,m,"GP1 ",rg(j),i,x(j,1),x(j,2),x(j,3),32,0.0,0,1.
      m=m+1
      o2counter=0
      endif
      if(ty(j).eq."C3' ".or.ty(j).eq."C4' ".or.ty(j).eq."C5' ") then
      xsum1=xsum1+x(j,1)
      ysum1=ysum1+x(j,2)
      zsum1=zsum1+x(j,3)
      endif 
      if(ty(j).eq."C3' ") then
      xsum1=0.3333*xsum1
      ysum1=0.3333*ysum1
      zsum1=0.3333*zsum1
      write(12,28)
     1 at,m,"GS1 ",rg(ianf(i)),i,xsum1,ysum1,zsum1,33,0.0,ima,1.
      m=m+1
      endif
      if(ty(j).eq."C2' ".or.ty(j).eq."C1' ".or.ty(j).eq."O4' ".or.ty(j)
     1.eq."O2' ") then
      xsum2=xsum2+x(j,1)
      ysum2=ysum2+x(j,2)
      zsum2=zsum2+x(j,3)
      o2counter=o2counter+1
      endif 
      if(ty(j).eq."C1' ") then
      xsum2=xsum2/o2counter
      ysum2=ysum2/o2counter
      zsum2=zsum2/o2counter
      write(12,28)
     1 at,m,"GS2 ",rg(ianf(i)),i,xsum2,ysum2,zsum2,34,0.0,ima,1.
      m=m+1
      o2counter=0
      endif
      if(ty(j).eq."N4  ") then 
      write(12,28)at,m,"GC3 ",rg(j),i,x(j,1),x(j,2),x(j,3),45,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."C6  ".or.ty(j).eq."C5  ") then
      xsum3=xsum3+x(j,1)
      ysum3=ysum3+x(j,2)
      zsum3=zsum3+x(j,3)
      endif 
      if(ty(j).eq."C6 ") then
      xsum3=0.5*xsum3
      ysum3=0.5*ysum3
      zsum3=0.5*zsum3
      write(12,28)
     1 at,m,"GC1 ",rg(ianf(i)),i,xsum3,ysum3,zsum3,43,0.0,ima,1.
      m=m+1
      endif
      if(ty(j).eq."N3  ".or.ty(j).eq."O2  ") then
      xsum4=xsum4+x(j,1)
      ysum4=ysum4+x(j,2)
      zsum4=zsum4+x(j,3)
      endif 
      if(ty(j).eq."N3 ") then
      xsum4=0.5*xsum4
      ysum4=0.5*ysum4
      zsum4=0.5*zsum4
      write(12,28)
     1 at,m,"GC2 ",rg(ianf(i)),i,xsum4,ysum4,zsum4,44,0.0,ima,1.
      m=m+1
      endif
  114 continue
      else if(rg(ianf(i)).eq."  U ") then
      do 116 j=ianf(i),ianf(i+1)-1
      if(ty(j).eq."P   ") then 
      write(12,28)at,m,"GP1 ",rg(j),i,x(j,1),x(j,2),x(j,3),32,0.0,0,1.
      m=m+1
      o2counter=0
      endif
      if(ty(j).eq."C3' ".or.ty(j).eq."C4' ".or.ty(j).eq."C5' ") then
      xsum1=xsum1+x(j,1)
      ysum1=ysum1+x(j,2)
      zsum1=zsum1+x(j,3)
      endif 
      if(ty(j).eq."C3' ") then
      xsum1=0.3333*xsum1
      ysum1=0.3333*ysum1
      zsum1=0.3333*zsum1
      write(12,28)
     1 at,m,"GS1 ",rg(ianf(i)),i,xsum1,ysum1,zsum1,33,0.0,ima,1.
      m=m+1
      endif
      if(ty(j).eq."C2' ".or.ty(j).eq."C1' ".or.ty(j).eq."O4' ".or.ty(j)
     1.eq."O2' ") then
      xsum2=xsum2+x(j,1)
      ysum2=ysum2+x(j,2)
      zsum2=zsum2+x(j,3)
      o2counter=o2counter+1
      endif 
      if(ty(j).eq."C1' ") then
      xsum2=xsum2/o2counter
      ysum2=ysum2/o2counter
      zsum2=zsum2/o2counter
      write(12,28)
     1 at,m,"GS2 ",rg(ianf(i)),i,xsum2,ysum2,zsum2,34,0.0,ima,1.
      m=m+1
      o2counter=0
      endif
      if(ty(j).eq."C6  ".or.ty(j).eq."C5  ") then
      xsum3=xsum3+x(j,1)
      ysum3=ysum3+x(j,2)
      zsum3=zsum3+x(j,3)
      endif 
      if(ty(j).eq."C6 ") then
      xsum3=0.5*xsum3
      ysum3=0.5*ysum3
      zsum3=0.5*zsum3
      write(12,28)
     1 at,m,"GU1 ",rg(ianf(i)),i,xsum3,ysum3,zsum3,46,0.0,ima,1.
      m=m+1
      endif
      if(ty(j).eq."N3  ".or.ty(j).eq."O2  ") then
      xsum4=xsum4+x(j,1)
      ysum4=ysum4+x(j,2)
      zsum4=zsum4+x(j,3)
      endif 
      if(ty(j).eq."N3 ") then
      xsum4=0.5*xsum4
      ysum4=0.5*ysum4
      zsum4=0.5*zsum4
      write(12,28)
     1 at,m,"GU2 ",rg(ianf(i)),i,xsum4,ysum4,zsum4,47,0.0,ima,1.
      m=m+1
      endif
      if(ty(j).eq."N3  ".or.ty(j).eq."O4  ") then
      xsum5=xsum5+x(j,1)
      ysum5=ysum5+x(j,2)
      zsum5=zsum5+x(j,3)
      endif 
      if(ty(j).eq."O4 ") then
      xsum5=0.5*xsum5
      ysum5=0.5*ysum5
      zsum5=0.5*zsum5
      write(12,28)
     1 at,m,"GU3 ",rg(ianf(i)),i,xsum5,ysum5,zsum5,48,0.0,ima,1.
      m=m+1
      endif
  116 continue
      else if(rg(ianf(i)).eq."  T ") then
      do 118 j=ianf(i),ianf(i+1)-1
      if(ty(j).eq."P   ") then 
      write(12,28)at,m,"GP1 ",rg(j),i,x(j,1),x(j,2),x(j,3),32,0.0,0,1.
      m=m+1
      o2counter=0
      endif
      if(ty(j).eq."C3' ".or.ty(j).eq."C4' ".or.ty(j).eq."C5' ") then
      xsum1=xsum1+x(j,1)
      ysum1=ysum1+x(j,2)
      zsum1=zsum1+x(j,3)
      endif 
      if(ty(j).eq."C3' ") then
      xsum1=0.3333*xsum1
      ysum1=0.3333*ysum1
      zsum1=0.3333*zsum1
      write(12,28)
     1 at,m,"GS1 ",rg(ianf(i)),i,xsum1,ysum1,zsum1,33,0.0,ima,1.
      m=m+1
      endif
      if(ty(j).eq."C2' ".or.ty(j).eq."C1' ".or.ty(j).eq."O4' ".or.ty(j)
     1.eq."O2' ") then
      xsum2=xsum2+x(j,1)
      ysum2=ysum2+x(j,2)
      zsum2=zsum2+x(j,3)
      o2counter=o2counter+1
      endif 
      if(ty(j).eq."C1' ") then
      xsum2=xsum2/o2counter
      ysum2=ysum2/o2counter
      zsum2=zsum2/o2counter
      write(12,28)
     1 at,m,"GS2 ",rg(ianf(i)),i,xsum2,ysum2,zsum2,34,0.0,ima,1.
      m=m+1
      o2counter=0
      endif
      if(ty(j).eq."C6  ".or.ty(j).eq."C5  ".or.ty(j).eq."C7  ") then
      xsum3=xsum3+x(j,1)
      ysum3=ysum3+x(j,2)
      zsum3=zsum3+x(j,3)
      endif 
      if(ty(j).eq."C6 ") then
      xsum3=0.333*xsum3
      ysum3=0.333*ysum3
      zsum3=0.333*zsum3
      write(12,28)
     1 at,m,"GT1 ",rg(ianf(i)),i,xsum3,ysum3,zsum3,49,0.0,ima,1.
      m=m+1
      endif
      if(ty(j).eq."N3  ".or.ty(j).eq."O2  ") then
      xsum4=xsum4+x(j,1)
      ysum4=ysum4+x(j,2)
      zsum4=zsum4+x(j,3)
      endif 
      if(ty(j).eq."N3 ") then
      xsum4=0.5*xsum4
      ysum4=0.5*ysum4
      zsum4=0.5*zsum4
      write(12,28)
     1 at,m,"GT2 ",rg(ianf(i)),i,xsum4,ysum4,zsum4,50,0.0,ima,1.
      m=m+1
      endif
      if(ty(j).eq."N3  ".or.ty(j).eq."O4  ") then
      xsum5=xsum5+x(j,1)
      ysum5=ysum5+x(j,2)
      zsum5=zsum5+x(j,3)
      endif 
      if(ty(j).eq."O4 ") then
      xsum5=0.5*xsum5
      ysum5=0.5*ysum5
      zsum5=0.5*zsum5
      write(12,28)
     1 at,m,"GT3 ",rg(ianf(i)),i,xsum5,ysum5,zsum5,51,0.0,ima,1.
      m=m+1
      endif
  118 continue
      else if(rg(ianf(i)).eq."GLU ") then
      do 120 j=ianf(i),ianf(i+1)-1
      if(ty(j).eq."N   ") then 
      write(12,28)at,m,"N   ",rg(j),i,x(j,1),x(j,2),x(j,3),30,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."CA  ") then 
      write(12,28)at,m,"CA  ",rg(j),i,x(j,1),x(j,2),x(j,3),99,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."C   ") then 
      write(12,28)at,m,"C   ",rg(j),i,x(j,1),x(j,2),x(j,3),99,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."O   ") then 
      write(12,28)at,m,"O   ",rg(j),i,x(j,1),x(j,2),x(j,3),31,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."CG  ") then
       if(mark(j).ne." ") then
       xsum1=0.0
       ysum1=0.0
       zsum1=0.0
       ima=ima+1
       endif
      write(12,28)at,m,"CB  ",rg(j),i,x(j,1),x(j,2),x(j,3),10,0.0,ima,1.
      m=m+1
       endif 
      if(ty(j).eq."CD  ".or.ty(j).eq."OE1 ".or.ty(j).eq."OE2 ") then
      xsum1=xsum1+x(j,1)
      ysum1=ysum1+x(j,2)
      zsum1=zsum1+x(j,3)
      endif
      if(ty(j).eq."OE2 ") then
      xsum1=0.333*xsum1
      ysum1=0.333*ysum1
      zsum1=0.333*zsum1
      write(12,28)
     1 at,m,"CO1 ",rg(ianf(i)),i,xsum1,ysum1,zsum1,11,-1.0,ima,1.     
      m=m+1
      endif
  120 continue
      else if(rg(ianf(i)).eq."GLN ") then
      do 130 j=ianf(i),ianf(i+1)-1
      if(ty(j).eq."N   ") then 
      write(12,28)at,m,"N   ",rg(j),i,x(j,1),x(j,2),x(j,3),30,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."CA  ") then 
      write(12,28)at,m,"CA  ",rg(j),i,x(j,1),x(j,2),x(j,3),99,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."C   ") then 
      write(12,28)at,m,"C   ",rg(j),i,x(j,1),x(j,2),x(j,3),99,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."O   ") then 
      write(12,28)at,m,"O   ",rg(j),i,x(j,1),x(j,2),x(j,3),31,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."CG  ") then
       if(mark(j).ne." ") then
       xsum1=0.0
       ysum1=0.0
       zsum1=0.0
       ima=ima+1
       endif
      write(12,28)at,m,"CB  ",rg(j),i,x(j,1),x(j,2),x(j,3),8,0.0,ima,1.
      m=m+1
       endif 
      if(ty(j).eq."CD  ".or.ty(j).eq."OE1 ".or.ty(j).eq."NE2 ") then
      xsum1=xsum1+x(j,1)
      ysum1=ysum1+x(j,2)
      zsum1=zsum1+x(j,3)
      endif
      if(ty(j).eq."NE2 ") then
      xsum1=0.333*xsum1
      ysum1=0.333*ysum1
      zsum1=0.333*zsum1
      write(12,28)
     1at,m,"CN1 ",rg(ianf(i)),i,xsum1,ysum1,zsum1,9,0.0,ima,1.
      m=m+1
      endif
  130 continue
      else if(rg(ianf(i)).eq."LYS ") then
      do 140 j=ianf(i),ianf(i+1)-1
      if(ty(j).eq."N   ") then 
      write(12,28)at,m,"N   ",rg(j),i,x(j,1),x(j,2),x(j,3),30,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."CA  ") then 
      write(12,28)at,m,"CA  ",rg(j),i,x(j,1),x(j,2),x(j,3),99,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."C   ") then 
      write(12,28)at,m,"C   ",rg(j),i,x(j,1),x(j,2),x(j,3),99,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."O   ") then 
      write(12,28)at,m,"O   ",rg(j),i,x(j,1),x(j,2),x(j,3),31,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."CG  ") then
       if(mark(j).ne." ")  ima=ima+1
       write(12,28)
     1 at,m,"CB  ",rg(j),i,x(j,1),x(j,2),x(j,3),16,0.0,ima,1.
      m=m+1
      endif
      if(ty(j).eq."CE  ") then
      write(12,28)at,m,"CE  ",rg(j),i,x(j,1),x(j,2),x(j,3),17,1.0,ima,1.
      m=m+1
      endif
  140 continue
      else if(rg(ianf(i)).eq."TRP ") then
      icount=0
      do 150 j=ianf(i),ianf(i+1)-1
      if(ty(j).eq."N   ") then 
      write(12,28)at,m,"N   ",rg(j),i,x(j,1),x(j,2),x(j,3),30,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."CA  ") then 
      write(12,28)at,m,"CA  ",rg(j),i,x(j,1),x(j,2),x(j,3),99,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."C   ") then 
      write(12,28)at,m,"C   ",rg(j),i,x(j,1),x(j,2),x(j,3),99,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."O   ") then 
      write(12,28)at,m,"O   ",rg(j),i,x(j,1),x(j,2),x(j,3),31,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."CG  ") then
       if(mark(j).ne." ") then
       xsum1=0.0
       ysum1=0.0
       zsum1=0.0
       ima=ima+1
       endif
      write(12,28)at,m,"CG  ",rg(j),i,x(j,1),x(j,2),x(j,3),25,0.0,ima,1.
      m=m+1
       endif 
      if(ty(j).eq."CD2 ".or.ty(j).eq."CE2 ".or.ty(j).eq."CE3 ".or.
     1ty(j).eq."CH2 ".or.ty(j).eq."CZ3 ".or.ty(j).eq."CZ2 ") then
      icount=icount+1
      xsum1=xsum1+x(j,1)
      ysum1=ysum1+x(j,2)
      zsum1=zsum1+x(j,3)
      endif
      if(icount.ge.6) then
      xsum1=xsum1/6
      ysum1=ysum1/6
      zsum1=zsum1/6
      icount=0
      write(12,28)
     1 at,m,"CSE ",rg(ianf(i)),i,xsum1,ysum1,zsum1,26,0.0,ima,1.
      m=m+1
      endif
  150 continue
      else if(rg(ianf(i)).eq."MET ") then
      do 155 j=ianf(i),ianf(i+1)-1
      if(ty(j).eq."N   ") then 
      write(12,28)at,m,"N   ",rg(j),i,x(j,1),x(j,2),x(j,3),30,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."CA  ") then 
      write(12,28)at,m,"CA  ",rg(j),i,x(j,1),x(j,2),x(j,3),99,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."C   ") then 
      write(12,28)at,m,"C   ",rg(j),i,x(j,1),x(j,2),x(j,3),99,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."O   ") then 
      write(12,28)at,m,"O   ",rg(j),i,x(j,1),x(j,2),x(j,3),31,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."CB  ") then
       if(mark(j).ne." ") then
       xsum1=0.0
       ysum1=0.0
       zsum1=0.0
       xsum2=0.0
       ysum2=0.0
       zsum2=0.0
       ima=ima+1
       endif
      endif
      if(ty(j).eq."CB  ".or.ty(j).eq."CG  ") then
      xsum1=xsum1+x(j,1)
      ysum1=ysum1+x(j,2)
      zsum1=zsum1+x(j,3)
      endif
      if(ty(j).eq."SD  ".or.ty(j).eq."CE  ") then
      xsum2=xsum2+x(j,1)
      ysum2=ysum2+x(j,2)
      zsum2=zsum2+x(j,3)
      endif
      if(ty(j).eq."CE  ") then
      xsum1=xsum1/2
      ysum1=ysum1/2
      zsum1=zsum1/2
      xsum2=xsum2/2
      ysum2=ysum2/2
      zsum2=zsum2/2
      write(12,28)
     1at,m,"CSE ",rg(ianf(i)),i,xsum1,ysum1,zsum1,18,0.0,ima,1.
      m=m+1
      write(12,28)
     1at,m,"CSE ",rg(ianf(i)),i,xsum2,ysum2,zsum2,19,0.0,ima,1.
      m=m+1
      endif
  155 continue
      else if(rg(ianf(i)).eq."PHE ") then
      icount=0
      do 160 j=ianf(i),ianf(i+1)-1
      if(ty(j).eq."N   ") then 
      write(12,28)at,m,"N   ",rg(j),i,x(j,1),x(j,2),x(j,3),30,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."CA  ") then 
      write(12,28)at,m,"CA  ",rg(j),i,x(j,1),x(j,2),x(j,3),99,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."C   ") then 
      write(12,28)at,m,"C   ",rg(j),i,x(j,1),x(j,2),x(j,3),99,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."O   ") then 
      write(12,28)at,m,"O   ",rg(j),i,x(j,1),x(j,2),x(j,3),31,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."CB  ") then
       if(mark(j).ne." ") then
       xsum1=0.0
       ysum1=0.0
       zsum1=0.0
       xsum2=0.0
       ysum2=0.0
       zsum2=0.0
       ima=ima+1
       endif
      endif
      if(ty(j).eq."CB  ".or.ty(j).eq."CG  ") then
      xsum1=xsum1+x(j,1)
      ysum1=ysum1+x(j,2)
      zsum1=zsum1+x(j,3)
      endif
      if(ty(j).eq."CD1 ".or.ty(j).eq."CD2 ".or.ty(j).eq."CE1 "
     1.or.ty(j).eq."CE2 ".or.ty(j).eq."CZ  ") then
      icount=icount+1
      xsum2=xsum2+x(j,1)
      ysum2=ysum2+x(j,2)
      zsum2=zsum2+x(j,3)
      endif
      if(icount.ge.5) then
      xsum1=xsum1/2
      ysum1=ysum1/2
      zsum1=zsum1/2
      xsum2=xsum2/5
      ysum2=ysum2/5
      zsum2=zsum2/5
      icount=0
      write(12,28)
     1at,m,"CSE ",rg(ianf(i)),i,xsum1,ysum1,zsum1,20,0.0,ima,1.
      m=m+1
      write(12,28)
     1at,m,"CSE ",rg(ianf(i)),i,xsum2,ysum2,zsum2,21,0.0,ima,1.
      m=m+1
      endif
  160 continue
      else if(rg(ianf(i)).eq."TYR ") then
      icount=0
      do 170 j=ianf(i),ianf(i+1)-1
      if(ty(j).eq."N   ") then 
      write(12,28)at,m,"N   ",rg(j),i,x(j,1),x(j,2),x(j,3),30,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."CA  ") then 
      write(12,28)at,m,"CA  ",rg(j),i,x(j,1),x(j,2),x(j,3),99,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."C   ") then 
      write(12,28)at,m,"C   ",rg(j),i,x(j,1),x(j,2),x(j,3),99,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."O   ") then 
      write(12,28)at,m,"O   ",rg(j),i,x(j,1),x(j,2),x(j,3),31,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."CB  ") then
       if(mark(j).ne." ") then
       xsum1=0.0
       ysum1=0.0
       zsum1=0.0
       xsum2=0.0
       ysum2=0.0
       zsum2=0.0
       ima=ima+1
       endif
      endif      
      if(ty(j).eq."CB  ".or.ty(j).eq."CG  ") then
      xsum1=xsum1+x(j,1)
      ysum1=ysum1+x(j,2)
      zsum1=zsum1+x(j,3)
      endif
      if(ty(j).eq."CD1 ".or.ty(j).eq."CD2 ".or.ty(j).eq."CE1 "
     1.or.ty(j).eq."CE2 ".or.ty(j).eq."CZ  ".or.ty(j).eq."OH  ") then
      icount=icount+1
      xsum2=xsum2+x(j,1)
      ysum2=ysum2+x(j,2)
      zsum2=zsum2+x(j,3)
      endif
      if(icount.ge.6) then
      xsum1=xsum1/2
      ysum1=ysum1/2
      zsum1=zsum1/2
      xsum2=xsum2/6
      ysum2=ysum2/6
      zsum2=zsum2/6
      icount=0
      write(12,28)
     1at,m,"CSE ",rg(ianf(i)),i,xsum1,ysum1,zsum1,27,0.0,ima,1.
      m=m+1
      write(12,28)
     1at,m,"CSE ",rg(ianf(i)),i,xsum2,ysum2,zsum2,28,0.0,ima,1.
      m=m+1
      endif
  170 continue
      else if(rg(ianf(i)).eq."HIS ") then
      icount=0
      do 180 j=ianf(i),ianf(i+1)-1
      if(ty(j).eq."N   ") then 
      write(12,28)at,m,"N   ",rg(j),i,x(j,1),x(j,2),x(j,3),30,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."CA  ") then 
      write(12,28)at,m,"CA  ",rg(j),i,x(j,1),x(j,2),x(j,3),99,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."C   ") then 
      write(12,28)at,m,"C   ",rg(j),i,x(j,1),x(j,2),x(j,3),99,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."O   ") then 
      write(12,28)at,m,"O   ",rg(j),i,x(j,1),x(j,2),x(j,3),31,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."CB  ") then
       if(mark(j).ne." ") then
       xsum1=0.0
       ysum1=0.0
       zsum1=0.0
       xsum2=0.0
       ysum2=0.0
       zsum2=0.0
       ima=ima+1
       endif
      endif
      if(ty(j).eq."CB  ".or.ty(j).eq."CG  ") then
      xsum1=xsum1+x(j,1)
      ysum1=ysum1+x(j,2)
      zsum1=zsum1+x(j,3)
      endif
      if(ty(j).eq."ND1 ".or.ty(j).eq."CD2 ".or.ty(j).eq."NE2 "
     1.or.ty(j).eq."CE1 ") then
      icount=icount+1
      xsum2=xsum2+x(j,1)
      ysum2=ysum2+x(j,2)
      zsum2=zsum2+x(j,3)
      endif
      if(icount.ge.4) then
      xsum1=xsum1/2
      ysum1=ysum1/2
      zsum1=zsum1/2
      xsum2=xsum2/4
      ysum2=ysum2/4
      zsum2=zsum2/4
      icount=0
      write(12,28)
     1at,m,"CSE ",rg(ianf(i)),i,xsum1,ysum1,zsum1,12,0.0,ima,1.
      m=m+1
      write(12,28)
     1at,m,"CSE ",rg(ianf(i)),i,xsum2,ysum2,zsum2,13,0.0,ima,1.
      m=m+1
      endif
  180 continue
      else if(rg(ianf(i)).eq."GLY ") then
      do 190  j=ianf(i),ianf(i+1)-1
      if(ty(j).eq."N   ") then 
      write(12,28)at,m,"N   ",rg(j),i,x(j,1),x(j,2),x(j,3),30,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."CA  ") then 
      write(12,28)at,m,"CA  ",rg(j),i,x(j,1),x(j,2),x(j,3),1,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."C   ") then 
      write(12,28)at,m,"C   ",rg(j),i,x(j,1),x(j,2),x(j,3),99,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."O   ") then 
      write(12,28)at,m,"O   ",rg(j),i,x(j,1),x(j,2),x(j,3),31,0.0,0,1.
      m=m+1
      endif
  190 continue
      else
      do 200 j=ianf(i),ianf(i+1)-1
      if(ty(j).eq."N   ") then 
      write(12,28)at,m,"N   ",rg(j),i,x(j,1),x(j,2),x(j,3),30,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."CA  ") then 
      write(12,28)at,m,"CA  ",rg(j),i,x(j,1),x(j,2),x(j,3),99,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."C   ") then 
      write(12,28)at,m,"C   ",rg(j),i,x(j,1),x(j,2),x(j,3),99,0.0,0,1.
      m=m+1
      endif
      if(ty(j).eq."O   ") then 
      write(12,28)at,m,"O   ",rg(j),i,x(j,1),x(j,2),x(j,3),31,0.0,0,1.
      m=m+1
      endif
      if (ty(j).ne."N   ".and.ty(j).ne."O   ".and.ty(j).ne."C   ") then
      xsum1=xsum1+x(j,1)
      ysum1=ysum1+x(j,2)
      zsum1=zsum1+x(j,3)
      kkk=kkk+1
      endif
  200 continue
      xsum1=xsum1/(kkk-1)
      ysum1=ysum1/(kkk-1)
      zsum1=zsum1/(kkk-1)
      if(rg(ianf(i)).eq."ALA ") then
      cha=0.0
      ittt=2
      else if (rg(ianf(i)).eq."ASN ") then
      cha=0.0
      ittt=5
      else if (rg(ianf(i)).eq."ASP ") then
      cha=-1.0
      ittt=6
      else if (rg(ianf(i)).eq."CYS ") then
      cha=0.0
      ittt=7
      else if (rg(ianf(i)).eq."ILE ") then
      cha=0.0
      ittt=14
      else if (rg(ianf(i)).eq."LEU ") then
      cha=0.0
      ittt=15
      else if (rg(ianf(i)).eq."PRO ") then
      cha=0.0
      ittt=22
      else if (rg(ianf(i)).eq."SER ") then
      cha=0.0
      ittt=23
      else if (rg(ianf(i)).eq."THR ") then
      cha=0.0
      ittt=24
      else if (rg(ianf(i)).eq."VAL ") then
      cha=0.0
      ittt=29
      endif
      write(12,28)at,m,"CSE ",rg(ianf(i)),i,xsum1,ysum1,zsum1,
     1 ittt,cha,0,1.
      m=m+1
      endif
  100 continue
   20 format(a80)
   26 format(a4,i7,2x,a4,a4,1x,i4,a1,3x,3f8.3)
   27 format(a12,i5,a4,4f13.6)
   28 format(a4,i7,2x,a4,a4,1x,i4,4x,3f8.3,i5,f8.3,i2,f5.2)
      end
c
