      program rotam 
c
c program to generate rotameric side chain states for a set of
c residues (in list.dat)a, rotamers are in rota.dat
c usage: $path/rotam structure.pdb rcut (cutoff radius to eliminate rotamers
c overlap with other atoms)
c author: Martin Zacharias, Jacobs University Bremen
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension  sta(9,6,2)
      common /coord/ x(6000),y(6000),z(6000),xn(6000),yn(6000),zn(6000)
      character*4 at,ty(2000),rg(2000),rgr(2000)
      character*80 name,name2,name3,b
      character*1 g(20)
      integer*4 ire(2000),irer(2000),in(500),im(500),ia(500),iu(500),
     1     il(500),ic(20),icm(20),ipair(500)
      data g/'A','B','C','D','E','F','G','H','I','J','K','L',
     1     'M','N','O','P','Q','R','S','T'/ 
      ang=57.29577951d0
      pi=3.141592654d0
      nrota=7
      call getarg(1,name)   
      call getarg(2,name2)   
      read(name2,*) rcut
      open(41,file='rota.dat')
      do i=1,3
      read(41,'(a4,2i4,12f7.1)') at,ii,jj,
     1 (sta(i,1,j),j=1,2),(sta(i,2,j),j=1,2),(sta(i,3,j),j=1,2),
     2 (sta(i,4,j),j=1,2),(sta(i,5,j),j=1,2),(sta(i,6,j),j=1,2)
      do m=1,6
      do j=1,3
      sta(i,m,j)=sta(i,m,j)/ang
      enddo
      enddo
      enddo
      do i=4,9
      read(41,'(a4,2i4,3f7.1)') at,ii,jj,
     1 sta(i,1,1),sta(i,2,1),sta(i,3,1)
      do m=1,3
      sta(i,m,1)=sta(i,m,1)/ang
      enddo
      enddo
      close(41)
      open(42,file=name)
      i=0
      ih=1
      iold=0
      iza=0
 40   read(42,20,end=30) b
 20   format(a80)
      if(b(:4).eq.'ATOM') then
      read(b,25)at,k,ty(i+1),rg(i+1),ire(i+1),x(i+1),y(i+1),z(i+1)
c     write(*,25)at,k,ty(i+1),rg(i+1),ire(i+1),x(i+1),y(i+1),z(i+1)
         i=i+1
         if(ire(i).ne.iold) then
            iold=ire(i)
            iza=iza+1
            rgr(iza)=rg(i)
         endif
         irer(i)=iza
         if(ty(i).eq.'OXT ') nlim=i
         at=rg(i)
         if(ty(i).eq.'CA  ') im(iza)=i
         if(ty(i).eq.'CB  ') then
            in(iza)=i
            il(iza)=i+1
         endif
         if(ty(i).eq.'N   ') iu(iza-1)=i-1
      endif 
      goto 40
 30   close(42)
      nlim=i
      natom=i
      nres=iza
      write(*,*)'natom,nre,nlim', natom,nres,nlim
      open(10,file='list.dat')
      read(10,*) nmax
      do 45 i=1,nmax
       read(10,*) ic(i),icm(i)
       write(*,*) ic(i),icm(i),rgr(ic(i))
 45   continue 
      close(10)
c check torsions winkel um CA-CB Bindung (N-CA-CB-CG)
      do 80 i=1,nmax
      ik=ic(i)
      call torcheck(im(ik)-1,im(ik),in(ik),il(ik),wink)
      write(*,*)'chi1',rgr(ik),ik,ang*wink
      if(rgr(ik).eq."LYS ".or.rgr(ik).eq."ARG ".or.rgr(ik).eq."MET ") 
     1 then
      call torcheck(im(ik),in(ik),il(ik),il(ik)+1,wink)
      write(*,*)'chi2',rgr(ik),ik,ang*wink
c regularize the torsion to 180:
      if(wink.le.0.0d0) then
      pdel=-pi-wink
      else
      pdel=pi-wink
      endif
      write(*,*)'pdel', pdel
      call torsion2(0,in(ik),il(ik),il(ik)+1,il(ik)+2,-pdel)
      call torcheck(in(ik),il(ik),il(ik)+1,il(ik)+2,wink)
      write(*,*)'chi3',rgr(ik),ik,ang*wink
      if(wink.le.0.0d0) then
      pdel=-pi-wink
      else
      pdel=pi-wink
      endif
      write(*,*)'pdel', pdel
      call torsion2(0,il(ik),il(ik)+1,iu(ik),il(ik)+2,-pdel)
      call torcheck(in(ik),il(ik),il(ik)+1,il(ik)+2,wink)
      write(*,*)'chi3n',rgr(ik),ik,ang*wink
      call torcheck(il(ik),il(ik)+1,il(ik)+2,il(ik)+3,wink)
      write(*,*)'chi4',rgr(ik),ik,ang*wink
      if(wink.le.0.0d0) then
      pdel=-pi-wink
      else
      pdel=pi-wink
      endif
      write(*,*)'pdel', pdel
      call torsion2(0,il(ik)+1,il(ik)+2,iu(ik),il(ik)+3,-pdel)
      call torcheck(il(ik),il(ik)+1,il(ik)+2,il(ik)+3,wink)
      write(*,*)'chi4n',rgr(ik),ik,ang*wink
      endif
   80 continue
c
c the file reg.pdb contains the structure with regularized torsions and
c residue numbering starting from 1
c
      open(11,file='reg.pdb')
      do 90 i=1,natom
      xn(i)=x(i)
      yn(i)=y(i)
      zn(i)=z(i)
      write(11,25)'ATOM',i,ty(i),rg(i),irer(i),x(i),y(i),z(i) 
   90 continue
      close(11)
      kpdb=0
      open(10,file='rota.pdb')
      k=1
      do 150 i=1,nlim
         if(irer(i).eq.ic(k)+1.and.ty(i).eq.'N   ') then
         ik=ic(k)
c     
c     for the selected residue a pairlist of neighbors needs to be generated
c     use all residues within 10.0 A of the CB atom
c     
         call pairl(ik,in(ik),irer,ipair,npair,natom)
         write(*,*)'pairs',in(ik),(ipair(kkj),kkj=1,npair)
c
c   this is to determine reference torsion angles
c
      if(rgr(ik).eq."LYS ".or.rgr(ik).eq."ARG ".or.
     1   rgr(ik).eq."MET ") then
         call torcheck(im(ik)-1,im(ik),in(ik),il(ik),win1)
         write(*,*)'chi1',rgr(ik),ik,ang*win1
         call torcheck(im(ik),in(ik),il(ik),il(ik)+1,win2)
         write(*,*)'chi2',rgr(ik),ik,ang*win2
c generate rotamers and check overlap
         do 160 ii=1,6 
         pdel1=sta(1,ii,1)-win1
         pdel2=sta(1,ii,2)-win2
         write(*,*)'pdel',ii,pdel1,pdel2
         call torsion2(0,im(ik),in(ik),iu(ik),il(ik),-pdel1)
         call torsion2(0,in(ik),il(ik),iu(ik),il(ik)+1,-pdel2)
         call overlap (il(ik)+1,iu(ik),ipair,npair,natom,ifl,rcut)
c        write(*,*)'here we are',ifl,natom,npair
c        write(*,*)'ifl0',i,k,irer(i),ic(k),ifl
         if(ifl.ne.1) then
         do 110 j=in(ik),iu(ik)
         kpdb=kpdb+1
         write(10,26)'ATOM',kpdb,
c        write(10,26)'ATOM',j+(ii-1)*(iu(ik)-in(ik)),
     1   ty(j),rg(j),g(ii),irer(j),g(ii),x(j),y(j),z(j)
 110     continue
         write(10,'(a3)')'TER'
         endif 
         do 170 j=in(ik),iu(ik)
         x(j)=xn(j)
         y(j)=yn(j)
         z(j)=zn(j)
 170     continue
 160     continue
         k=k+1
       else if(rgr(ik).eq."GLU ".or.rgr(ik).eq."GLN ".or.
     1         rgr(ik).eq."HIS ".or.rgr(ik).eq."TRP ".or.
     2         rgr(ik).eq."PHE ".or.rgr(ik).eq."TYR ") then
         call torcheck(im(ik)-1,im(ik),in(ik),il(ik),win1)
         write(*,*)'chi1',rgr(ik),ik,ang*win1
         do 165 ii=1,3
         pdel1=sta(4,ii,1)-win1
         write(*,*)'pdel',ii,pdel1
         call torsion2(0,im(ik),in(ik),iu(ik),il(ik),-pdel1)
         call overlap (il(ik)+1,iu(ik),ipair,npair,natom,ifl,rcut)
c        write(*,*)'ifl',i,k,irer(i),ic(k),ifl
         if(ifl.ne.1) then
         do 185 j=in(ik),iu(ik)
         kpdb=kpdb+1
         write(10,26)'ATOM',kpdb,
c        write(10,26)'ATOM',j+(ii-1)*(iu(ik)-in(ik)),
     1   ty(j),rg(j),g(ii),irer(j),g(ii),x(j),y(j),z(j)
 185     continue
         write(10,'(a3)')'TER'
         endif
         do 175 j=in(ik),iu(ik)
         x(j)=xn(j)
         y(j)=yn(j)
         z(j)=zn(j)
 175     continue
 165     continue
         k=k+1
         endif
c  reset coordinates for next rotamer
         endif
c do not write out the original copy if copies are written
c        write(*,*)'ifl2',i,k,irer(i),ic(k),ifl
         if(irer(i).ne.ic(k).or.ty(i).eq."N   ".
     1   or.ty(i).eq."CA  ".or.ty(i).eq."C   ".or.ty(i).eq."O   ") then
         kpdb=kpdb+1
         write(10,25)'ATOM',kpdb,
c        write(10,25)'ATOM',i+icm(k-1)*(iu(ik)-in(ik)),
     1   ty(i),rg(i),irer(i),xn(i),yn(i),zn(i)
         endif
 150  continue
      write(10,'(a3)') 'TER'
      close(10)
 25   format(a4,i7,2x,a4,a4,2x,i3,4x,3f8.3)
 26   format(a4,i7,2x,a4,a4,a1,1x,i3,a1,3x,3f8.3)
      end
c
      subroutine torcheck(ii,jj,kk,ll,winkel)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      common /coord/ x(6000),y(6000),z(6000),xn(6000),yn(6000),zn(6000)
      dimension dij(3),djk(3),dkl(3),dik(3),djl(3),rik(3),rjl(3),ril(3)
      dij(1)=x(jj)-x(ii)
      djk(1)=x(jj)-x(kk)
      dkl(1)=x(ll)-x(kk)
      dik(1)=x(kk)-x(ii)
      djl(1)=x(ll)-x(jj)
      dij(2)=y(jj)-y(ii)
      djk(2)=y(jj)-y(kk)
      dkl(2)=y(ll)-y(kk)
      dik(2)=y(kk)-y(ii)
      djl(2)=y(ll)-y(jj)
      dij(3)=z(jj)-z(ii)
      djk(3)=z(jj)-z(kk)
      dkl(3)=z(ll)-z(kk)
      dik(3)=z(kk)-z(ii)
      djl(3)=z(ll)-z(jj)
      call cross(dij,djk,rik)
      call cross(djk,dkl,rjl)
      rikn=0.0d0
      rjln=0.0d0
      do 130 j=1,3
      rikn=rikn+rik(j)**2 
      rjln=rjln+rjl(j)**2 
  130 continue
c     write(*,*) rikn,rjln
      rikn=1.0d0/dsqrt(rikn)
      rjln=1.0d0/dsqrt(rjln)
      cijkl=0.0d0
      do 140 j=1,3
      rik(j)=rik(j)*rikn
      rjl(j)=rjl(j)*rjln
      cijkl=cijkl+rik(j)*rjl(j) 
  140 continue
      if(cijkl.ge.1.0d0) cijkl=0.99999999999d0
      if(cijkl.le.-1.0d0) cijkl=-0.99999999999d0
c     write(*,*)'cijkl',cijkl
      call cross(rik,rjl,ril)
      skjil=djk(1)*ril(1)+djk(2)*ril(2)+djk(3)*ril(3)
c     write(*,*)'cijkl',cijkl,skjil
      if(skjil.le.0.0d0) then
      winkel=-1.0d0*dacos(cijkl)
      else
      winkel=1.0d0*dacos(cijkl)
      endif
c     skjil=skjil/sqrt(djk(1)**2+djk(2)**2+djk(3)**2)
c     write(*,*)'cijkl',cijkl,skjil
c     winkel=skjil*dacos(cijkl)
c     write(*,*)'winkel',winkel
      if(winkel.ge.0.0d0.and.winkel.le.0.0000005d0) winkel=0.0000005d0
      if(winkel.lt.0.0d0.and.winkel.ge.-0.0000005d0)winkel=-0.0000005d0
c     write(*,*)'winkel',winkel
      return
      end
c
      subroutine cross(r1,r2,cr)
      real*8 r1(3),r2(3),cr(3)
      cr(1)=r1(3)*r2(2)-r1(2)*r2(3)
      cr(2)=r1(1)*r2(3)-r1(3)*r2(1)
      cr(3)=r1(2)*r2(1)-r1(1)*r2(2)
      return 
      end
c
      subroutine overlap (ilow,iup,ipair,npair,natom,ifl,rcut)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      common /coord/ x(6000),y(6000),z(6000),xn(6000),yn(6000),zn(6000)
      dimension ipair(500)
      ifl=0
      rlim=rcut*rcut
      do j=ilow,iup
c     write(*,*) j,ilow,iup,npair,x(j),y(j),z(j)
      do i=1,npair
      k=ipair(i)
c     write(*,*)'k',k
      if((x(j)-xn(k))**2+(y(j)-yn(k))**2+(z(j)-zn(k))**2.
     1   le.rlim) then
      ifl=1
      return
      endif
      enddo
      enddo
      return
      end
c
      subroutine pairl (ik,ibet,irer,ipair,npair,natom)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      common /coord/ x(6000),y(6000),z(6000),xn(6000),yn(6000),zn(6000)
      integer*4 ipair(500),irer(1900)
      j=0
      do i=1,natom
      if((x(ibet)-x(i))**2+(y(ibet)-y(i))**2+(z(ibet)-z(i))**2.
     1   le.70.0.and.irer(i).ne.ik) then
      j=j+1
      ipair(j)=i
      endif
      enddo
      npair=j
      return
      end
c
      subroutine torsion (le,hh,kk,iend,ianf,pdel)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      common /coord/ x(6000),y(6000),z(6000),xn(6000),yn(6000),zn(6000)
      integer le,hh,kk,ianf,iend
c     torsion angle around bond between atoms hh and kk is
c     rotated by an angle alpha. zcs=cos(alpha), zss=sin(alpha)
c     le is the number of atoms in the system, ianf is the 
c     first atom which moves upon rotating around the bond,
c     iend is the last atom to be moved. All atoms in between
c     these two numbers will be moved so the atoms must be
c     appropriately ordered. 
      zcs=dcos(pdel)
      zss=-dsin(pdel)
      rhk1=x(kk)-x(hh)
      rhk2=y(kk)-y(hh)
      rhk3=z(kk)-z(hh)
      rhkb=rhk1**2+rhk2**2+rhk3**2
      rhkb=sqrt(rhkb)
      rhk1=rhk1/rhkb
      rhk2=rhk2/rhkb
      rhk3=rhk3/rhkb
c
c     go through all atoms between ianf and iend
c
      do 100,i=ianf,iend
      rd1=x(i)-x(kk)
      rd2=y(i)-y(kk)
      rd3=z(i)-z(kk)
c
c     projection on vector rhk
c
      rl3=rhk1*rd1+rhk2*rd2+rhk3*rd3
      rl1=rl3*rhk1
      rl2=rl3*rhk2
      rl3=rl3*rhk3
c
c     define dr
c
      dr1=rd1-rl1
      dr2=rd2-rl2
      dr3=rd3-rl3
c
c     define rn
c
      rn1=rhk2*dr3-dr2*rhk3
      rn2=rhk3*dr1-dr3*rhk1
      rn3=rhk1*dr2-dr1*rhk2
c
c     calculate drn
c
      drn1=zcs*dr1+zss*rn1
      drn2=zcs*dr2+zss*rn2
      drn3=zcs*dr3+zss*rn3
c
c     new position
c
      xn(i)=x(kk)+rl1+drn1
      yn(i)=y(kk)+rl2+drn2
      zn(i)=z(kk)+rl3+drn3
  100 continue
      return
      end
c
      subroutine torsion2 (le,hh,kk,iend,ianf,pdel)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      common /coord/ x(6000),y(6000),z(6000),xn(6000),yn(6000),zn(6000)
      integer le,hh,kk,ianf,iend
c     torsion angle around bond between atoms hh and kk is
c     rotated by an angle alpha. zcs=cos(alpha), zss=sin(alpha)
c     le is the number of atoms in the system, ianf is the 
c     first atom which moves upon rotating around the bond,
c     iend is the last atom to be moved. All atoms in between
c     these two numbers will be moved so the atoms must be
c     appropriately ordered. 
      zcs=dcos(pdel)
      zss=-dsin(pdel)
      rhk1=x(kk)-x(hh)
      rhk2=y(kk)-y(hh)
      rhk3=z(kk)-z(hh)
      rhkb=rhk1**2+rhk2**2+rhk3**2
      rhkb=sqrt(rhkb)
      rhk1=rhk1/rhkb
      rhk2=rhk2/rhkb
      rhk3=rhk3/rhkb
c
c     go through all atoms between ianf and iend
c
      do 100,i=ianf,iend
      rd1=x(i)-x(kk)
      rd2=y(i)-y(kk)
      rd3=z(i)-z(kk)
c
c     projection on vector rhk
c
      rl3=rhk1*rd1+rhk2*rd2+rhk3*rd3
      rl1=rl3*rhk1
      rl2=rl3*rhk2
      rl3=rl3*rhk3
c
c     define dr
c
      dr1=rd1-rl1
      dr2=rd2-rl2
      dr3=rd3-rl3
c
c     define rn
c
      rn1=rhk2*dr3-dr2*rhk3
      rn2=rhk3*dr1-dr3*rhk1
      rn3=rhk1*dr2-dr1*rhk2
c
c     calculate drn
c
      drn1=zcs*dr1+zss*rn1
      drn2=zcs*dr2+zss*rn2
      drn3=zcs*dr3+zss*rn3
c
c     new position
c
      x(i)=x(kk)+rl1+drn1
      y(i)=y(kk)+rl2+drn2
      z(i)=z(kk)+rl3+drn3
  100 continue
      return
      end
