      program rmsca
c
c  program to superimpose two structures with respect to CA-atoms
c  usage: $path/rmsca structure1.pdb structure2.pdb
c  author: Martin Zacharias, Jacobs University Bremen
c 
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      integer iflag
      dimension cor(2,19000,3),aa(6),cc(3,3),r(3,3),
     & u(3,3),d(3),corn(19000,3),dn(3),ires(8)
      character*12 file(100)
      character*80,buff,name1,name2,name3,name4,name5
      character*4 atom,typ(19000),re(19000),sucr(1),ttyp(19000),
     1            rre(19000)
c     data sucr/'N   ','CA  ','C   ','O   ','H    ','CB  '/
      data sucr/'CA  '/ 
      data ires/1,2,3,4,8,9,10,11/
      integer ie(19000),nu(19000),ire(19000),nnu(19000),iie(19000),
     1        jre(19000)
      nsurc=11
      call getarg(1,name1)
      call getarg(2,name2)
      iflag=1
      open(unit=42,file=name1,status='old')
      do m=1,3
      d(m)=0.
      enddo
      i=1
      j=1
40    read(42,20,end=30) buff
      if(buff(:4).eq.'ATOM') then
      read(buff,25) atom,ie(i),typ(i),re(i),nu(i),
     1              (cor(1,i,m),m=1,3) 
      do 35 k=1,1
      if(sucr(k).eq.typ(i)) then
      ire(j)=i
      do m=1,3
      d(m)=d(m)+cor(1,i,m)
      enddo
      j=j+1
      endif
   35 continue
      i=i+1
      endif
      goto 40
   30 natm=i-1
      matm=j-1
      close(42)
      do m=1,3
      d(m)=d(m)/float(matm)
      enddo
      do k=1,natm
      do m=1,3
      cor(1,k,m)=cor(1,k,m)-d(m)
      enddo
      enddo
      do m=1,3
      dn(m)=d(m)
      d(m)=0.
      enddo
      open(unit=42,file=name2,status='old')
      i=1
      j=1
140   read(42,20,end=130) buff
      if(buff(:4).eq.'ATOM') then
      read(buff,25) atom,iie(i),ttyp(i),rre(i),nnu(i),
     1    (cor(2,i,m),m=1,3)
      do 135 k=1,1
      if(sucr(k).eq.ttyp(i)) then
      jre(j)=i
      j=j+1
      do m=1,3
      d(m)=d(m)+cor(2,i,m)
      enddo
      endif
  135 continue
      i=i+1
      endif
      goto 140
  130 natm2=i-1
      close(42)
c     write(*,*)'natom' ,natm,natm2,matm,j-1
c     do 145 k=1,matm
c     write(*,*) k,ire(k),typ(ire(k)),jre(k),ttyp(jre(k))
c 145 continue
c     if(natm.ne.natm2.or.matm.ne.j-1) goto 1001
      do m=1,3
      d(m)=d(m)/float(matm)
      enddo
      do k=1,natm2
      do m=1,3
      cor(2,k,m)=cor(2,k,m)-d(m)
      enddo
      enddo
      do m=1,3
      do n=1,3
      r(m,n)=0.0
      do k=1,matm
      r(m,n)=r(m,n)+cor(1,ire(k),m)*cor(2,jre(k),n)
      enddo
      enddo
      enddo
      l=0
      do m=1,3
      do n=1,m
      l=l+1
      aa(l)=0.
      do k=1,3
      aa(l)=aa(l)+r(k,m)*r(k,n)
      enddo
      enddo
      enddo
      call eigen(aa,cc,3)
      d(1)=sqrt(1./aa(1))
      d(2)=sqrt(1./aa(3))
      d(3)=sqrt(1./aa(6))  
      do m=1,3
      do n=1,3
      u(m,n)=0.
      do k=1,3
      do l=1,3
      u(m,n)=u(m,n)+r(m,k)*cc(k,l)*cc(n,l)*d(l)
      enddo
      enddo
      enddo
      enddo
      rms=0.
      klm=1
      do k=1,natm2
      do m=1,3
      dd=cor(1,k,m)
      corn(k,m)=0.
      do n=1,3
      corn(k,m)=corn(k,m)+u(m,n)*cor(2,k,n)
      dd=dd-u(m,n)*cor(2,k,n)
      enddo
      enddo
      enddo
      dd=0.0
      open(13,file='rmsd.dat')
      do k=1,matm
      ddd=(corn(jre(k),1)-cor(1,ire(k),1))**2 
     1     +(corn(jre(k),2)-cor(1,ire(k),2))**2
     2     +(corn(jre(k),3)-cor(1,ire(k),3))**2
      dd=dd+ddd
      write(13,'(2f8.3)') float(k),sqrt(ddd)
      enddo
      close(13)
      rms=sqrt(dd/float(matm))
      write(*,'(f9.4,2i5)') rms
      if(iflag.eq.1) then
      open(12,file='out.pdb')
      write(12,'(a11)')'MODEL     1'
      do i=1,natm
      write(12,25) atom,ie(i),typ(i),re(i),nu(i),
     1             (cor(1,i,m)+dn(m),m=1,3)
      end do
      write(12,'(a3)')'TER'
      write(12,'(a3)')'ENDMDL'
      write(12,'(a11)')'MODEL     2'
      do i=1,natm2
      write(12,25) atom,iie(i),ttyp(i),rre(i),nnu(i),
     1             (corn(i,m)+dn(m),m=1,3)
      end do
      write(12,'(a3)')'TER'
      write(12,'(a3)')'END'
      close(12)
      endif
200   format(1x,a12,' / ',a12,'    ...rms= ',f8.3,2i5)
   20 format(a80)
   25 format(a4,i7,2x,a4,a4,2x,i3,4x,3f8.3)
   28 format(a4,i7,3x,a4,a3,2x,i3,4x,3f8.3)
   26 format(30x,3f8.3)
      goto 2001
1001  write(*,*)'number of atoms in the two structures differ',
     1           natm,natm2,matm,j-1
2001  end
c
      subroutine eigen(a,r,n)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension a(1),r(1)
      range=1.0D-12
      mv=0
      iq=-n
      do 20 j=1,n
      iq=iq+n
      do 20 i=1,n
      ij=iq+i
      r(ij)=0.0
      if(i-j) 20,15,20
15    r(ij)=1.0
20    continue
      anorm=0.0
      do 35 i=1,n
      do 35 j=i,n
      if(i-j) 30,35,30
30    ia=i+(j*j-j)/2
      anorm=anorm+a(ia)*a(ia)
35    continue
      if(anorm) 165,165,40
40    anorm=sqrt(2.0*anorm)
      anrmx=anorm*range/float(n)
      ind=0
      thr=anorm
45    thr=thr/float(n)
50    l=1
55    m=l+1
60    mq=(m*m-m)/2
      lq=(l*l-l)/2
      lm=l+mq
62    if(abs(a(lm))-thr) 130,65,65
65    ind=1
      ll=l+lq
      mm=m+mq
      x=0.5*(a(ll)-a(mm))
68    y=-a(lm)/sqrt(a(lm)*a(lm)+x*x)
      if(x) 70,75,75
70    y=-y
75    sinx=y/sqrt(2.0*(1.0+(sqrt(1.0-y*y))))
      sinx2=sinx*sinx
78    cosx=sqrt(1.0-sinx2)
      cosx2=cosx*cosx
      sincs=sinx*cosx
      ilq=n*(l-1)
      imq=n*(m-1)
      do 125 i=1,n
      iq=(i*i-i)/2
      if(i-l) 80,115,80
80    if(i-m) 85,115,90
85    im=i+mq
      go to 95
90    im=m+iq
95    if(i-l) 100,105,105
100   il=i+lq
      go to 110
105   il=l+iq
110   x=a(il)*cosx-a(im)*sinx
      a(im)=a(il)*sinx+a(im)*cosx
      a(il)=x
115   if(mv-1) 120,125,120
120   ilr=ilq+i
      imr=imq+i
      x=r(ilr)*cosx-r(imr)*sinx
      r(imr)=r(ilr)*sinx+r(imr)*cosx
      r(ilr)=x
125   continue
      x=2.0*a(lm)*sincs
      y=a(ll)*cosx2+a(mm)*sinx2-x
      x=a(ll)*sinx2+a(mm)*cosx2+x
      a(lm)=(a(ll)-a(mm))*sincs+a(lm)*(cosx2-sinx2)
      a(ll)=y
      a(mm)=x
130   if(m-n) 135,140,135
135   m=m+1
      go to 60
140   if(l-(n-1)) 145,150,145
145   l=l+1
      go to 55
150   if(ind-1) 160,155,160
155   ind=0
      go to 50
160   if(thr-anrmx) 165,165,45
165   iq=-n
      do 185 i=1,n
      iq=iq+n
      ll=i+(i*i-i)/2
      jq=n*(i-2)
      do 185 j=i,n
      jq=jq+n
      mm=j+(j*j-j)/2
      if(a(ll)-a(mm)) 170,185,185
170   x=a(ll)
      a(ll)=a(mm)
      a(mm)=x
      if(mv-1) 175,185,175
175   do 180 k=1,n
      ilr=iq+k
      imr=jq+k
      x=r(ilr)
      r(ilr)=r(imr)
180   r(imr)=x
185   continue
      return
      end
