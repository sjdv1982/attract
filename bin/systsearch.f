      program systsearch

      implicit none
      integer ntheta, nrot, ntrans, ktrans
      real*8 theta
      integer nphi
      dimension theta(20), nphi(20)
      real*8 ddx,ddy,ddz
      integer i,iii,jjj,kkk, count
      character*100 fmt
      real*8 phi,phii,ssi,rot, realcount
      
      real*8 pi, zwopi, ten, one
      pi=3.141592654d0
      zwopi=2.0d0*pi
      
      write (*,'(a11)') '#pivot auto'
      write (*,'(a25)') '#centered receptor: false'
      write (*,'(a23)') '#centered ligands: true'
      open(11,file='rotation.dat')
      read(11,*) ntheta,nrot
      do 593 i=1,ntheta
      read(11,*) theta(i),nphi(i)
      theta(i)=zwopi*theta(i)/360.0d0
c      write(*,*) theta(i),nphi(i)
  593 continue
      close(11)
      
      count = 0
      
      open(12,file='translate.dat')
      read(12,*) ntrans
      do 597 ktrans=1,ntrans
c read starting placements
      read(12,31) ddx,ddy,ddz
c generate various initial ligand protein orientations
      do 600 kkk=1,ntheta
      ssi=theta(kkk)
      phii=zwopi/nphi(kkk)
      do 620 jjj=1,nphi(kkk)
      phi=jjj*phii
      do 640 iii=1,nrot
      rot=iii*zwopi/nrot
      count = count + 1
      realcount  = count
      ten = 10.0
      one = 1.0
      write(fmt,'(a5,i1,a1)') '(a1,i', 
     1 int((log(realcount)+0.1)/log(ten)+one), ')'
      write (*,fmt) "#",count
      write (*,*) 0,0,0,0,0,0
      write (*,*) phi,ssi,rot,ddx,ddy,ddz
640   continue
620   continue
600   continue
597   continue

   31 format(30x,3f8.3)      
      end
