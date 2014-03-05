      subroutine monte(cartstatehandle,ministatehandle,
     1 nhm, nihm, nlig,
     2 ens, phi, ssi, rot, xa, ya, za, morph, dlig,
     3 locrests, has_locrests,
     4 seed, label,
     5 gesa, energies, lablen)
c
c  variable metric minimizer (Harwell subroutine lib.  as in Jumna with modifications)
c     minimizes a single structure

      implicit none

c     Parameters
      integer cartstatehandle,ministatehandle
      include 'max.fin'
      integer nlig, seed
      real*8 locrests
      dimension locrests(3,maxlig)
      integer has_locrests      
      dimension has_locrests(maxlig)
      real *8 gesa, energies
      dimension energies(6)
      integer lablen
      character label
      dimension label(lablen)
      
      integer nhm
      dimension nhm(maxlig)
      integer nihm
      dimension nihm(maxlig)
      integer ens
      dimension ens(maxlig)
      real*8 phi, ssi, rot, dlig, xa, ya, za, morph
      dimension phi(maxlig), ssi(maxlig), rot(maxlig)
      dimension dlig(maxmode+maxindexmode, maxlig)
      dimension xa(maxlig), ya(maxlig), za(maxlig)
      dimension morph(maxlig)      

c     Local variables      
      real*8 enew, energies0
      real*8 rrot1,rrot2,rrot3,rrot4,sphi,sssi,srot
      dimension energies0(6)
c     integer dseed,i,ii,j,jj,k,kk,itr,nfun
      integer i,ii,j,jj,k,kk,itr,nfun
      integer itra, ieig, iindex, iori, fixre, gridmode,iscore,ivmax
      integer ju,ju0,jl,jb,nmodes,nimodes, jn, jn0
      integer iab,ijk,iaccept
      real*8 xnull
      real*8 scalecenter,scalemode,ensprob,scalerot,rr
      real*8 rotmat,randrot,newrot,sum 
      real*8 xaa,delta,deltamorph,bol,pi,mctemp
      integer ensaa
      real*8 dseed
      dimension xaa(maxdof)
      dimension ensaa(maxlig)      
      dimension delta(maxdof), deltamorph(maxlig)
      dimension rr(maxdof),randrot(0:9),rotmat(0:9),newrot(0:9)
      integer nrens
      dimension nrens(maxlig)
      pointer(ptr_nrens,nrens)
      real*8 neomorph
      integer, parameter :: ERROR_UNIT = 0
      pi=3.141592654d0

      do i=1, maxlig      
      ensaa(i) = 0
      enddo
      
      call ministate_f_monte(ministatehandle,
     1 iscore,ivmax,iori,itra,ieig,iindex,fixre,gridmode,mctemp,
     2 scalerot,scalecenter,scalemode,ensprob)
     
c     always calculate only energies
      iab = 0

c
c  all variables without lig-hm
c
      jb=3*iori*(nlig-fixre)+3*itra*(nlig-fixre)
c  all variables including lig-hm
      nmodes = 0
      nimodes = 0
      do 5 i=fixre, nlig
      nmodes = nmodes + nhm(i)
      nimodes = nimodes + nihm(i)
    5 continue
      ju=jb+ieig*nmodes
      jn = ju + iindex*nimodes
      ju0=ju
      jn0 = jn
      do i=1,nlig
      if (morph(i).ge.0) then
      jn = jn + 1
      endif
      enddo      
      
c  only trans or ori
      jl=3*iori*(nlig-fixre)

      call ministate_calc_pairlist(ministatehandle,cartstatehandle)
      call cartstate_get_nrens(cartstatehandle, ptr_nrens)      
     
      xnull=0.0d0
c     dseed=seed
      dseed=12345
      scalerot=pi*scalerot
      nfun=0
      itr=0
      if (iscore.eq.1) then
        iori = 1
        itra = 1
      endif
c intial energy evaluation      
c
      call energy(cartstatehandle,ministatehandle,
     1 iab,iori,itra,ieig,iindex,fixre,gridmode,
     2 ens,phi,ssi,rot,xa,ya,za,morph,dlig,
     3 locrests, has_locrests, seed,
     4 gesa,energies,delta,deltamorph)
       
      if (iscore.eq.2) then
        call print_struc2(seed,label,gesa,energies,nlig,
     1  ens,phi,ssi,rot,xa,ya,za,locrests,morph,
     2  nhm,nihm,dlig,has_locrests,lablen)
      endif     
      
c   start Monte Carlo
      iaccept=1
      do 4000 ijk=1,ivmax 
c store old Euler angle, position and ligand and receptor coordinates
c
c phi,ssi,rot for first molecule are fixed!
      if(iaccept.eq.1) then
      ensaa(i)=ens(i)      
      if(iori.eq.1) then 
      do 118 i=1+fixre,nlig
      ii=3*(i-fixre-1)
      xaa(ii+1)=phi(i)
      xaa(ii+2)=ssi(i)
      xaa(ii+3)=rot(i)
  118 continue
      endif
      if(itra.eq.1) then
      do 122 i=1+fixre,nlig
      ii=jl+3*(i-fixre-1)
      xaa(ii+1)=xa(i)
      xaa(ii+2)=ya(i)
      xaa(ii+3)=za(i)
  122 continue
      endif
      jj = jn0
      do i=1,nlig
       if (morph(i).ge.0) then
        xaa(jj+1) = morph(i)
        jj = jj + 1
       endif      
      enddo
      
c if ligand flex is included store deformation factor in every mode in dlig(j)
 
      if(ieig.eq.1) then
      jj = 0
      do 130 j=1,nlig
      do 131 i=1,nhm(j)
      xaa(jb+jj+i)=dlig(i,j)
  131 continue
      jj = jj + nhm(j)
  130 continue
      endif

      if(iindex.eq.1) then
      jj = 0
      do 140 j=1,nlig
      do 141 i=1,nihm(j)
      xaa(ju0+jj+i)= dlig(ju0+i,j)
  141 continue
      jj = jj + nihm(j)
  140 continue
      endif

      endif
c old Cartesians are not stored!
c generate a total of ju random numbers
c     write(*,*)'random rr(1),rr(2)..', rr(1),rr(2),rr(3)
c make an Euler rotation
      do i=1,nlig
        if (nrens(i).gt.0.and.morph(i).lt.0) then
	  call GGUBS(dseed,2,rr)
	  if (rr(1).lt.ensprob) then
c	    ens(i) = int(rr(2)*nrens(i))+1
            call enstrans(cartstatehandle,i-1,ens(i),rr(2), ens(i))
	  endif
        endif
      enddo
      if(iori.eq.1) then
        do 190 i=1+fixre,nlig
        ii=3*(i-fixre-1)
c        write(*,*)'before',phi(i),ssi(i),rot(i)
c       call crand(dseed,5,rr)
        call GGUBS(dseed,5,rr)
c       write(*,*)'rand',rr(1),rr(2),rr(3),rr(4),rr(5)
c       dseed = int(10000*rr(5))
c       write(*,*)'dseed', dseed
        rr(1)=rr(1)-0.5d0
        rr(2)=rr(2)-0.5d0
        rr(3)=rr(3)-0.5d0
        sum=1d0/sqrt(rr(1)**2+rr(2)**2+rr(3)**2)
        rr(1)=sum*rr(1)
        rr(2)=sum*rr(2)
        rr(3)=sum*rr(3)
        rr(4)= scalerot*(rr(4)-0.5d0) 
        rrot1=rr(1) 
        rrot2=rr(2) 
        rrot3=rr(3) 
        rrot4=rr(4) 
        sphi=phi(i)
        sssi=ssi(i)
        srot=rot(i)
c       write(*,*)'phi(i),ssi(i),rot(i)',i,phi(i),ssi(i),rot(i)
        call axisrot(rr,randrot)
        call euler2rotmat(phi(i),ssi(i),rot(i),rotmat)
        call matmult(rotmat,randrot,newrot)
        phi(i) = atan2(newrot(5),newrot(2))
        ssi(i) = acos(newrot(8))
        rot(i) = atan2(-newrot(7),-newrot(6))       
        if (abs(newrot(8)) >= 0.9999) then 
c         gimbal lock        
          phi(i) = 0.0d0
          if (abs(newrot(0)) >= 0.9999) then
            if (newrot(0) * newrot(8) < 0) then
              rot(i) = pi
            else
              rot(i) = 0.0d0      
            endif     
            if (newrot(8) < 0) then
              ssi(i) = pi
            else
              ssi(i) = 0.0d0
            endif               
            ssi(i) = 0.0d0
            rot(i) = 0.0d0
          else          
            if (newrot(8) < 0) then
              ssi(i) = pi
              rot(i) = -acos(-newrot(0))
            else 
              ssi(i) = 0.0d0
              rot(i) = acos(newrot(0))
            endif
          endif
          if (newrot(1) < 0) then
            rot(i) = -rot(i)
          endif          
        endif      
c       write(*,*)'new ii,c,phi,ssi,rot',i,ii,c,phi(i),ssi(i),rot(i)
  190   continue
      endif
c make a move in HM direction and update x, y(1,i) and y(2,i) and dlig(j)
c     call crand(dseed,ju+1,rr) 
      call GGUBS(dseed,jn+1,rr)
c     dseed = int(10000*rr(ju+1))
      if(ieig.eq.1) then
      kk = 0
      do 1180 k=1,nlig
      do 1200 i=1,nhm(k)
      dlig(i,k)=xaa(i+jb+kk)+scalemode*(rr(i+jb+kk)-0.5d0)
 1200 continue
      kk = kk + nhm(k)
 1180 continue
      endif
      if(iindex.eq.1) then
      kk = 0
      do 1280 k=1,nlig
      do 1300 i=1,nihm(k)
      dlig(ju0+i,k)=xaa(i+ju0+kk)+scalemode*(rr(i+ju0+kk)-0.5d0)
 1300 continue
      kk = kk+ nihm(k)
 1280 continue
      endif
c make a translation of the ligand center
      if(itra.eq.1) then
      do 1220 i=1+fixre,nlig
      ii=jl+3*(i-fixre-1)
      xa(i)=xa(i)+scalecenter*(0.5d0-rr(ii+1))
      ya(i)=ya(i)+scalecenter*(0.5d0-rr(ii+2))
      za(i)=za(i)+scalecenter*(0.5d0-rr(ii+3))
c     write(*,*)'trans-step',i,ii,
c    1 0.5d0-rr(ii+1),0.5d0-rr(ii+2),0.5d0-rr(ii+3)
c    1 rr(ii+1),rr(ii+2),rr(ii+3),xaa(ii+1),xaa(ii+2),xaa(ii+3)
 1220 continue
      endif
      
      jj = jn0
      do i=1,nlig
       if (morph(i).ge.0) then
        neomorph = morph(i)+scalemode*(0.5d0-rr(ii+1))
        if (neomorph.lt.0) neomorph = 0
        if (neomorph.gt.nrens(i)-1.001) neomorph = nrens(i)-1.001
        morph(i) = neomorph
        jj = jj + 1
       endif      
      enddo
      
      call energy(cartstatehandle,ministatehandle,
     1 iab,iori,itra,ieig,iindex,fixre,gridmode,
     2 ens,phi,ssi,rot,xa,ya,za,morph,dlig,
     3 locrests, has_locrests, seed,
     4 enew,energies0,delta,deltamorph)
c  new energy 
c      write (ERROR_UNIT,*),'Energy2', enew 
      bol=enew-gesa
      if (mctemp.eq.0) then
      bol=sign(1.0d0,-bol)
      else
      bol=exp(-bol/mctemp)
      endif
c      write(*,*)'exp(bol)',enew,gesa,enew-gesa,bol
c     call crand(dseed,2,rr)
      call GGUBS(dseed,2,rr)
c     dseed = int(10000*rr(2))
      if(bol.gt.rr(1)) then
c      write(ERROR_UNIT,*)'accept the step', bol, rr(1)
c     write(*,*)
c    1 'rrot1,rrot2,rrot3,rrot4,sphi,phi(i),sssi,ssi(i),srot,rot(i)',
c    2 rrot1,rrot2,rrot3,rrot4,sphi,phi(2),sssi,ssi(2),srot,rot(2)
      gesa=enew
      energies(:)=energies0(:)
      iaccept=1
      if (iscore.eq.2) then
        call print_struc2(seed,label,gesa,energies,nlig,
     1  ens,phi,ssi,rot,xa,ya,za,locrests,morph,
     2  nhm,nihm,dlig,has_locrests,lablen)
      endif           
c overwrite old xaa variables, see above
      else
c do not overwrite xaa variables
c      write(ERROR_UNIT,*)' step rejected'
      iaccept=0
      ens(i)=ensaa(i)
      if(iori.eq.1) then
      do 1118 i=1+fixre,nlig
      ii=3*(i-fixre-1)
      phi(i)=xaa(ii+1)
      ssi(i)=xaa(ii+2)
      rot(i)=xaa(ii+3)
 1118 continue
      endif
      if(itra.eq.1) then
      do 1122 i=1+fixre,nlig
      ii=jl+3*(i-fixre-1)
      xa(i)=xaa(ii+1)
      ya(i)=xaa(ii+2)
      za(i)=xaa(ii+3)
 1122 continue
      endif
      
c if ligand flex is included store deformation factor in every mode in dlig(j)

      if(ieig.eq.1) then
      jj = 0
      do 230 j=1,nlig
      do 231 i=1,nhm(j)
      dlig(i,j)=xaa(jb+jj+i)
  231 continue
      jj = jj + nhm(j)
  230 continue
      endif

      if(iindex.eq.1) then
      jj = 0
      do 240 j=1,nlig
      do 241 i=1,nihm(j)
      dlig(nhm(j)+i,j)=xaa(ju0+jj+i)
  241 continue
      jj = jj + nihm(j)
  240 continue
      endif
      endif

      jj = jn0
      do i=1,nlig
       if (morph(i).ge.0) then
        morph(i) = xaa(jj+1)
        jj = jj + 1
       endif      
      enddo


 4000 continue

c     Clean up
      call ministate_free_pairlist(ministatehandle)      
      end
