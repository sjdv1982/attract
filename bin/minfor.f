      subroutine minfor(maxatom,totmaxatom,maxres,
     1 maxlig, maxdof, maxmode, maxmolpair,
     2 cartstatehandle,ministatehandle,
     3 nhm, nlig, 
     4 ens, phi, ssi, rot, xa, ya, za, dlig, seed, label,
     5 gesa, energies, lablen)
c
c  variable metric minimizer (Harwell subroutine lib.  as in Jumna with modifications)
c     minimizes a single structure

      implicit none

c     Parameters
      integer cartstatehandle,ministatehandle
      integer maxlig,maxatom,totmaxatom,maxres,maxdof,maxmode,
     1 maxmolpair
      integer nlig, seed
      real*8 gesa, energies
      dimension energies(6)
      integer lablen
      character label
      dimension label(lablen)
      
      integer nhm
      dimension nhm(maxlig)
      integer ens
      dimension ens(maxlig)      
      real*8 phi, ssi, rot, dlig, xa, ya, za
      dimension phi(maxlig), ssi(maxlig), rot(maxlig)
      dimension dlig(maxmode, maxlig)
      dimension xa(maxlig), ya(maxlig), za(maxlig)

c     Local variables      
      integer i,ii,iii,n,j,jj,k,kk,ir,isfv,itr,nfun,np     
      integer itra, ieig, iori, fixre, gridmode,iscore,ivmax
      integer ju,jl,jb,nmodes
      integer iab
      real*8 c,acc,dff,dgb,enlig,f,fa,fmin,gl1,gl2,gmin,dga,xnull,w
      real*8 fb,dnorm
      
      real*8 h,g,ga,gb,xaa,xbb,d,delta,step,stepbd,steplb,stmin
      dimension h(maxdof*maxdof)
      dimension g(maxdof),ga(maxdof),gb(maxdof),w(maxdof)
      dimension xaa(maxdof), xbb(maxdof), d(maxdof)
      dimension delta(maxdof)

      call ministate_f_minfor(ministatehandle,
     1 iscore,ivmax,iori,itra,ieig,fixre,gridmode)
     
c     always calculate both forces and energies      
      iab = 1

c
c  all variables without lig-hm
c
      jb=3*iori*(nlig-fixre)+3*itra*(nlig-fixre)
c  all variables including lig-hm
      nmodes = 0
      do i=fixre, nlig
      nmodes = nmodes + nhm(i)
      enddo
      ju=jb+ieig*nmodes
c  only trans or ori
      jl=3*iori*(nlig-fixre)

      call ministate_calc_pairlist(ministatehandle,cartstatehandle)
     
      xnull=0.0d0
      nfun=0
      itr=0
      np=ju+1
      acc=0.000000001D0
c     set the hessian to a diagonal matrix 
      c = 1.0D0
      k=(ju*np)/2
      do i=1,k
       h(i)=xnull
      enddo
      k=1
      do i=1,ju
       h(k)=0.01d0*c
       k=k+np-i
      enddo
c     set some variables for the first iteration
      dff=xnull
      if (iscore.eq.1) then
        iori = 1
	itra = 1
      endif
      
      call energy(maxdof,maxmolpair,
     1 maxlig, maxatom,totmaxatom,maxmode,maxres,
     2 cartstatehandle,ministatehandle,
     3 iab,iori,itra,ieig,fixre,gridmode,
     4 ens,phi,ssi,rot,xa,ya,za,dlig,seed,
     5 gesa,energies,delta)
       
      if (iscore.eq.1) then
       write(*,*), 'Energy:', gesa 
       write(*,'(6f12.3)'), 
     1  energies(1),energies(2),energies(3),
     2  energies(4),energies(5),energies(6)
       j = jb
       do i=fixre,nlig-1
        ii = 3 * (i-fixre)
	iii = jl + 3 * (i-fixre)
	write(*,*),'Gradients:', delta(ii+1),delta(ii+2),delta(ii+3),
     1	 delta(iii+1),delta(iii+2),delta(iii+3)
        if ((ieig.eq.1).AND.nhm(i).gt.0) then
	 write(*,*), 'Mode gradients:'
	 do n=1,nhm(i)
          write(*,*),delta(j+1)
	 enddo
	endif		
	j = j + nhm(i)
       enddo      
       go to 256
      else if (iscore.eq.2) then
        call print_struc2(seed,label,gesa,energies,nlig,
     1  ens,phi,ssi,rot,xa,ya,za,nhm,dlig,lablen)	
      endif     
      
110   fa=gesa
      isfv=1
c store forces, Euler angle, position and ligand and receptor coordinates
c     write(*,*)'nlig,ju,jl,jb',nlig,ju,jl,jb
c     write(*,*)'delta',(delta(i),i=1,ju)
      do i=1,ju
       g(i)=-delta(i)
       ga(i)=g(i)
      enddo
c
c phi,ssi,rot for first molecule are fixed!
      if(iori.eq.1) then 
       do i=1+fixre,nlig
        ii=3*(i-fixre-1)
        xaa(ii+1)=phi(i)
        xaa(ii+2)=ssi(i)
        xaa(ii+3)=rot(i)
       enddo
      endif
      if(itra.eq.1) then
       do i=1+fixre,nlig
        ii=jl+3*(i-fixre-1)
        xaa(ii+1)=xa(i)
        xaa(ii+2)=ya(i)
        xaa(ii+3)=za(i)
       enddo
      endif
c if ligand flex is included store deformation factor in every mode in dlig(j)
 
      if(ieig.eq.1) then
       jj = 0
       do j=1,nlig
        do i=1,nhm(j)
         xaa(jb+jj+i)=dlig(i,j)
        enddo
        jj = jj + nhm(j)
       enddo
      endif
      
c     begin the iteration by giving the required printing
135   itr=itr+1
c     calculate the search direction of the iteration
      do i=1,ju
       d(i)=-ga(i)
      enddo
      call mc11e (h,ju,d,w,ju)
      do i=1,ju
c      write(*,*)'d(i)',i,d(i)
      enddo

c     calculate a lower bound on the step-length
c     and the initial directional derivative
      c=xnull
      dga=xnull
      do i=1,ju
      c=max(c,abs(d(i)))
      dga=dga+ga(i)*d(i)
c      write(*,*)'ga(i),d(i),dga,c',ga(i),d(i),dga,c
      enddo
c     test if the search direction is downhill
      if (dga.ge.xnull) go to 240
c     set the initial step-length of the line search
      stmin=xnull
      stepbd=xnull
      steplb=acc/c
      fmin=fa
      gmin=dga
      step=1.0d0
      if (dff.le.xnull) step=min(step,1.0d0/c)
      if (dff.gt.xnull) step=min(step,(dff+dff)/(-dga))
170   c=stmin+step
c     write(*,*)'dff,dga,stmin,step,c,nfun,ivmax',
c    1           dff,dga,stmin,step,c,nfun,ivmax
c     test whether func has been called ivmax times
      if (nfun.ge.ivmax) go to 220
      nfun=nfun+1
c     calculate another function value and gradient
      do i=1,ju
      g(i)=-delta(i)
c      write(*,*)'g(i)',g(i),d(i)
      enddo
c make an Euler rotation
      if(iori.eq.1) then
       do i=1+fixre,nlig
        ii=3*(i-fixre-1)
        phi(i)=xaa(ii+1)+c*d(ii+1)
        ssi(i)=xaa(ii+2)+c*d(ii+2)
        rot(i)=xaa(ii+3)+c*d(ii+3)
c        write(*,*)'new ii,c,phi,ssi,rot',i,ii,c,phi(i),ssi(i),rot(i)
       enddo
      endif
c make a move in HM direction and update x, y(1,i) and y(2,i) and dlig(j)
      if(ieig.eq.1) then
       kk = 0
       do k=1,nlig
        do i=1,nhm(k)
         dlig(i,k)=xaa(i+jb+kk)-d(i+jb+kk)*c
        enddo
        kk = kk + nhm(k)
       enddo
      endif
c make a translation of the ligand center
      if(itra.eq.1) then
       do i=1+fixre,nlig
        ii=jl+3*(i-fixre-1)
        xa(i)=xaa(ii+1)+c*d(ii+1)
        ya(i)=xaa(ii+2)+c*d(ii+2)
        za(i)=xaa(ii+3)+c*d(ii+3)
c       write(*,*)'new ii,c,xa ya za',i,ii,c,xa(i),ya(i),za(i),d(ii+1)
       enddo
      endif

      do i=1,3
c      write (*,*),'lig',i,phi(i),ssi(i),rot(i),xa(i),ya(i),za(i)
      enddo
      call energy(maxdof,maxmolpair,
     1 maxlig, maxatom,totmaxatom,maxmode,maxres,
     2 cartstatehandle,ministatehandle,
     3 iab,iori,itra,ieig,fixre,gridmode,
     4 ens,phi,ssi,rot,xa,ya,za,dlig,seed,
     5 fb,energies,delta)
      dnorm=xnull
      do i=1,ju
       dnorm=dnorm+delta(i)**2
      enddo
      dnorm=sqrt(dnorm)
c      write (*,*),'Energy2', fb,dnorm
      
      do i=1,ju
      gb(i)=-delta(i)
      enddo
      
c     store this function value if it is the smallest so far      
      if(iori.eq.1) then
       do i=1+fixre,nlig
        ii=3*(i-fixre-1)
        xbb(ii+1)=phi(i)
        xbb(ii+2)=ssi(i)
        xbb(ii+3)=rot(i)
       enddo
      endif
      if(itra.eq.1) then
       do i=1+fixre,nlig
        ii=jl+3*(i-fixre-1)
        xbb(ii+1)=xa(i)
        xbb(ii+2)=ya(i)
        xbb(ii+3)=za(i)
       enddo
      endif

      if(ieig.eq.1) then
       jj = 0
       do j=1,nlig
        do i=1,nhm(j)
         xbb(jb+jj+i)=dlig(i,j)
        enddo
        jj = jj + nhm(j)
       enddo
      endif

      isfv=min(2,isfv)
      if (fb.gt.gesa) go to 220
      if (fb.lt.gesa) go to 200
      gl1=xnull
      gl2=xnull
      do i=1,ju
       gl1=gl1+(g(i))**2
       gl2=gl2+(gb(i))**2
      enddo
      if (gl2.ge.gl1) go to 220
200   isfv=3
      gesa=fb
      if (iscore.eq.2) then
        call print_struc2(seed,label,gesa,energies,nlig,
     1  ens,phi,ssi,rot,xa,ya,za,nhm,dlig,lablen)	
      endif           
      do i=1,ju
       g(i)=gb(i)
      enddo
      
220   if(iori.eq.1) then
       do i=1+fixre,nlig
        ii=3*(i-fixre-1)
        phi(i)=xbb(ii+1)
        ssi(i)=xbb(ii+2)
        rot(i)=xbb(ii+3)
       enddo
      endif
      if(itra.eq.1) then
       do i=1+fixre,nlig
        ii=jl+3*(i-fixre-1)
        xa(i)=xbb(ii+1)
        ya(i)=xbb(ii+2)
        za(i)=xbb(ii+3)
       enddo
      endif 
      if (nfun.ge.ivmax) go to 250
      dgb=xnull
      do i=1,ju
       dgb=dgb+gb(i)*d(i)
      enddo
      
c     branch if we have found a new lower bound on the step-length
      if (fb-fa.le.0.1d0*c*dga) go to 280
c     finish the iteration if the current step is steplb
      if (step.gt.steplb) go to 270
240   if (isfv.ge.2) go to 110
c     at this stage the whole calculation is complete
250   if(nfun.lt.ivmax) then
      nfun=nfun+1
      endif
      
      call energy(maxdof,maxmolpair,
     1 maxlig, maxatom,totmaxatom,maxmode,maxres,
     2 cartstatehandle,ministatehandle,
     3 iab,iori,itra,ieig,fixre,gridmode,
     4 ens,phi,ssi,rot,xa,ya,za,dlig,seed,
     5 gesa,energies,delta)
      if (iscore.eq.2) then
        call print_struc2(seed,label,gesa,energies,nlig,
     1  ens,phi,ssi,rot,xa,ya,za,nhm,dlig,lablen)	
      endif     

      do 255 i=1,ju
255   g(i)=-delta(i)
      
256   call ministate_free_pairlist(ministatehandle)      

c      write (*,*),'Final energy:', gesa

      return
c     calculate a new step-length by cubic interpolation
270   stepbd=step
      c=gmin+dgb-3.0d0*(fb-fmin)/step
      c=gmin/(c+gmin-sqrt(c*c-gmin*dgb))
      step=step*max(0.1d0,c)
      go to 170
c     set the new bounds on the step-length
280   stepbd=stepbd-step
      stmin=c
      fmin=fb
      gmin=dgb
c     calculate a new step-length by extrapolation
      step=9.0d0*stmin
      if (stepbd.gt.xnull) step=0.5d0*stepbd
      c=dga+3.0d0*dgb-4.0d0*(fb-fa)/stmin
      if (c.gt.xnull) step=min(step,stmin*max(1.0d0,-dgb/c))
      if (dgb.lt.0.7d0*dga) go to 170
c     test for convergence of the iterations
      isfv=4-isfv
      if (stmin+step.le.steplb) go to 240
c     revise the second derivative matrix
      ir=-ju
      do i=1,ju
       xaa(i)=xbb(i)
       xbb(i)=ga(i)
       d(i)=gb(i)-ga(i)
       ga(i)=gb(i)
      enddo
      call mc11a(h,ju,xbb,1.0d0/dga,w,ir,1,xnull)
      ir=-ir
      call mc11a (h,ju,d,1.0d0/(stmin*(dgb-dga)),d,ir,0,xnull)
c     branch if the rank of the new matrix is deficient
      if (ir.lt.ju) go to 250
c     begin another iteration
      dff=fa-fb
      fa=fb
      go to 135
      end
