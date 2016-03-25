      subroutine minfor(
     1 smug, ivmax, x
     2 )
c
c  variable metric minimizer (Harwell subroutine lib.  as in Jumna with modifications)
c     minimizes a single structure


      use, intrinsic :: iso_c_binding, only: c_ptr
      implicit none

c     Parameters
      type(c_ptr) smug
      integer ivmax
      real*8 x
      dimension x(6)

c     Local variables
      real*8  gesa
      integer i,k,ir,isfv,itr,nfun,np,jn
      real*8 c,acc,dff,dgb,fa,fmin,gl1,gl2,gmin,dga,xnull,w
      real*8 fb
      
      real*8 h,g,ga,gb,xaa,xbb,d,step,stepbd,steplb,stmin
      dimension h(6*6)
      dimension g(6),ga(6),gb(6),w(6)
      dimension xaa(6), xbb(6), d(6)

      integer, parameter:: ERROR_UNIT = 0

      xnull=0.0d0
      dga=xnull
      jn = 6

      do i=1,jn
       ga(i)=xnull
       d(i)=xnull
       xaa(i)=x(i)
c       write(ERROR_UNIT,*)'xaa',i,xaa(i)
      enddo



      nfun=0
      itr=0
      np=jn+1
      acc=0.000000001D0
c     set the hessian to a diagonal matrix 
      c = 1.0D0
      k=(jn*np)/2
      do i=1,k
       h(i)=xnull
      enddo
      k=1
      do i=1,jn
       h(k)=0.01d0*c
       k=k+np-i
      enddo
c     set some variables for the first iteration
      dff=xnull
      call energy_for_fortran_to_call(smug, xaa, gesa, g)
      
110   fa=gesa
      isfv=1
c store forces, Euler angle, position and ligand and receptor coordinates
      do i=1,jn
       ga(i)=g(i)
      enddo
      
      
c     begin the iteration by giving the required printing
135   itr=itr+1
c     calculate the search direction of the iteration
      do i=1,jn
       d(i)=-ga(i)
      enddo
      call mc11e (h,jn,d,w,jn)
      do i=1,jn
c      write(ERROR_UNIT,*)'d(i)',i,d(i)
      enddo

c     calculate a lower bound on the step-length
c     and the initial directional derivative
      c=xnull
      dga=xnull
      do i=1,jn
      c=max(c,abs(d(i)))
      dga=dga+ga(i)*d(i)
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
c     test whether func has been called ivmax times
      if (nfun.ge.ivmax) go to 250
      nfun=nfun+1
c     calculate another function value and gradient
      
c     make an Euler rotation + tranlation of ligand center
      do i=1,jn
       xbb(i)=xaa(i)+c*d(i)
      enddo

      call energy_for_fortran_to_call(smug, xbb, fb, gb)
c     store this function value if it is the smallest so far      
      isfv=min(2,isfv)
      if (fb.gt.gesa) go to 220
      if (fb.lt.gesa) go to 200
      gl1=xnull
      gl2=xnull
      do i=1,jn
       gl1=gl1+(g(i))**2
       gl2=gl2+(gb(i))**2
      enddo
      if (gl2.ge.gl1) go to 220
200   isfv=3
      gesa=fb
      do i=1,jn
       x(i) = xbb(i)
       g(i)=gb(i)
      enddo
      
220   dgb=xnull
      do i=1,jn
       dgb=dgb+gb(i)*d(i)
      enddo
      
c     branch if we have found a new lower bound on the step-length
      if (fb-fa.le.0.1d0*c*dga) go to 280
c     finish the iteration if the current step is steplb
      if (step.gt.steplb) go to 270
240   if (isfv.ge.2) go to 110

c     at this stage the whole calculation is complete
250   return

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
      ir=-jn
      do i=1,jn
       xaa(i)=xbb(i)
       xbb(i)=ga(i)
       d(i)=gb(i)-ga(i)
       ga(i)=gb(i)
      enddo
      call mc11a(h,jn,xbb,1.0d0/dga,w,ir,1,xnull)
      ir=-ir
      call mc11a (h,jn,d,1.0d0/(stmin*(dgb-dga)),d,ir,0,xnull)
c     branch if the rank of the new matrix is deficient
      if (ir.lt.jn) go to 250
c     begin another iteration
      dff=fa-fb
      fa=fb
      go to 135
      end
