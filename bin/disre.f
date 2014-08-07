      subroutine disre(maxlig,cartstatehandle,ministatehandle,
     1 iab,iori,itra,fixre,xa,ya,za,
     2 locrests, has_locrests,
     3 delta,erest)
c  adds a restraining contribution due to attraction between 
c centers of the ligand
c
      implicit none

c     Parameters
      integer maxlig,cartstatehandle,ministatehandle
      integer iab,iori,itra,fixre
      real*8 delta(maxlig),xa(maxlig),ya(maxlig),za(maxlig)
      real*8 erest
      real*8 locrests
      dimension locrests(3,maxlig)
      integer has_locrests      
      dimension has_locrests(maxlig)
      
c     Handle variables
      integer nlig
      real*8 pivot(maxlig,3)
      integer gravity
      real*8 rstk
      pointer (ptr_pivot,pivot)
            
c     Local variables
      integer n,nn,i,jl
      real*8 et,et2,xd,yd,zd,xf,yf,zf,k
      
      erest = 0.0d0
      call cartstate_f_disre(cartstatehandle,nlig,ptr_pivot)      
      call ministate_f_disre(ministatehandle,gravity,rstk)
           
      jl=3*iori*(nlig-fixre)
      
      do 100 n=1+fixre,nlig
      i = jl + 3 * (n-fixre-1)
      xf = 0.0d0
      yf = 0.0d0
      zf = 0.0d0
      
      if (has_locrests(n).gt.0) then
      xd=xa(n)+pivot(n,1)-locrests(1,n)
      yd=ya(n)+pivot(n,2)-locrests(2,n)
      zd=za(n)+pivot(n,3)-locrests(3,n)
      et=(xd**2+yd**2+zd**2)

      erest=erest+rstk*et
      xf=xf-2.0d0*rstk*xd
      yf=yf-2.0d0*rstk*yd
      zf=zf-2.0d0*rstk*zd      

      endif

      if (gravity.eq.3.or.gravity.eq.4.or.gravity.eq.5) then
c     gravity=3 => to all other centers
c     gravity=4 => away from all other centers
c     gravity=5 => away from all other centers
            
      do nn=1,nlig           
      xd=xa(n)+pivot(n,1)-xa(nn)-pivot(nn,1)
      yd=ya(n)+pivot(n,2)-ya(nn)-pivot(nn,2)
      zd=za(n)+pivot(n,3)-za(nn)-pivot(nn,3)
      et=(xd**2+yd**2+zd**2)

      et2 = et
      k = rstk
      if (gravity.eq.4) then
c     gravity=4 => away from all other centers
        if (et.lt.(200*200)) then
        et2 = 200*200 - et
        k = -rstk
        else
        et2 = 0
        k = 0
        endif
      endif
      if (gravity.eq.5) then
c     gravity=5 => somewhat away from all other centers
        if (et.lt.(20*20)) then
        et2 = 20*20 - et
        k = -10*rstk
        else
        et2 = 0
        k = 0
        endif
      endif
      
      erest=erest+rstk*et2
      xf=xf-2.0d0*k*xd
      yf=yf-2.0d0*k*yd
      zf=zf-2.0d0*k*zd      
      
      end do
      
      endif
      
      if (gravity.eq.1.or.gravity.eq.2.or.gravity.eq.4.
     1 or.gravity.eq.5) then

      if (gravity.eq.1.or.gravity.eq.4.or.gravity.eq.5) then
      if (gravity.eq.1) k = rstk
c     gravity=1 => to global origin      
      if (gravity.eq.4.or.gravity.eq.5) k = 0.5 * rstk
c     gravity=4,5 => to global origin, two times as weak
      xd=xa(n)+pivot(n,1)
      yd=ya(n)+pivot(n,2)
      zd=za(n)+pivot(n,3)
      else if (gravity.eq.2) then
c     gravity=2 => to receptor origin
      k = rstk
      xd=xa(n)+pivot(n,1)-xa(1)-pivot(1,1)
      yd=ya(n)+pivot(n,2)-ya(1)-pivot(1,2)
      zd=za(n)+pivot(n,3)-za(1)-pivot(1,3)
      endif
      
      et=(xd**2+yd**2+zd**2)

c     Fourth order restraints, maybe a bit too steep...
c      erest=erest+k*et*et
c      xf=xf-4.0d0*k*et*xd
c      yf=yf-4.0d0*k*et*yd
c      zf=zf-4.0d0*k*et*zd       

c     Harmonic restraints
      erest=erest+k*et
      xf=xf-2.0d0*k*xd
      yf=yf-2.0d0*k*yd
      zf=zf-2.0d0*k*zd       

      endif

      if ((iab.eq.1).and.(n.gt.1.or.fixre.eq.0)) then
      
      delta(i+1) = delta(i+1) + xf
      delta(i+2) = delta(i+2) + yf
      delta(i+3) = delta(i+3) + zf
      end if
      
100   continue      
      
      return
      end
