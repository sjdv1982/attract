      subroutine disre(maxlig,cartstatehandle,ministatehandle,
     1 iab,iori,itra,fixre,xa,ya,za,delta,erest)
c  adds a restraining contribution due to attraction between 
c centers of the ligand
c
      implicit none

c     Parameters
      integer maxlig,cartstatehandle,ministatehandle
      integer iab,iori,itra,fixre
      real*8 delta(maxlig),xa(maxlig),ya(maxlig),za(maxlig)
      real*8 erest
      
c     Handle variables
      integer nlig
      real*8 pivot(maxlig,3)
      integer gravity
      real*8 rstk
      pointer (ptr_pivot,pivot)
            
c     Local variables
      integer n,nn,i,jl
      real*8 et,xd,yd,zd,xf,yf,zf
      
      call cartstate_f_disre(cartstatehandle,nlig,ptr_pivot)      
      call ministate_f_disre(ministatehandle,gravity,rstk)
      
      if (gravity.eq.0) return
      
      jl=3*iori*(nlig-fixre)
      
      do 100 n=1+fixre,nlig
      i = jl + 3 * (n-fixre-1)
      xf = 0.0d0
      yf = 0.0d0
      zf = 0.0d0

      if (gravity.eq.3) then
c     gravity=3 => to all other centers
      
      do 90 nn=1,nlig           
      xd=xa(n)+pivot(n,1)-xa(nn)-pivot(nn,1)
      yd=ya(n)+pivot(n,2)-ya(nn)-pivot(nn,2)
      zd=za(n)+pivot(n,3)-za(nn)-pivot(nn,3)
      et=(xd**2+yd**2+zd**2)
      erest=erest+rstk*et*et
      xf=xf-4.0d0*rstk*et*xd
      yf=yf-4.0d0*rstk*et*yd
      zf=zf-4.0d0*rstk*et*zd      
 90   continue      
      
      else 

      if (gravity.eq.1) then
c     gravity=1 => to global origin      
      xd=xa(n)+pivot(n,1)
      yd=ya(n)+pivot(n,2)
      zd=za(n)+pivot(n,3)
      else
c     gravity=2 => to receptor origin
      xd=xa(n)+pivot(n,1)-xa(1)-pivot(1,1)
      yd=ya(n)+pivot(n,2)-ya(1)-pivot(1,2)
      zd=za(n)+pivot(n,3)-za(1)-pivot(1,3)
      endif
      
      et=(xd**2+yd**2+zd**2)

c     Fourth order restraints, maybe a bit too steep...
c      erest=erest+rstk*et*et
c      xf=xf-4.0d0*rstk*et*xd
c      yf=yf-4.0d0*rstk*et*yd
c      zf=zf-4.0d0*rstk*et*zd       

c     Harmonic restraints
      erest=erest+rstk*et
      xf=xf-2.0d0*rstk*xd
      yf=yf-2.0d0*rstk*yd
      zf=zf-2.0d0*rstk*zd       

      endif

      if ((iab.eq.1).and.(n.gt.1.or.fixre.eq.0)) then
      
      delta(i+1) = delta(i+1) + xf
      delta(i+2) = delta(i+2) + yf
      delta(i+3) = delta(i+3) + zf
      end if
      
100   continue      
      
      return
      end
