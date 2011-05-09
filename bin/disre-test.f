      subroutine disre
c  adds a restraining contribution due to attraction between 
c center of receptor and clostest CA atom of the ligand protein
c
      include "CADD.t"
      erest=0.0d0
c     jj=3*(jrst-1)
      jj = 3*(1455-1137-1)
      jjj=3*(903-1)
      xd=y(iconf,jj+1)- yc(iconf,jjj+1)
      yd=y(iconf,jj+2)- yc(iconf,jjj+2)
      zd=y(iconf,jj+3)- yc(iconf,jjj+3)
      et=(xd**2+yd**2+zd**2)

      limsq = 4.0d0
      violation = et - limsq
      write(*,*),'VIOLATION',violation
      if (et.lt.limsq) return
      cforce = 1.0d0
      erest = 0.5 * cforce * violation
      factor = sqrt(violation/et)
      if(iab.eq.1) then
      write(*,*),'FORCES',-cforce*factor*xd,
     1 -cforce*factor*yd,-cforce*factor*zd
      f(jj+1)=f(jj+1)-cforce*factor*xd
      f(jj+2)=f(jj+2)-cforce*factor*yd
      f(jj+3)=f(jj+3)-cforce*factor*zd
      endif

c      erest=rstk*et*et
c      write(*,*),'EREST',erest, et, xd,y(iconf,jj+1),ce(1)
c      if(iab.eq.1) then
c      f(jj+1)=f(jj+1)-4.0d0*rstk*et*xd
c      f(jj+2)=f(jj+2)-4.0d0*rstk*et*yd
c      f(jj+3)=f(jj+3)-4.0d0*rstk*et*zd
c      endif

      
      return
      end
