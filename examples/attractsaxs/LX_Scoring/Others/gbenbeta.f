c
c calculate born radii
      subroutine gbrad(corm, rsolv, fs, delta, kam, rborn)
CF2PY INTENT(OUT) :: rborn
CF2PY INTENT(HIDE) :: kam
      real*8 corm(kam,3)
      real*8 rsolv(kam)
      real*8 fs(kam)
      INTEGER kam
      real*8 r(kam)
      real*8 rborn(kam)
      real*8 sam(kam)
      real*8 delta
c
      do 50 i = 1, kam
      r(i) = rsolv(i)-delta
   50 continue
c
      do 60 i = 1,kam
       sam(i) = 1.0d0 / r(i)
   60 continue
       do 70 i = 1,kam-1
        xi = corm(i,1)
        yi = corm(i,2)
        zi = corm(i,3)
        ri = r(i)
        si = ri * fs(i)
        si2 = si * si
        do 80 k = i+1, kam
         rk = r(k)
         sk = rk * fs(k)
         sk2 = sk * sk
         dik = sqrt((corm(k,1)-xi)**2+(corm(k,2)-yi)**2+
     1              (corm(k,3)-zi)**2)
          rdik=1.0d0/dik
          if (ri .lt. dik+sk) then
           wik = 1.0d0 / max(ri,dik-sk)
           uik = 1.0d0 / (dik+sk)
           wik2 = wik * wik
           uik2 = uik * uik
           duw2=uik2-wik2
           term=wik-uik+0.25d0*dik*duw2+(0.5d0*rdik)*log(uik/wik)
     1                -(0.25d0*sk2*rdik)*duw2
           sam(i)=sam(i)-0.5d0*term
          end if
          if (rk .lt. dik+si) then
           wik = 1.0d0 / max(rk,dik-si)
           uik = 1.0d0 / (dik+si)
           wik2 = wik * wik
           uik2 = uik * uik
           duw2=uik2-wik2
           term=wik-uik+0.25d0*dik*duw2+(0.5d0*rdik)*log(uik/wik)
     1                -(0.25d0*si2*rdik)*duw2 
           sam(k) = sam(k) - 0.5d0*term
          end if
   80   continue
   70  continue
        do 90 i = 1, kam
            rborn(i) = 1.0d0 / sam(i)
            rborn(i) = max(r(i),rborn(i))
   90   continue
	end
c
c
c calculate born radii of complexes
      subroutine gbradcomplex(corm, recsolv, ligsolv,
     1 recfs, ligfs, recborn, ligborn, delta, rkam, lkam, kam, rborn)
CF2PY INTENT(OUT) :: rborn
CF2PY INTENT(HIDE) :: rkam
CF2PY INTENT(HIDE) :: lkam
CF2PY INTENT(HIDE) :: kam
      INTEGER rkam
      INTEGER lkam
      integer kam,i,k
      real*8 corm(kam,3)
      real*8 recsolv(rkam)
      real*8 ligsolv(lkam)
      real*8 recfs(rkam)
      real*8 ligfs(lkam)
      real*8 recborn(rkam)
      real*8 ligborn(lkam)
      real*8 r(kam)
      real*8 rborn(kam)
      real*8 sam(kam)
      real*8 xi,yi,zi,ri,si,si2,rk,sk,sk2,dik,rdik,wik,uik,wik2,uik2,
     1 duw2, term, delta
c
      do 50 i = 1, rkam
      r(i) = recsolv(i)-delta
      sam(i) = 1.0d0/recborn(i)
   50 continue
c
      do 60 i = 1,lkam
      r(i+rkam) = ligsolv(i)-0.12d0
      sam(i+rkam) = 1.0d0/ligborn(i)
   60 continue
       do 70 i = 1,rkam
        xi = corm(i,1)
        yi = corm(i,2)
        zi = corm(i,3)
        ri = r(i)
        si = ri * recfs(i)
        si2 = si * si
        do 80 k = 1, lkam
         rk = r(k+rkam)
         sk = rk * ligfs(k)
         sk2 = sk * sk
         dik = sqrt((corm(k+rkam,1)-xi)**2+(corm(k+rkam,2)-yi)**2+
     1              (corm(k+rkam,3)-zi)**2)
          rdik=1.0d0/dik
          if (ri .lt. dik+sk) then
           wik = 1.0d0 / max(ri,dik-sk)
           uik = 1.0d0 / (dik+sk)
           wik2 = wik * wik
           uik2 = uik * uik
           duw2=uik2-wik2
           term=wik-uik+0.25d0*dik*duw2+(0.5d0*rdik)*log(uik/wik)
     1                -(0.25d0*sk2*rdik)*duw2
           sam(i)=sam(i)-0.5d0*term
          end if
          if (rk .lt. dik+si) then
           wik = 1.0d0 / max(rk,dik-si)
           uik = 1.0d0 / (dik+si)
           wik2 = wik * wik
           uik2 = uik * uik
           duw2=uik2-wik2
           term=wik-uik+0.25d0*dik*duw2+(0.5d0*rdik)*log(uik/wik)
     1                -(0.25d0*si2*rdik)*duw2 
           sam(k+rkam) = sam(k+rkam) - 0.5d0*term
          end if
   80   continue
   70  continue
        do 90 i = 1, kam
            rborn(i) = 1.0d0 / sam(i)
            rborn(i) = max(r(i),rborn(i))
   90   continue
	end
c
c     calculate born energy wth given born radi
      subroutine gbenergy(corm, rborn, dmon, ensol,kam)
CF2PY INTENT(HIDE) :: kam
CF2PY INTENT(OUT) :: ensol
      real*8 corm(kam,3)
      INTEGER kam
      real*8 rborn(kam)
      real*8 ensol
      real*8 dmon(kam)
      ensol = 0.0d0
      convk=332.053986d0
      dwater = 78.3d0
      f = -convk * (1.0d0 - 1.0d0/dwater)
c      ecoulomb=0.0d0
      do 110 i = 1,kam
         xi = corm(i,1)
         yi = corm(i,2)
         zi = corm(i,3)
         fi = f * dmon(i)
c         fcc = convk*dmon(i)
	 fact=0.5d0*fi*dmon(i)
	 ensol= ensol + fact/rborn(i)
         do 120 k = i+1,kam 
            xr = xi - corm(k,1)
            yr = yi - corm(k,2)
            zr = zi - corm(k,3)
            r2 = xr*xr + yr*yr + zr*zr
            fik = fi * dmon(k)
c            fccik= fcc* dmon(k)
            rb2 = rborn(i) * rborn(k)
            expterm = exp(-0.25d0*r2/rb2)
            rfgb = 1.0d0/sqrt(r2 + rb2*expterm)
c            ecoulomb=ecoulomb+fccik/sqrt(r2)
            e = fik*rfgb
            ensol = ensol + e
  120    continue
  110 continue
c      write(*,*)'ensol',ensol
      return
      end
c
c
c     calculate born energy wth given born radi
      subroutine gbenfixed(cormr,corml,rborn,dmonr,dmonl,
     1 ensol,kamr,kaml)
CF2PY INTENT(HIDE) :: kamr
CF2PY INTENT(HIDE) :: kaml
CF2PY INTENT(OUT) :: ensol
      real*8 cormr(kamr,3)
      real*8 corml(kaml,3)
      INTEGER kaml
      integer kamr
      real*8 rborn
      real*8 ensol
      real*8 dmonr(kamr)
      real*8 dmonl(kaml)
      real*8 dwater, f, convk, rb2, xi, yi, zi, fi
      real*8 expterm, rfgb, e, xr, yr, zr, r2, fik
      ensol = 0.0d0
      convk=332.053986d0
      dwater = 78.3d0
      f = -convk * (1.0d0 - 1.0d0/dwater)
c      ecoulomb=0.0d0
      rb2 = rborn * rborn
      do 110 i = 1,kamr
         xi = cormr(i,1)
         yi = cormr(i,2)
         zi = cormr(i,3)
         fi = f * dmonr(i)
         do 120 k = 1,kaml 
            xr = xi - corml(k,1)
            yr = yi - corml(k,2)
            zr = zi - corml(k,3)
            r2 = xr*xr + yr*yr + zr*zr
            fik = fi * dmonl(k)
c            fccik= fcc* dmon(k)
            expterm = exp(-0.25d0*r2/rb2)
            rfgb = 1.0d0/sqrt(r2 + rb2*expterm)
c            ecoulomb=ecoulomb+fccik/sqrt(r2)
            e = fik*rfgb
            ensol = ensol + e
  120    continue
  110 continue
c      write(*,*)'ensol',ensol
      return
      end
c
c
c     calculate born energy wth given born radi
      subroutine elec(cormr,corml, dmonr,dmonl,shift,ensol,kamr,kaml)
CF2PY INTENT(HIDE) :: kamr
CF2PY INTENT(HIDE) :: kaml
CF2PY INTENT(OUT) :: ensol
      real*8 cormr(kamr,3)
      real*8 corml(kaml,3)
      INTEGER kaml
      integer kamr
      real*8 ensol, shift
      real*8 dmonr(kamr)
      real*8 dmonl(kaml)
      real*8 dwater, f, convk, xi, yi, zi, fi
      real*8 rfgb, e, xr, yr, zr, r2, fik
      ensol = 0.0d0
      convk=332.053986d0/10.0d0
      dwater = 78.3d0
      f = convk
      shiftoff=shift*shift
c      *(1.0d0 - 1.0d0/dwater)
      do 110 i = 1,kamr
         xi = cormr(i,1)
         yi = cormr(i,2)
         zi = cormr(i,3)
         fi = f * dmonr(i)
         do 120 k = 1,kaml 
            xr = xi - corml(k,1)
            yr = yi - corml(k,2)
            zr = zi - corml(k,3)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .lt. shiftoff) then
	     r2=shiftoff
	    endif
            fik = fi * dmonl(k)
            rfgb = 1.0d0/sqrt(r2)
            e = fik*rfgb
            ensol = ensol + e
  120    continue
  110 continue
c      write(*,*)'ensol',ensol
      return
      end
c
c recalculate born radii
      subroutine gbradcmpx(corm, lenr, rsolv, fs, rbornin, kam, rborn)
CF2PY INTENT(OUT) :: rborn
CF2PY INTENT(HIDE) :: kam
      real*8 corm(kam,3)
      real*8 rsolv(kam)
      real*8 fs(kam)
      INTEGER kam
      real*8 r(kam)
      real*8 rborn(kam)
      real*8 rbornin(kam)
      real*8 sam(kam)
c
      do 60 i = 1,kam
       sam(i) = 1.0d0 / rbornin(i)
       r(i) = rsolv(i)-0.12d0
   60 continue
       do 70 i = 1,lenr
        xi = corm(i,1)
        yi = corm(i,2)
        zi = corm(i,3)
        ri = r(i)
        si = ri * fs(i)
        si2 = si * si
        do 80 k = lenr, kam
         rk = r(k)
         sk = rk * fs(k)
         sk2 = sk * sk
         dik = sqrt((corm(k,1)-xi)**2+(corm(k,2)-yi)**2+
     1              (corm(k,3)-zi)**2)
	 if (dik .lt. 30) then
          rdik=1.0d0/dik
          if (ri .lt. dik+sk) then
           wik = 1.0d0 / max(ri,dik-sk)
           uik = 1.0d0 / (dik+sk)
           wik2 = wik * wik
           uik2 = uik * uik
           duw2=uik2-wik2
           term=wik-uik+0.25d0*dik*duw2+(0.5d0*rdik)*log(uik/wik)
     1                -(0.25d0*sk2*rdik)*duw2
           sam(i)=sam(i)-0.5d0*term
          end if
          if (rk .lt. dik+si) then
           wik = 1.0d0 / max(rk,dik-si)
           uik = 1.0d0 / (dik+si)
           wik2 = wik * wik
           uik2 = uik * uik
           duw2=uik2-wik2
           term=wik-uik+0.25d0*dik*duw2+(0.5d0*rdik)*log(uik/wik)
     1                -(0.25d0*si2*rdik)*duw2 
           sam(k) = sam(k) - 0.5d0*term
          end if
         end if
   80   continue
   70  continue
        do 90 i = 1, kam
            rborn(i) = 1.0d0 / sam(i)
            rborn(i) = max(r(i),rborn(i))
   90   continue
	end
c
c recalculate born radii
      subroutine gbencmpx(corm,lenr,rsolv,fs, dmon,rbornin,kam,ensol)
CF2PY INTENT(OUT) :: ensol
CF2PY INTENT(HIDE) :: kam
      real corm(kam,3)
      real rsolv(kam)
      real dmon(kam)
      real fs(kam)
      INTEGER kam
      real r(kam)
      real rborn(kam)
      real rbornin(kam)
      real sam(kam)
c
      do 60 i = 1,kam
       sam(i) = 1.0d0 / rbornin(i)
       r(i) = rsolv(i)-0.12d0
   60 continue
       do 70 i = 1,lenr
        xi = corm(i,1)
        yi = corm(i,2)
        zi = corm(i,3)
        ri = r(i)
        si = ri * fs(i)
        si2 = si * si
        do 80 k = lenr, kam
         rk = r(k)
         sk = rk * fs(k)
         sk2 = sk * sk
         dik = sqrt((corm(k,1)-xi)**2+(corm(k,2)-yi)**2+
     1              (corm(k,3)-zi)**2)
	 if (dik .lt. 35) then
          rdik=1.0d0/dik
          if (ri .lt. dik+sk) then
           wik = 1.0d0 / max(ri,dik-sk)
           uik = 1.0d0 / (dik+sk)
           wik2 = wik * wik
           uik2 = uik * uik
           duw2=uik2-wik2
           term=wik-uik+0.25d0*dik*duw2+(0.5d0*rdik)*log(uik/wik)
     1                -(0.25d0*sk2*rdik)*duw2
           sam(i)=sam(i)-0.5d0*term
          end if
          if (rk .lt. dik+si) then
           wik = 1.0d0 / max(rk,dik-si)
           uik = 1.0d0 / (dik+si)
           wik2 = wik * wik
           uik2 = uik * uik
           duw2=uik2-wik2
           term=wik-uik+0.25d0*dik*duw2+(0.5d0*rdik)*log(uik/wik)
     1                -(0.25d0*si2*rdik)*duw2 
           sam(k) = sam(k) - 0.5d0*term
          end if
         end if
   80   continue
   70  continue
        do 90 i = 1, kam
            rborn(i) = 1.0d0 / sam(i)
            rborn(i) = max(r(i),rborn(i))
   90   continue
c      
      ensol = 0.0d0
      convk=332.053986d0
      dwater = 78.3d0
      f = -convk * (1.0d0 - 1.0d0/dwater)
c      ecoulomb=0.0d0
      do 110 i = 1,kam
         xi = corm(i,1)
         yi = corm(i,2)
         zi = corm(i,3)
         fi = f * dmon(i)
c         fcc = convk*dmon(i)
	 fact=0.5d0*fi*dmon(i)
	 ensol= ensol + fact/rborn(i)
         do 120 k = i+1,kam 
            xr = xi - corm(k,1)
            yr = yi - corm(k,2)
            zr = zi - corm(k,3)
            r2 = xr*xr + yr*yr + zr*zr
            fik = fi * dmon(k)
c            fccik= fcc* dmon(k)
            rb2 = rborn(i) * rborn(k)
            expterm = exp(-0.25d0*r2/rb2)
            rfgb = 1.0d0/sqrt(r2 + rb2*expterm)
c            ecoulomb=ecoulomb+fccik/sqrt(r2)
            e = fik*rfgb
            ensol = ensol + e
  120    continue
  110 continue
c      write(*,*)'ensol',ensol
      return
      end