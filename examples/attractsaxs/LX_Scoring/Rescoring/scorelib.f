      subroutine score(coorrec, rlen, coorlig, llen, rshift, rcut,
     1 power1, power2, epsi, sigi, ivori, hexp, rmax,
     2 atypsrec, atypslig, atyps, energy)
cf2py intent(out) :: energy
cf2py intent(hide) :: rlen
cf2py intent(hide) :: llen
cf2py intent(hide) :: atyps
      real*8 rshift2, rcut, power1, power2, rshift
      integer rlen
      integer llen, atyps
      real*8 epsi(atyps, atyps), sigi(atyps, atyps)
      real*8 rmini(atyps, atyps), hexp(atyps, atyps)
      real*8 rmax(atyps, atyps), emini(atyps, atyps)
      integer ivori(atyps, atyps)
      integer i, j
      real*8 coorrec(rlen,3)
      real*8 coorlig(llen,3)
      integer atypsrec(rlen)
      integer atypslig(llen)
      real*8 energy
      real*8 xr, yr, zr, xl, yl, zl, dist2
      integer ar, al, ivor
      real*8 dist, sig, eps, rmin, emin, rm, rma, he
c
      pree = (power2/power1)**(power1/(power1-power2))
     1       -(power2/power1)**(power2/(power1-power2))
      prer = (power1/power2)**(1.0d0/(power1-power2))
      do 300, i =1, atyps
       do 400, j=1, atyps
        emini(i,j) = pree*epsi(i,j)
        rmini(i,j) = prer*sigi(i,j)
  400  continue
  300 continue
      energy=0.0d0 
      rshift2=rshift**2
      do 100, i=1, rlen
       xr=coorrec(i,1)
       yr=coorrec(i,2)
       zr=coorrec(i,3)
       ar=atypsrec(i)
       do 200, j=1, llen
        xl=coorlig(j,1)
        yl=coorlig(j,2)
        zl=coorlig(j,3)
        al=atypslig(j)
        dist2=(xr-xl)**2+(yr-yl)**2+(zr-zl)**2
        if (dist2 .le. rcut**2) then
         if (dist2 .lt. rshift2) then
          dist2=rshift2
         endif
         dist = sqrt(dist2) 
         sig = sigi(ar,al)
         eps = epsi(ar,al)
         ivor = ivori(ar,al)
         rmin = rmini(ar,al)
         emin = emini(ar,al)
         rm = rmax(ar,al)
         rma = rmin+rm
         he = hexp(ar,al)
c         write(*,*) i, j, dist, ar, al, sig, eps, ivor, rmin, emin, rm, he
         if (dist .le. rmin) then
c          write(*,*) '1'
	  energy = energy + eps*(((sig**power1)/(dist**power1))
     1                    -((sig**power2)/(dist**power2)))
     2                    +(ivor-1)*emin     
     3      +(ivor*he*exp(-(dist-2.0d0*rmin+sig)**2/(rmin-sig)**2))
	 else if (dist .le. rma) then
	  if (ivor .eq. 1) then
c	   write(*,*) '2a'
	   energy = energy+ emin
	  else
c	   write(*,*) '2b'
	   energy= energy +ivor*(eps*(((sig**power1)/(dist**power1))
     1      -((sig**power2)/(dist**power2)))
     2      +(he*exp(-(dist-2.0d0*rmin+sig)**2/(rmin-sig)**2)))
          endif
         else
	  if (ivor .eq. 1) then
c	  write(*,*) '3a'
	  energy=energy+ivor*(eps*(((sig**power1)/((dist-rm)**power1))
     1      -((sig**power2)/((dist-rm)**power2)))
     2      +(he*exp(-((dist-rm)-2.0d0*rmin+sig)**2/(rmin-sig)**2)))
	  else
c	   write(*,*) '3b'
	   energy=energy +ivor*(eps*(((sig**power1)/(dist**power1))
     1      -((sig**power2)/(dist**power2)))
     2      +(he*exp(-(dist-2.0d0*rmin+sig)**2/(rmin-sig)**2)))
          endif
	 endif
c	 write(*,*) i, j, dist, energy
        endif
  200  continue
  100 continue
      return
      end
c      
      subroutine stepscore(rcut, coorrec, rlen, coorlig, llen,
     1 atypsrec, atypslig, param, bins, atyps, score)
cf2py intent(out) :: score
cf2py intent(hide) :: rlen
cf2py intent(hide) :: llen
cf2py intent(hide) :: bins
cf2py intent(hide) :: atyps
      integer bins
      real rcut(bins)
      integer rlen
      integer llen
      integer i, j, atyps, b
      real coorrec(rlen,3)
      real coorlig(llen,3)
      real*8 param(bins, atyps,atyps)
      real*8 score
      integer atypsrec(rlen)
      integer atypslig(llen)
      real xr, yr, zr, xl, yl, zl, dist2, rend2, r2
      integer ar, al, s
c gtypes is the number of parameter: gtypes = atyps*(atpys-1)/2 +atyps
c modatyp is atyps calculated from gtypes
      score=0.0d0
      rend2=rcut(bins)
      do 100, i=1, rlen
       xr=coorrec(i,1)
       yr=coorrec(i,2)
       zr=coorrec(i,3)
       ar=atypsrec(i)
       do 200, j=1, llen
        xl=coorlig(j,1)
        yl=coorlig(j,2)
        zl=coorlig(j,3)
        al=atypslig(j)
        dist2=(xr-xl)**2+(yr-yl)**2+(zr-zl)**2
        if (dist2 .le. rend2) then
	 do 300, b=1, bins
         r2=rcut(b)
         if (dist2 .le. r2) then
          score= score+param(b,ar,al)
          exit
         endif
  300   continue
        endif
  200  continue
  100 continue
      return
      end
c      
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