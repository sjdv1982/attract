       subroutine globalenergy(
     1	cartstatehandle, ministatehandle,
     2  ens,phi,ssi,rot,xa,ya,za,morph,dlig,
     3  locrests, has_locrests, seed,
     4  iab,iori,itra,ieig,iindex,fixre,
     5  energies, delta, deltamorph)

       implicit none
       
c      Parameters
       include "max.fin"
       integer cartstatehandle,ministatehandle
       integer iab,iori,itra,ieig,iindex,fixre,seed
       real*8 energies(6)
c      energies(1): vdw
c      energies(2): elec
c      energies(3): restraints (distance and symmetry)
c      energies(4): flexible displacement (modes and morphing)
c      energies(5): origin displacement (gravity and location restraints)
c      energies(6): atom density grid
       real*8 delta(maxdof)

       real*8 locrests
       dimension locrests(3,maxlig)
       integer has_locrests      
       dimension has_locrests(maxlig)

       integer ens
       dimension ens(maxlig)
       real*8 phi, ssi, rot, morph, dlig, xa, ya, za
       dimension phi(maxlig), ssi(maxlig), rot(maxlig)
       dimension dlig(maxmode+maxindexmode, maxlig)
       dimension xa(maxlig), ya(maxlig), za(maxlig)
       dimension morph(maxlig)
       real*8 deltamorph(maxlig)       

c      Handle variables: cartstate ligand
       real*8 xl(3*maxatom), fl(3*maxatom)
       pointer(ptr_xl,xl)
       pointer(ptr_fl,fl)

c      Handle variables: full coordinates and modes
       integer nlig,nall,nall3
       real*8 xb,xold,x,xori,xori0,f,eig, val,pivot
       integer index_eig
       real*8 index_val
       real*8 ff, frot
       real*8 morph_fconstant
       integer ieins,nhm,nihm,natom,iaci_old
       dimension xb(3*totmaxatom),x(3*totmaxatom),f(3*totmaxatom)
       dimension xold(3*totmaxatom)
       dimension xori(3*totmaxatom),xori0(3*totmaxatom)
       dimension eig(3*maxatom,maxmode,maxlig)
       dimension index_eig(maxlenindexmode,maxindexmode,maxlig)
       dimension index_val(maxlenindexmode,maxindexmode,maxlig)
       dimension val(maxlig,maxmode)
       dimension ieins(maxlig),nhm(maxlig), nihm(maxlig)
       dimension pivot(maxlig,3)
       dimension ff(maxlig)
       dimension frot(9,maxlig)
       dimension natom(maxlig)
       dimension iaci_old(maxlig)
       pointer(ptr_xb,xb)
       pointer(ptr_x,x)
       pointer(ptr_xori,xori)
       pointer(ptr_xori0,xori0)       
       pointer(ptr_f, f)
       pointer(ptr_eig,eig)
       pointer(ptr_val,val)
       pointer(ptr_index_eig,index_eig)
       pointer(ptr_index_val,index_val)
       pointer(ptr_ieins,ieins)
       pointer(ptr_nhm,nhm)
       pointer(ptr_nihm,nihm)
       pointer(ptr_pivot,pivot)
       pointer(ptr_ff, ff)
       pointer(ptr_frot, frot)
       pointer(ptr_natom,natom)
       pointer(ptr_iaci_old,iaci_old)
       real*8 ensd
       dimension ensd(3*maxatom)
       pointer(ptr_ensd, ensd)
       
c      Local variables   
       integer i,ii,n,jl,jb, ju, jn
       real*8 cdelta(6+maxmode+maxindexmode), rotmat(9)
       real*8 enligl,xnull
       real*8 pm2(3,3,3)
       real*8 rotmatinv(9)
       real*8 flcopy
       dimension flcopy(3*maxatom)
       real*8 dmmy1 
       integer dmmy2
       real*8 dmmy3
       real*8 dmmy4
       real*8 zero
       real*8 pivotnull
       dimension pivotnull(3)
       real*8 energy0
       integer has_globalenergy
       integer nmodes, nimodes
       integer, parameter :: ERROR_UNIT = 0
       dmmy3 = -1.0d0
       zero=0.0d0
       pivotnull(1) = zero
       pivotnull(2) = zero
       pivotnull(3) = zero
c      Is global energy used at all?
       call ministate_has_globalenergy(ministatehandle, 
     1  has_globalenergy)

       if (has_globalenergy.gt.0) then
c reset forces
       call reset_forces(cartstatehandle)
  
c get parameters      
       call cartstate_f_globalenergy(cartstatehandle,
     1  nlig, nall, nall3, morph_fconstant, ptr_nhm,ptr_nihm,ptr_ieins,
     2  ptr_eig, ptr_val, ptr_index_eig, ptr_index_val,
     3  ptr_xb, ptr_x, ptr_xori, ptr_xori0,
     4  ptr_f,ptr_pivot, ptr_natom, ptr_iaci_old)
       call cartstate_get_forcerot(cartstatehandle,ptr_ff,ptr_frot)
      
       jl=3*iori*(nlig-fixre)
       jb=3*iori*(nlig-fixre)+3*itra*(nlig-fixre)
       nmodes = 0
       nimodes = 0
       do i=1,nlig
       nmodes = nmodes + nhm(i)
       nimodes = nimodes + nihm(i)
       enddo
       ju = jb + ieig*nmodes
       jn = ju + iindex*nimodes
c apply ensemble/normal mode/index mode deformations
       if (has_globalenergy.eq.1) then
       do 5 i=1, nlig
       call cartstate_get_ensd(cartstatehandle,i-1,ens(i),ptr_ensd,
     1  -1.0d0,dmmy1,dmmy2)
       call deform(ens(i),ensd,dmmy3,dmmy4,dlig(:,i),
     1  nhm,nihm,i-1,ieins,eig,index_eig, index_val,xb,x,xori,xori0,0)
5      continue
c TODO Add index modes
c apply symmetry restraints

       
       xold(1:nall3) = x(1:nall3)
c       call memcpy(xold,x,nall3*8)
       
c      rotate to the global frame
       do 6 i=1, nlig
       call euler2rotmat(phi(i),ssi(i),rot(i),rotmat)
       call rotate(maxlig,3*totmaxatom,rotmat,xa(i),ya(i),za(i),pivot,
     1  i-1,ieins,x)
6      continue

c       write(*,'(a2,3f8.3)'), 'G',x(3*902+1),x(3*902+2),x(3*902+3)
c       write(*,'(a2,3f8.3)'), 'G',x(3*1454+1),x(3*1454+2),x(3*1454+3)
c       write(*,*),(x(3*902+1)-x(3*1454+1))* (x(3*902+1)-x(3*1454+1))+ 
c     1  (x(3*902+2)-x(3*1454+2))* (x(3*902+2)-x(3*1454+2)) +
c     2  (x(3*902+3)-x(3*1454+3))* (x(3*902+3)-x(3*1454+3))
c       stop

       call restrain(ministatehandle,cartstatehandle,seed,
     1  iab,energy0)
       energies(3) = energies(3) + energy0
c      WRITE(ERROR_UNIT,*) "Restrain energy", energy0
       call sym(cartstatehandle, iab, energy0)
       energies(3) = energies(3) + energy0
       call atomdensitygridenergy(energy0,nall,x,xori,iaci_old,
     1  nlig,natom,f,iab)
       energies(6) = energies(6) + energy0

       x(1:nall3) = xold(1:nall3)
c       call memcpy(x,xold,nall3*8)
c endif has_globalenergy.eq.1
       endif 
c calculate DOF deltas
       do 10 i=1, nlig

       xnull = 0.0d0   
       do 2 ii=1,6+maxmode+maxindexmode
       cdelta(ii) = xnull
2      continue  
       
       call cartstate_select_ligand2(cartstatehandle,i-1,ptr_xl,ptr_fl)  

       call forcerotscale(3*maxatom,
     1  ff(i),frot(:,i),natom(i),f(i))
       if (has_globalenergy.eq.1)then
       if ((iori.eq.1).AND.(i.gt.fixre)) then
c       write(*,'(a4,i3,f8.3,f8.3,f8.3,f8.3,f8.3,f8.3)'),
c     1  'DOFS',i,phi(i),ssi(i),rot(i),xa(i),ya(i),za(i)
       call euler2torquemat(phi(i),ssi(i),rot(i),pm2)
       call rota(xl,fl,cdelta,pm2,natom(i))
       ii = 3 * (i-fixre-1)
       delta(ii+1) = delta(ii+1) + cdelta(1)
       delta(ii+2) = delta(ii+2) + cdelta(2)
       delta(ii+3) = delta(ii+3) + cdelta(3)
       endif             

       if ((itra.eq.1).AND.(i.gt.fixre)) then
       call trans(fl,cdelta,natom(i))
       ii = jl + 3 * (i-fixre-1)
       delta(ii+1) = delta(ii+1) + cdelta(4)
       delta(ii+2) = delta(ii+2) + cdelta(5)
       delta(ii+3) = delta(ii+3) + cdelta(6)
       endif
       endif
       if (ieig.eq.1) then
       if (has_globalenergy.eq.1)then
c      rotate forces into ligand frame
       call euler2rotmat(phi(i),ssi(i),rot(i),rotmat)
       call matcopy(rotmat, rotmatinv)
       call matinv(rotmatinv)
       flcopy(:) = fl(:)
       call rotate1(3*maxatom,
     1  rotmatinv,zero,zero,zero,
     2   pivotnull, natom(i),flcopy)
       call ligmin(flcopy,natom(i),i,eig,nhm(i),cdelta)
       endif
       call moderest(maxdof,maxmode,dlig(:,i),nhm(i),val(:,i),
     1  cdelta, energies(3))      

      call moderest(maxdof,maxmode,dlig(:,i),nhm(i),val(i,:),
     1  cdelta, energy0)
       energies(4) = energies(4) + energy0
c      WRITE(ERROR_UNIT,*) "Moderest energy", energy0
       	
       ii = jb
       do 23 n=1,i-1
       ii = ii + nhm(n)
23     continue
       do 24 n=1,nhm(i)
       delta(ii+n) = delta(ii+n) + cdelta(6+n)
24     continue
       endif
       

     
       if ((iindex.eq.1).and.(has_globalenergy.eq.1)) then
c      rotate forces into ligand frame
       call euler2rotmat(phi(i),ssi(i),rot(i),rotmat)
       call matcopy(rotmat, rotmatinv)
       call matinv(rotmatinv)
       flcopy(:) = fl(:)
       call rotate1(3*maxatom,
     1  rotmatinv,zero,zero,zero,
     2   pivotnull, natom(i),flcopy)
       call ligmin_index(flcopy,natom(i),i,index_eig,index_val,
     1                    nhm(i),nihm(i),cdelta)

     
       ii = ju
       do 13 n=1,i-1
       ii = ii + nihm(n)
13     continue
       do 14 n=1,nihm(i)
       delta(ii+n) = delta(ii+n) + cdelta(6+nhm(i)+n)
14     continue
       endif
       
10     continue      
c      end if (has_globalenergy.gt.0)
       endif 

       call disre(maxlig,cartstatehandle,ministatehandle,
     1  iab,iori,itra,fixre,xa,ya,za,
     2  locrests, has_locrests,
     3  delta,energy0)
       energies(5) = energies(5) + energy0

       call ene_morph(morph_fconstant, morph, deltamorph, 
     1  nlig, energy0)
       energies(4) = energies(4) + energy0     

       end
