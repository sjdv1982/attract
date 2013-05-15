c     NOTE: Make sure that 3*npsmgroups is not larger than maxmode! If necessary, change MAXMODE in max.h and recompile
c     NOTE: These arrays can get pretty big; If short on memory, reduce MAXATOM in max.h and recompile

      subroutine prepare_springs(max3atom,maxatom,maxmode,
    2  natom,nres,xori0,iaci,iei,springfile,eigfile,
    3  npsmgroups,psmgroups,pseudomodes,spring,nhm0,eig,
    4  nrw,rw,rr,rrv)
      implicit none 

c  prepare spring arrays for a single RNA/protein molecule (called "ligand")
c   to be used with the kick routine
c  input variables:
c    natom: number of atoms in the ligand
c    nres: number of residues in the ligand
c    xori0: the original PDB coordinates, before any deformation or pivot centering applied
c    iaci: atom type of each atom in the ligand
c    iei: residue index of each atom in the ligand (starting at 1)
c    springfile: filename of the elastic network (with <natom> x <natom> values)
c    eigfile: filename of the *base* harmonic modes (<nhm> modes of <natom x 3> mode coordinates)
c  output variables:
c    npsmgroups: number of pseudomode groups that will control the deformation
c     npsmgroups can be <natom>, <nres> or any other number
c    psmgroups: a <npsmgroups> x <natom> array of booleans (mask), 
c     indicating whether or not an atom belongs to that pseudomode group
c    pseudomodes: (3 * <npsmgroups> modes of <natom x 3> mode coordinates (double) )
c     typically, these will be masks of the underlying atomic x y z components, 
c      closely related to the psmgroups above
c     For example, there can be one pseudomode group per residue
c      In that case, there will be three pseudomodes per residue, indicating residue displacement in x,y and z
c      Then, the x-pseudomode for residue 10 will have (1.0, 0.0, 0.0) for all atoms in residue 10, 
c      and (0.0, 0.0, 0.0) for all other atoms
c      In comparison, the 10th pseudomode group will have 1 for all atoms in residue 10, and 0 for all other atoms
c     Alternatively, there could be one pseudomode group (i.e 3 pseudomodes: x/y/z) per atom 
c      In that case, the z-pseudomode for atom 13 will have (0.0, 0.0, 1.0) for atom 13 
c       and (0.0, 0.0, 0.0) otherwise
c    spring: an array of <npsmgroups> x <npsmgroups> spring constants
c     This will typically be the contents of <springfile>, but with the size reduced from <natom> to <npseudomodes>
c      for example, by averaging over all atoms in a residue, or by selecting a single atom per residue
c    nhm0: number of base harmonic modes, see below
c    eig: an array of <nhm0> x <npsmgroups x 3> base harmonic modes
c     This will typically be the contents of <eigfile>, but with the size reduced from <natom x 3> to <npsmgroups x 3>
c      for example, by averaging over all atoms in a residue, or by selecting a single atom per residue
c     These numbers, multiplied with a local filter, are used by the kick routine to generate values for the pseudomode displacements
c    nrw: number of local filters, see below. In the current implementation, <nrw> is <natom>, 
c      and the filters are Gaussian
c    rw: a <nrw> x <npsmgroups x 3> array describing the local Gaussian filters
c    rrv: a <npsmgroups> x <npsmgroups x 3> array describing the distance matrix 
c      of the pseudomode groups
c     rrv(i,j) = atom(j) - atom(i)
c    rr: a <npsmgroups> x <npsmgroups> array describing the scalar sizes of rrv

c     parameters:      input
      integer max3atom,maxatom,maxmode
      integer natom,nres
      
      real*8 xori0
      dimension xori0(max3atom)
      integer iaci, iei
      dimension iaci(maxatom),iei(maxatom)

      character*100 springfile
      character*100 eigfile
      
c     parameters:      output
      integer npsmgroups,nhm0,nrw
      
      integer psmgroups
      dimension psmgroups(maxmode,maxatom)
      real*8 pseudomodes
      dimension pseudomodes(maxmode,max3atom)
      real*8 spring
      dimension spring(maxatom,maxatom)
      real*8 eig
      dimension eig(maxmode,max3atom)
      real*8 rw
      dimension rw(maxatom,max3atom)
      real*8 rrv
      dimension rrv(maxmode,maxmode,3)
      real*8 rr
      dimension rr(maxmode,maxmode)

c     local variables
      real*8 spring0
      dimension spring0(maxatom,maxatom)      
      real*8 eig0
      dimension eig0(maxmode,max3atom)      
      double val_eig0
      dimension val_eig0(maxmode)
      integer neig0
      
      integer i,j,k, is_protein, atomtype
      integer is_selected
      integer representative
      dimension representative(maxmode)
      real*8 groupavg
      dimension groupavg(3*maxmode)
      integer curratom,curratom2
      integer nselected
      dimension nselected(maxmode)
      integer selected
      dimension selected(maxatom)
      real*8 avgx,avgy,avgz
      real*8 cx,cy,cz
      real*8 dx,dy,dz,rd
      real*8 cl
      real*8 springscale
           
c     Read contents of spring file       
      open(51,file=springfile)
      do 100 i=1,natom
      read(51,*)(spring0(i,j),j=1,natom)
  100 continue

c     Read base harmonic modes
      call read_hm(eigfile, "spring ligand", 1, 
     2 natom, neig0, val_eig0, eig0, 0)
      
c     Create pseudomode atom groups, one group per residue    
      call psmgroup_by_residue(maxatom,maxmode,
    2  natom,nres,iaci,iei,npsmgroups,psmgroups)
c     (Alternatively, you could invoke psmgroup_by_atom, or write your own subroutine)
    
c  selection loop: select the relevant atoms belonging to each group
      
      do i=1,npsmgroups

c     Are we dealing with a protein residue group?
      is_protein=1
      do j=1,natom
       if (psmgroups(i,j).gt.0) then
        atomtype = iei(j)
        if ((atomtype.gt.32).AND.(atomtype.lt.65)) then
         is_protein=0
         break
        endif
       endif  
      enddo

c     Obtain selected atoms
c       For protein (contains no atomtypes 32 < x < 65): the CA (atom type 1)
c       For non-protein (contains atomtypes 32 < x < 65): any atom
      nselected(i)=0
      do j=1,natom
       if (psmgroups(i,j).gt.0) then
        atomtype = iei(j)
        is_selected = 0
        if (is_protein.gt.0) then
          if (atomtype.eq.1) then
            is_selected = 1
          endif
        else !non-protein
          is_selected = 1
        endif
        if (is_selected.gt.0) then
         nselected(i) = nselected(i) + 1
         selected(nselected(i)) = j 
        endif
       endif  
      enddo
      
c     TODO: assert that at least one atom has been selected (Fortran stop if not?)

c     Reducing: we want to keep only the spring/base mode data that corresponds to the selected atoms

c     Calculate the average group coordinates    
      avgx = 0.0d0
      avgy = 0.0d0
      avgz = 0.0d0
      do j=1,nselected(i)
       curratom = selected(k)
       avgx = avgx + xori0(3*(curratom-1)+1)
       avgy = avgy + xori0(3*(curratom-1)+2)
       avgz = avgz + xori0(3*(curratom-1)+3)
      enddo
      avgx = avgx / nselected(i)
      avgy = avgy / nselected(i)
      avgz = avgz / nselected(i)          
      groupavg(3*(i-1)+1) = avgx
      groupavg(3*(i-1)+2) = avgy
      groupavg(3*(i-1)+3) = avgz

c     Take the first selected atom as representative for springs reducing
      representative(i) = selected(1)

c     Reduce the base harmonic modes, selecting the average mode value for the selected atoms
      do j=1,neig0
       avgx = 0.0d0
       avgy = 0.0d0
       avgz = 0.0d0
       do k=1,nselected(i)
        curratom = selected(k)
        avgx = avgx + eig0(j,3*(curratom-1)+1)
        avgy = avgy + eig0(j,3*(curratom-1)+2)
        avgz = avgz + eig0(j,3*(curratom-1)+3)
       enddo
       avgx = avgx / nselected(i)
       avgy = avgy / nselected(i)
       avgz = avgz / nselected(i)    
       eig(j,3*(i-1)+1) = avgx
       eig(j,3*(i-1)+2) = avgy
       eig(j,3*(i-1)+3) = avgz
      enddo
      
c  END of selection loop
      enddo      
      
c     Reduce the springs based on the representative atoms
      do i=1,npsmgroups
      curratom = representative(i)
      do j=1,npsmgroups
      curratom2 = representative(j)
      spring(i,j) = spring0(curratom,curratom2)
c     Since we are representing many distances, we may want to scale with group size
      springscale = nselected(i) * nselected(j)
      spring(i,j) = spring(i,j) * springscale
c
      enddo
      enddo

c     Compute matrix of inter-group distance vectors (j - i) (rrv) 
c      and distances (rr)
c      using the group averages
      do i=1,npsmgroups
      cx = groupavg(3*(i-1)+1)
      cy = groupavg(3*(i-1)+2)
      cz = groupavg(3*(i-1)+3)
      do j=1,npsmgroups
      dx = groupavg(3*(j-1)+1) - cx
      dy = groupavg(3*(j-1)+2) - cy
      dz = groupavg(3*(j-1)+3) - cz
      rrv(i,j,1) = dx
      rrv(i,j,2) = dy
      rrv(i,j,3) = dz
      rr(i,j) = sqrt(dx**2+dy**2+dz**2)
      enddo
      enddo

c     Compute local Gaussian filters for the kick routine
c      using the group averages
c     There will be <natom> local filters
      nrw = natom
c     Alternatively, you could define <npsmgroups> filters
c
c     correlation length for kick range:
      cl=0.001
c     note: if cl->0 then "kick" becomes a casual normal mode deformation
      do i=1,natom
      cx = xori0(3*(i-1)+1)
      cy = xori0(3*(i-1)+2)
      cz = xori0(3*(i-1)+3)
      do j=1,npsmgroups
      dx = groupavg(3*(j-1)+1) - cx
      dy = groupavg(3*(j-1)+2) - cy
      dz = groupavg(3*(j-1)+3) - cz
      rd = dx**2+dy**2+dz**2
      rw(i,j)=exp(-cl*rd)
      rw(j,i)=rw(i,j)      
      enddo
      enddo
      
c     Build pseudomode vectors from the groups
      call psmgroup_to_psm(max3atom,maxatom,maxmode,
    2  natom,npsmgroups,psmgroups,pseudomodes)

      nhm0 = neig0 ! Reducing does not reduce the *number* of base modes
      
      return
      end
            
