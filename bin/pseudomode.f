      subroutine psmgroup_by_residue(maxatom,maxmode,
    2  natom,nres,iaci,iei,npsmgroups,psmgroups)
      implicit none 
      
c     Prepares pseudomode groups by assigning each residue to a different group

c     parameters:      input
      integer maxatom,maxmode
      integer natom,nres
      
      integer iaci, iei
      dimension iaci(maxatom),iei(maxatom)
      
c     parameters:      output
      integer npsmgroups
      
      integer psmgroups
      dimension psmgroups(maxmode,maxatom)

c     local variables      
      integer i,j,resindex
      
      do i=1, nres
      do j=1, maxatom
      psmgroups(i,j) = 0
      enddo      
      enddo
      
      do i=1, natom
      resindex = iei(i)
      psmgroups(resindex,i) = 1
      enddo
      
      npsmgroups = nres
      return

      subroutine psmgroup_by_atom(maxatom,maxmode,
    2  natom,nres,iaci,iei,npsmgroups,psmgroups)
      implicit none 
      
c     Prepares pseudomode groups by assigning each atom to a different group

c     parameters:      input
      integer maxatom,maxmode
      integer natom,nres
      
      integer iaci, iei
      dimension iaci(maxatom),iei(maxatom)
      
c     parameters:      output
      integer npsmgroups
      
      integer psmgroups
      dimension psmgroups(maxmode,maxatom)

c     local variables      
      integer i,j
      
      do i=1, natom
      do j=1, maxatom
      psmgroups(i,j) = 0
      enddo      
      enddo
      
      do i=1, natom
      psmgroups(i,i) = 1
      enddo

      npsmgroups = natom
      return

      end
      
      subroutine psmgroup_to_psm(max3atom,maxatom,maxmode,
    2  natom,npsmgroups,psmgroups,pseudomodes)
      implicit none 
      
c     Converts pseudomode groups to Cartesian pseudomodes, assigning 1.0 and 0.0 to the mode vectors
c     The number of pseudomodes will the three times the number of groups (x y and z)

c     parameters:      input
      integer max3atom,maxatom,maxmode            
      integer natom,npsmgroups
      
      integer psmgroups
      dimension psmgroups(maxmode,maxatom)

c     parameters:      output
      double*8 pseudomodes
      dimension pseudomodes(maxmode,max3atom)

c     local variables      
      integer i,j,k
      real*8 v

      do i=1, 3 * npsmgroups
      do j=1, max3atom
      pseudomodes(i,j) = 0
      enddo      
      enddo

      do i=1, npsmgroups
       do j=1, natom
       v = psmgroups(i,j)
       if (v.gt.0) then
        do k=1, 3 !x,y,z
        pseudomodes(3*(i-1)+k,3*(j-1)+k) = 1.0
        enddo      
       endif
       enddo
      enddo
