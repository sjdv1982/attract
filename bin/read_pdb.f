c     Variables for which I don't know what they do:
c      ialst, ialpst, ceest, dice, ce, disa, xrst, yrst, zrst, jrst
c     Changed variables: nlig => natom, n3l => n3atom, ntotlig => nlig, 
c       nresl => nres, nalllig => nall, nalllig3 => nall3

      subroutine read_one_pdb(maxlig, totmaxres, totmaxatom,
     1 pdbfile,kai,tyi,rgi,iei,x,iaci,xlai,
     2 icop,we,chai,ncop,nmaxco,natco,
     3 nlig, nres, natom, n3atom, nall, nall3, ieins, ieins3)
      
c      
c  eingabe1:
c       nlig
c       file1.pdb
c       file2.pdb
c       ....
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)      
      character*100 pdbfile
      integer totmaxatom, totmaxres, maxlig
      dimension kai(totmaxatom),
     1 iei(totmaxatom), x(totmaxatom), iaci(totmaxatom),
     2 xlai(totmaxatom), icop(totmaxatom), we(totmaxatom),
     3 chai(totmaxatom)
      dimension ncop(0:10,0:20,totmaxres),nmaxco(totmaxres),
     1 natco(totmaxres)
      dimension ieins(0:maxlig-1),ieins3(0:maxlig-1),
     1 natom(0:maxlig-1),n3atom(0:maxlig-1)
      dimension nres(0:maxlig-1)
      character*4 tyi(totmaxatom),rgi(totmaxatom)
      
      i=0
      irs=0

      open(42,file=pdbfile)
      read(42,*) nlig
      do 60 ijk=0,nlig-1            
      call read_pdb(42, maxlig, totmaxatom, totmaxres, ijk,nlig,
     1 kai,tyi,rgi,iei,x,iaci,xlai,
     2 icop,we,chai,ncop,nmaxco,natco,
     3 nres, natom, n3atom, i, irs, ieins, ieins3)
   60 continue
      close(42)
      nall=i
      nall3=3*i    
c     write(*,*)'all ligand atoms',nall,(natom(j),j=0,nlig-1)
      end

      subroutine read_single_pdb(maxlig, totmaxres, totmaxatom,
     1 pdbfile,kai,tyi,rgi,iei,x,iaci,xlai,
     2 icop,we,chai,ncop,nmaxco,natco,
     3 nlig, nres, natom, n3atom, nall, nall3, ieins, ieins3,i,irs)
      
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)      
      character*100 pdbfile
      integer totmaxatom, totmaxres, maxlig
      dimension kai(totmaxatom),
     1 iei(totmaxatom), x(totmaxatom), iaci(totmaxatom),
     2 xlai(totmaxatom), icop(totmaxatom), we(totmaxatom),
     3 chai(totmaxatom)
      dimension ncop(0:10,0:20,totmaxres),nmaxco(totmaxres),
     1 natco(totmaxres)
      dimension ieins(0:maxlig-1),ieins3(0:maxlig-1),
     1 natom(0:maxlig-1),n3atom(0:maxlig-1)
      dimension nres(0:maxlig-1)
      character*4 tyi(totmaxatom),rgi(totmaxatom)

      open(42,file=pdbfile)
      call read_pdb(42, maxlig, totmaxatom, totmaxres, 0,nlig,
     1 kai,tyi,rgi,iei,x,iaci,xlai,
     2 icop,we,chai,ncop,nmaxco,natco,
     3 nres, natom, n3atom, i, irs, ieins, ieins3)
      close(42)
      nall=i
      nall3=3*i    
c     write(*,*)'all ligand atoms',nall,(natom(j),j=0,nlig-1)
      end

      subroutine read_two_pdbs(maxlig, totmaxres, totmaxatom,
     1 pdbfile1,pdbfile2,kai,tyi,rgi,iei,x,iaci,xlai,
     2 icop,we,chai,ncop,nmaxco,natco,
     3 nres, natom, n3atom, nall, nall3, ieins, ieins3)
      
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)      
      character*100 pdbfile1,pdbfile2
      integer totmaxatom, totmaxres, maxlig
      dimension kai(totmaxatom),
     1 iei(totmaxatom), x(totmaxatom), iaci(totmaxatom),
     2 xlai(totmaxatom), icop(totmaxatom), we(totmaxatom),
     3 chai(totmaxatom)
      dimension ncop(0:10,0:20,totmaxres),nmaxco(totmaxres),
     1 natco(totmaxres)
      dimension ieins(0:maxlig-1),ieins3(0:maxlig-1),
     1 natom(0:maxlig-1),n3atom(0:maxlig-1)
      dimension nres(0:maxlig-1)
      character*4 tyi(totmaxatom),rgi(totmaxatom)
      
      
      
      i=0
      irs=0

      open(42,file=pdbfile1)
      call read_pdb(42, maxlig, totmaxatom, totmaxres, 0,2,
     1 kai,tyi,rgi,iei,x,iaci,xlai,
     2 icop,we,chai,ncop,nmaxco,natco,
     3 nres, natom, n3atom, i, irs, ieins, ieins3)
      close(42)
      open(42,file=pdbfile2)
      call read_pdb(42, maxlig, totmaxatom, totmaxres, 1,2,
     1 kai,tyi,rgi,iei,x,iaci,xlai,
     2 icop,we,chai,ncop,nmaxco,natco,
     3 nres, natom, n3atom, i, irs, ieins, ieins3)
      close(42)
      nall=i
      nall3=3*i    
c      write(*,*)'all ligand atoms',nall,(natom(j),j=0,2-1)
      end
            
      subroutine read_pdb(filehandle,maxlig,totmaxatom,totmaxres,      
     1 ijk,nlig,kai,tyi,rgi, iei,x,iaci,xlai,
     2 icop,we,chai,ncop,nmaxco,natco,
     3 nres, natom, n3atom, i, irs, ieins, ieins3)

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)      
      character*100 b
      integer maxlig, totmaxatom, totmaxres, filehandle
      character*4 at, tyi(totmaxatom),rgi(totmaxatom)
      dimension kai(totmaxatom),
     1 iei(totmaxatom), x(totmaxatom), iaci(totmaxatom),
     2 xlai(totmaxatom), icop(totmaxatom), we(totmaxatom),
     3 chai(totmaxatom)
      dimension ncop(0:10,0:20,totmaxres),nmaxco(totmaxres),
     1 natco(totmaxres)
      dimension ieins(0:maxlig-1),ieins3(0:maxlig-1),
     1 natom(0:maxlig-1),n3atom(0:maxlig-1)
      dimension nres(0:maxlig-1)
      

      iold=0
      irso=0
      izz=1

c open and read ligand protein file
  200 read(filehandle,20,end=226) b
      if(b(:3).eq.'TER') goto 226
      if(b(:4).eq.'ATOM') then
       ii=3*i
       if (i.gt.totmaxatom) then
         write (*,*), "TOTMAXATOM exceeded:", totmaxatom
         stop
       endif
       read(b,26) at,kai(i+1),tyi(i+1),rgi(i+1),iei(i+1),x(ii+1),
     1  x(ii+2),x(ii+3),iaci(i+1),xlai(i+1),icop(i+1),we(i+1)
c       chai(i+1)=felec*xlai(i+1) #done later
       chai(i+1) = xlai(i+1)
       kai(i+1)=i+1
c
c this for renumbering of residues
c and identification of conformational copies of side chains or loop segments
c
        if(iei(i+1).ne.irso) then
         irso=iei(i+1)
c
c irs: residue counter
c
         irs=irs+1
         iold=0
        endif
c renumbering of residues
       iei(i+1)=irs
c check if copy
        if(icop(i+1).ne.0) then
c izz: copy counter
        izz=izz+1
        else
        izz=0
        endif
       if(icop(i+1).ne.iold) then
        izz=1
        iold=icop(i+1)
       endif
c array for the atoms (izz) of each copy (iold) of a residue (irs)
       ncop(iold,izz,irs)=i+1
c max number of copies of residue irs 
       nmaxco(irs)=icop(i+1)
c max number of copies on residue irs 
       natco(irs)=izz
       i=i+1
      endif
      goto 200

  226 nres(ijk)=iei(i)
      ieins(ijk)=i
      ieins3(ijk)=3*i
      natom(ijk)=i
      if (ijk.ge.0) natom(ijk) = natom(ijk)-ieins(ijk-1)
      if (natom(ijk).gt.maxatom) then
        write (*,*), "MAXATOM exceeded:", maxatom
        stop      
      endif
      n3atom(ijk)=3*natom(ijk)
c      write(*,*)'nres,ieins',ijk,nres(ijk),ieins(ijk),natom(ijk),
c     1           ieins3(ijk),nlig

   20 format(a100)
   26 format(a4,i7,2x,a4,a4,2x,i3,4x,3f8.3,i5,f8.3,i2,f5.2)   

      end
      
      subroutine apply_permi(totmaxatom,nall,chai,permi)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)      
      integer totmaxatom
      dimension chai(totmaxatom)

c     Electric permittivity constant (epsilon); 1 = vacuum

      felec=332.053986d0

c     The felec constant is the electrostatic energy, in kcal/mol,      
c      between two electrons, at 1 A distance, in vacuum
c     The formula is:
c      e**2 * NA * KC * Ang * kcal 
c     where:
c      e = charge of the electron, in Coulomb
c      NA = Avogadro's number
c      KC = Coulomb force constant
c      Ang = the size of an Angstrom (10**-10 meter)
c      kcal = the amount of Joules per kcal, 4184
c

      felec=sqrt(felec/permi)
c      write(*,*), "FELEC", felec/permi
c      stop
c      This is because the charges are multiplied with each other
c      Multiplying every charge with sqrt(felec) allows removing felec
c       from the calculation
       
      do i=1,nall
      chai(i) = felec*chai(i)
      end do
      

      end
