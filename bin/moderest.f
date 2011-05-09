      subroutine moderest(maxdof,maxmode,
     1                  dlig,nhm,val,delta,
     2                  enlig)
      implicit none
      
c     Parameters      
      real*8 enlig
      integer maxlig,maxdof,maxmode
      real*8 dlig,val,delta
      integer nhm
      dimension delta(maxdof),val(maxmode),dlig(maxdof)


c     Local variables
      real*8 deloldl,en
      integer i

c
c  add the implicit part
c
      do 50 i=1,nhm
      en=val(i)*dlig(i)**4
      enlig=enlig+en
      deloldl=4.0d0*val(i)*dlig(i)**3
      delta(6+i)=delta(6+i)+deloldl
   50 continue
      return 
      end
