      program modesca
c
c this program calculates the Hinsen modes (Proteins, 1998) using only CA atoms of a protein
c usage: $path/modesca structureca.pdb (should contain only CA atoms)
c author: Martin Zacharias, Jacosb University Bremen
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      character*80 name,b
      dimension xc(4000),x(4000),d2f(4000,4000),dx(3),
     1          a(9000000),r(16000000)
      dimension ka(2000),ie(2000),ipoint(2000)
      character*4 ty(2000),rg(2000),at
      call getarg(1,name)
      open(42,file=name)
      i=0
      k=0
      write(*,*)'start'
   40 read(42,20,end=30) b
      if(b(:4).eq.'ATOM') then
      read(b,26) at,ka(i+1),ty(i+1),rg(i+1),ie(i+1),x(3*i+1),
     1       x(3*i+2),x(3*i+3)
      if(ty(i+1).eq."CA  ") then
      ipoint(k+1)=i+1
      k=k+1
      endif
      i=i+1
      endif
      goto 40
   30 close(42)
      ires=k
      ipoint(k)=i+1
      natom=i
      n3=3*natom
      fac=1.0d0
      write(*,*)'ires,natom',ires,natom,n3
c
      do 60 i=1,n3
      do 70 j=1,n3
      d2f(i,j)=0.0d0
   70 continue
   60 continue
c generate second derivative matrix
      do 100 i=1,natom
      ii=3*(i-1)
      do 110 j=1,natom
      if(j.ne.i) then
      jj=3*(j-1)
      r2=xnull
      do 120 k=1,3
      dx(k)=x(ii+k)-x(jj+k)
      r2=r2+dx(k)**2
  120 continue
      fb=fac*exp(-r2/16.0d0)
      b0=2.0d0*fb/r2
      do 140 k=1,3
      do 150 l=1,3
      d2f(ii+k,jj+l)=d2f(ii+k,jj+l)-b0*dx(k)*dx(l)
  150 continue
  140 continue
      do 240 k=1,3
      do 250 l=k,3
      d2f(ii+k,ii+l)=d2f(ii+k,ii+l)+b0*dx(k)*dx(l)
c     d2f(jj+k,jj+l)=d2f(jj+k,jj+l)+b0*dx(k)*dx(l)
  250 continue
  240 continue
      endif
  110 continue
  100 continue
      k=0
      do 480 i=1,n3
c     write(*,'(500f10.5)') (d2f(j,i),j=1,i)
      do 490 j=1,i
      k=k+1
      a(k)=d2f(j,i)
  490 continue
  480 continue
      do 500 i=1,n3*n3
      r(i)=0.0
  500 continue
      itest=n3/6
      irest=n3-6*itest
      write(*,*)'itest,irest',itest,irest
      mv=0
      call eigen(a,r,n3,mv)
      open(90,file='eigen.out')
      do 550 i=1,n3
      ii=i-1
      write(90,*) i,a(i*(i+1)/2)
      write(*,*) i,a(i*(i+1)/2)
      do k=1,itest
      kk=6*(k-1)
      write(90,'(6f15.10)') (r(n3*ii+kk+j),j=1,6)
      enddo
      if(irest.gt.0)write(90,'(6f15.10)') (r(n3*ii+6*itest+j),j=1,irest)
  550 continue 
      close(90)
c generate diagonal elements
   20 format(a80)
   26 format(a4,i7,2x,a4,a4,2x,i3,4x,3f8.3)
      end
c
      SUBROUTINE EIGEN (A,R,N,MV)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
C
C                                                                      C
C     SUBROUTINE EIGEN (A,R,N,MV)                                      C
C                                                                      C
COMMENT   EIGEN COMPUTES EIGENVALUES AND EIGENVECTORS OF THE REAL      C
C     SYMMETRIC N*N MATRIX A, USING THE DIAGONALIZATION METHOD         C
C     DESCRIBED IN "MATHEMATICAL METHODS FOR DIGITAL COMPUTERS", EDS.  C
C     A.RALSTON AND H.S.WILF, WILEY, NEW YORK, 1962, CHAPTER 7.        C
C     IT HAS BEEN COPIED FROM THE IBM SCIENTIFIC SUBROUTINE PACKAGE.   C
C                                                                      C
C     A(1..N*(N+1)/2) = MATRIX TO BE DIAGONALIZED, STORED IN SYMMETRIC C
C                       STORAGE MODE, VIZ. THE I,J-TH ELEMENT (I.GE.J) C
C                       IS STORED AT THE LOCATION K=I*(I-1)/2+J IN A;  C
C                       THE EIGENVALUES ARE DELIVERED IN DESCENDING    C
C                       ORDER ON THE DIAGONAL, VIZ. AT THE LOCATIONS   C
C                       K=I*(I+1)/2                                    C
C     R(1..N,1..N) = DELIVERED WITH THE CORRESPONDING EIGENVECTORS     C
C                    STORED COLUMNWISE                                 C
C     N = ORDER OF MATRICES A AND R                                    C
C     MV = 0 : EIGENVALUES AND EIGENVECTORS ARE COMPUTED               C
C        = 1 : ONLY EIGENVALUES ARE COMPUTED                           C
C                                                                      C
      DIMENSION A(1),R(1)
C
C*****GENERATE IDENTITY MATRIX
    5 RANGE=1.D-12
      IF (MV-1) 10,25,10
   10 IQ=-N
      DO 20 J=1,N
      IQ=IQ+N
      DO 20 I=1,N
      IJ=IQ+I
      R(IJ)=0.E0
      IF (I-J) 20,15,20
   15 R(IJ)=1.E0
   20 CONTINUE
C
C*****COMPUTE INITIAL AND FINAL NORMS (ANORM AND ANRMX)
   25 ANORM=0.E0
      DO 35 I=1,N
      DO 35 J=I,N
      IF (I-J) 30,35,30
   30 IA=I+(J*J-J)/2
      ANORM=ANORM+A(IA)*A(IA)
   35 CONTINUE
      IF (ANORM) 165,165,40
   40 ANORM=1.414E0* SQRT(ANORM)
      ANRMX=ANORM*RANGE/FLOAT(N)
C
C*****INITIALIZE INDICATORS AND COMPUTE THRESHOLD, THR
      IND=0
      THR=ANORM
   45 THR=THR/FLOAT(N)
   50 L=1
   55 M=L+1
C
C*****COMPUT SIN AND COS
   60 MQ=(M*M-M)/2
      LQ=(L*L-L)/2
      LM=L+MQ
   62 IF ( ABS(A(LM))-THR) 130,65,65
   65 IND=1
      LL=L+LQ
      MM=M+MQ
      X=0.5E0*(A(LL)-A(MM))
   68 Y=-A(LM)/ SQRT(A(LM)*A(LM)+X*X)
      IF (X) 70,75,75
   70 Y=-Y
   75 SINX=Y/ SQRT(2.E0*(1.E0+( SQRT(1.E0-Y*Y))))
      SINX2=SINX*SINX
   78 COSX= SQRT(1.E0-SINX2)
      COSX2=COSX*COSX
      SINCS=SINX*COSX
C
C*****ROTATE L AND M COLUMNS
      ILQ=N*(L-1)
      IMQ=N*(M-1)
      DO 125 I=1,N
      IQ=(I*I-I)/2
      IF (I-L) 80,115,80
   80 IF (I-M) 85,115,90
   85 IM=I+MQ
      GOTO 95
   90 IM=M+IQ
   95 IF (I-L) 100,105,105
  100 IL=I+LQ
      GOTO 110
  105 IL=L+IQ
  110 X=A(IL)*COSX-A(IM)*SINX
      A(IM)=A(IL)*SINX+A(IM)*COSX
      A(IL)=X
  115 IF (MV-1) 120,125,120
  120 ILR=ILQ+I
      IMR=IMQ+I
      X=R(ILR)*COSX-R(IMR)*SINX
      R(IMR)=R(ILR)*SINX+R(IMR)*COSX
      R(ILR)=X
  125 CONTINUE
      X=2.E0*A(LM)*SINCS
      Y=A(LL)*COSX2+A(MM)*SINX2-X
      X=A(LL)*SINX2+A(MM)*COSX2+X
      A(LM)=(A(LL)-A(MM))*SINCS+A(LM)*(COSX2-SINX2)
      A(LL)=Y
      A(MM)=X
C
C*****TESTS FOR COMPLETION
C
C*****TEST FOR M = LAST COLUMN
  130 IF (M-N) 135,140,135
  135 M=M+1
      GOTO 60
C
C*****TEST FOR L = SECOND FROM LAST COLUMN
  140 IF (L-(N-1)) 145,150,145
  145 L=L+1
      GOTO 55
  150 IF (IND-1) 160,155,160
  155 IND=0
      GOTO 50
C
C*****COMPARE THRESHOLD WITH FINAL NORM
  160 IF (THR-ANRMX) 165,165,45
C
C*****SORT EIGENVALUES AND EIGENVECTORS
  165 IQ=-N
      DO 185 I=1,N
      IQ=IQ+N
      LL=I+(I*I-I)/2
      JQ=N*(I-2)
      DO 185 J=I,N
      JQ=JQ+N
      MM=J+(J*J-J)/2
      IF (A(LL)-A(MM)) 170,185,185
  170 X=A(LL)
      A(LL)=A(MM)
      A(MM)=X
      IF (MV-1) 175,185,175
  175 DO 180 K=1,N
      ILR=IQ+K
      IMR=JQ+K
      X=R(ILR)
      R(ILR)=R(IMR)
  180 R(IMR)=X
  185 CONTINUE
C
      RETURN
      END
