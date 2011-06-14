      SUBROUTINE GGUBS (DSEED,NR,R)                                     GGUS0390
C                                  SPECIFICATIONS FOR ARGUMENTS         GGUS0400
      INTEGER            NR                                             GGUS0410
      REAL*8             R(NR)                                          GGUS0420
      DOUBLE PRECISION   DSEED                                          GGUS0430
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   GGUS0440
      INTEGER            I                                              GGUS0450
      DOUBLE PRECISION   D2P31M,D2P31                                   GGUS0460
C                                  D2P31M=(2**31) - 1                   GGUS0470
C                                  D2P31 =(2**31)(OR AN ADJUSTED VALUE) GGUS0480
      DATA               D2P31M/2147483647.D0/                          GGUS0490
      DATA               D2P31/2147483711.D0/                           GGUS0500
C                                  FIRST EXECUTABLE STATEMENT           GGUS0510
      DO 5 I=1,NR                                                       GGUS0520
         DSEED = DMOD(16807.D0*DSEED,D2P31M)                            GGUS0530
      R(I) = DSEED / D2P31                                              GGUS0540
c     write(*,*)'pseudo', r(i)
    5 continue
      RETURN                                                            GGUS0550
      END                                                               GGUS0560
