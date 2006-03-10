      SUBROUTINE MYDUMP(I,J)
*
*
*       COMMON save or read.
*       --------------------
*
      INCLUDE 'params.h'
      PARAMETER  (NA=56*NMAX,NB=111*KMAX+7,
     &            NC=NMAX+LMAX+(2*LMAX+3)*KMAX+177,
     &            NF=44*MMAX,NG=17*NMAX+5*KMAX+410,
     &            NH=2*NMAX+84,NM=32*NTMAX,NQ=20*MCL+20,
     &            NR=31*MMAX)
      REAL*4  A,B,C,E,F,G,H,Q,R,S
      INTEGER  IC,ID,IR
*
      COMMON/NBODY/  A(NA)
      COMMON/PAIRS/  B(NB)
      COMMON/NAMES/  IC(NC)
      COMMON/COUNTS/ ID(76)
      COMMON/PARAMS/ E(326)
      COMMON/BINARY/ F(NF)
      COMMON/STARS/  G(NG)
      COMMON/BLOCKS/ H(NH)
      COMMON/MODES/  C(NM)
      COMMON/RAND2/  IR(99)
      COMMON/CLOUDS/ Q(NQ)
      COMMON/RCHE/   R(NR)
      COMMON/GALAXY/ S(40)
*
*
*       Open unit #J by reading dummy and rewinding.
      REWIND J
      READ (J,ERR=10,END=10)  DUMMY
   10 REWIND J
*
*       Read or save all COMMON variables (valid for tape or disc).
      IF (I.EQ.0) THEN
          READ (J)   A, B, IC, ID, E, F, G, H, IR, C, Q, R, S
      ELSE
          WRITE (J)  A, B, IC, ID, E, F, G, H, IR, C, Q, R, S
          CLOSE (UNIT=J)
      END IF
*
      RETURN
*
      END
