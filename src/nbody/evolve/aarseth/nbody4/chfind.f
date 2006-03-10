      SUBROUTINE CHFIND
*
*
*       Identification of global chain index.
*       -------------------------------------
*
      INCLUDE 'common4.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX),LISTCM(LMAX)
*
*
*       Re-determine global chain index & perturber list after switching.
      DO 10 I = IFIRST,N
          IF (NAME(I).EQ.0.AND.BODY(I).GT.0.0) THEN
              ICH = I
*       Form new neighbour list for chain c.m.
              CALL CMLIST
*       Update the chain perturber list consistently with new sequence.
              CALL CHLIST(ICH)
*       Resolve global coordinates & velocities of subsystem.
              CALL XCPRED(0)
              GO TO 20
          END IF
   10 CONTINUE
*
   20 RETURN
*
      END
