      SUBROUTINE SUBSYS(NSYS,CM)
*
*
*       Initialization of subsystem.
*       ----------------------------
*
      INCLUDE 'common4.h'
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      REAL*8  CM(10)
*
*
*       Increase subsystem counter and set current index.
      NSUB = NSUB + 1
      ISUB = NSUB
*
*       Set zero name & mass in last component to distinguish triple case.
      IF (NSYS.EQ.3) THEN
          NAMES(4,ISUB) = 0
          BODYS(4,ISUB) = 0.0D0
      END IF
*
*       Save global masses & names of new subsystem components.
      DO 10 L = 1,NSYS
          J = JLIST(L)
          BODYS(L,ISUB) = BODY(J)
          NAMES(L,ISUB) = NAME(J)
   10 CONTINUE
*
*       Form ghosts and initialize integration variables (NSYS = 3 or 4).
      DO 20 L = 1,NSYS
          J = JLIST(L)
          BODY(J) = 0.0D0
          T0(J) = 1.0D+06
          DO 15 K = 1,3
              X0DOT(K,J) = 0.0D0
              XDOT(K,J) = 0.0D0
              F(K,J) = 0.0D0
              FDOT(K,J) = 0.0D0
              D2(K,J) = 0.0D0
              D3(K,J) = 0.0D0
   15     CONTINUE
*       Set large X0 & X to avoid perturber selection (no escape).
          X0(1,J) = 1.0D+06
          X(1,J) = 1.0D+06
   20 CONTINUE
*
*       Place c.m. of subsystem in first location (ICOMP may switch!).
      ICOMP = JLIST(1)
      T0(ICOMP) = TIME
      BODY(ICOMP) = CM(7)
*       Define zero name for identification purpose (only for chain c.m.).
      IF (ISYS(ISUB).EQ.3) NAME(ICOMP) = 0
      DO 30 K = 1,3
          X(K,ICOMP) = CM(K)
          X0(K,ICOMP) = CM(K)
          XDOT(K,ICOMP) = CM(K+3)
          X0DOT(K,ICOMP) = CM(K+3)
   30 CONTINUE
*
*       Predict coordinates & velocities of all other particles (order FDOT).
      DO 40 J = IFIRST,NTOT
          CALL XVPRED(J,0)
   40 CONTINUE
*
*       Construct force polynomial for c.m. motion (NLIST mod not needed).
      CALL FPOLY1(ICOMP,ICOMP,0)
*
*       Initialize decision-making variables for multiple regularization.
      T0S(ISUB) = TIME
      TS(ISUB) = TIME
      STEPS(ISUB) = STEP(ICOMP)
*
*       Obtain maximum unperturbed separation based on dominant neighbour.
      CALL EXTEND(ISUB)
*
      RETURN
*
      END
