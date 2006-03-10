      SUBROUTINE CORE
*
*
*       Density centre & core radius.
*       -----------------------------
*
      INCLUDE 'common4.h'
      REAL*8 RLIST(10*LMAX),RHOS
      REAL*8 R2,RHO
      COMMON/WORK1/ R2(NMAX),RHO(NMAX)
      INTEGER KLIST(10*LMAX)
*
*
*       Introduce N-dependent limit for small systems.
      K = 1 + FLOAT(N)**0.3333
      NBCRIT = MIN(K,6)
*       Estimate length scale for sufficient membership.
      RCORE2 = 0.6d0*RSCALE**2
      IF (TIME.GT.0.0) RCORE2 = MAX(4.d0*RC**2,RCORE2)
*
*       Select central particles inside RCORE2.
    5 NC = 0
      DO 10 I = IFIRST,NTOT
          RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2
     &                                 + (X(3,I) - RDENS(3))**2
*       Save square radius of central members in temporary array for NBLIST.
          IF (RI2.LT.RCORE2) THEN
              NC = NC + 1
              R2(NC) = RI2
              JLIST(NC) = I
          END IF
   10 CONTINUE
*
*       Check membership in range > MIN(N/5,20) and < N/2 (N > 50 on call).
      IF (NC.LT.N/5.OR.NC.LT.20) THEN
          RCORE2 = 1.5d0*RCORE2
          GO TO 5
      END IF
      IF (NC.GT.N/2) THEN
          RCORE2 = 0.66d0*RCORE2
          GO TO 5
      END IF
*
*       Sort NC square radii sequentially in JLIST.
      CALL SORT1(NC,R2,JLIST)
*
*       Define basic neighbour radius squared.
      RS2 = 0.04d0*RSCALE**2/FLOAT(N)**0.66667
*
*       Obtain individual densities.
      RHOS = 0.d0
      DO 50 L = 1,NC
          I = JLIST(L)
          RI2 = R2(L)
          IF (RS2.GT.0.16*RI2) THEN
              RS2 = 0.16d0*RI2
          END IF
*       Determine neighbour list of #I (note: RS2 is modified for next I).
          RS2I = RS2
          CALL NCLIST(I,RS2,RI2,RS2I)
          NNB = ILIST(1)
*       Estimate D6 to limit the candidates above NNB = 10.
          D6 = 1.6d0*(6.d0/FLOAT(NNB))**0.66*RS2I
          XI = X(1,I)
          YI = X(2,I)
          ZI = X(3,I)
*
   20     N6 = 0
          DO 25 LJ = 2,NNB+1
              J = ILIST(LJ)
              RIJ2 = (XI - X(1,J))**2 + (YI - X(2,J))**2 +
     &                                  (ZI - X(3,J))**2
              IF (RIJ2.LT.D6) THEN
                  N6 = N6 + 1
                  KLIST(N6) = J
                  RLIST(N6) = RIJ2
              END IF
   25     CONTINUE
*
*       Make another iteration in rare cases.
          IF (N6.LT.NBCRIT) THEN
              D6 = 1.5d0*D6
              IF (N6.LT.NNB) GO TO 20
          END IF
*
*       Sort list of square distances.
          DO 40 II = 1,N6
              DO 35 JJ = II+1,N6
                  IF (RLIST(JJ).LT.RLIST(II)) THEN
                      RDUM = RLIST(II)
                      IDUM = KLIST(II)
                      RLIST(II) = RLIST(JJ)
                      KLIST(II) = KLIST(JJ)
                      RLIST(JJ) = RDUM
                      KLIST(JJ) = IDUM
                  END IF
   35         CONTINUE
   40     CONTINUE
*
          I6 = MIN(N6,6)
*       Calculate total mass of five nearest neighbours.
          XMASS = 0.d0
          DO 45 LJ = 1,I6-1
              XMASS = XMASS + BODY(KLIST(LJ))
   45     CONTINUE
*
          RHO(I) = XMASS/(RLIST(I6)*SQRT(RLIST(I6)))
*       Assign zero weight if not enough neighbours.
          IF (N6.LT.6) RHO(I) = 0.0D0
          RHOS = MAX(RHOS,RHO(I))
   50 CONTINUE
*
*       Determine density centre.
      DO 60 K = 1,3
          RDENS(K) = 0.d0
   60 CONTINUE
*
      RHO1 = 0.d0
      DO 70 L = 1,NC
          I = JLIST(L)
          DO 65 K = 1,3
              RDENS(K) = RDENS(K) + RHO(I)*X(K,I)
   65     CONTINUE
          RHO1 = RHO1 + RHO(I)
   70 CONTINUE
*
*       Set current density centre based on improved determination.
      DO 75 K = 1,3
          RDENS(K) = RDENS(K)/RHO1
          IF (KZ(5).EQ.2.AND.TTOT.LT.TCRIT) THEN
              RDENS(K) = 0.0
          END IF
   75 CONTINUE
*
*       Obtain density radius & averaged density.
      RC = 0.d0
      RHO2 = 0.d0
      DO 80 L = 1,NC
          I = JLIST(L)
          XID = X(1,I) - RDENS(1)
          YID = X(2,I) - RDENS(2)
          ZID = X(3,I) - RDENS(3)
          RID2 = XID**2 + YID**2 + ZID**2
          RC = RC + RHO(I)**2*RID2
          RHO2 = RHO2 + RHO(I)**2
   80 CONTINUE
*
*       Form core radius and average & maximum density (scaled in OUTPUT).
      RC = SQRT(RC/RHO2)
      RHOD = RHO2/RHO1
      RHOD = (3.d0/12.566d0)*RHOD
      RHOM = (3.d0/12.566d0)*RHOS
*       NB! Scaling factor 5 of Casertano & Hut eq. (V.3) replaced by XMASS.
      IF (RC.EQ.0.0) THEN
          RC = RSCALE
          NC = N/2
      END IF
*
*       Sum particles & mass inside the core radius and set rms velocity.
      NC1 = 0
      VC = 0.d0
      ZMC = 0.d0
      RC2 = RC**2
      DO 85 L = 1,NC
          I = JLIST(L)
          XID = X(1,I) - RDENS(1)
          YID = X(2,I) - RDENS(2)
          ZID = X(3,I) - RDENS(3)
          RID2 = XID**2 + YID**2 + ZID**2
          IF (RID2.LT.RC2) THEN
              NC1 = NC1 + 1
              VC = VC + XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
              ZMC = ZMC + BODY(I)
          END IF
   85 CONTINUE
*
*       Set core membership & rms velocity.
      NC = MAX(NC1,2)
      VC = SQRT(VC/FLOAT(NC))
*
      RETURN
*
      END
