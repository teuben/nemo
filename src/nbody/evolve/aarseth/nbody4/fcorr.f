      SUBROUTINE FCORR(I,DM,KW)
*
*
*       Total force corrections due to mass loss.
*       -----------------------------------------
*
      INCLUDE 'common4.h'
      SAVE VP,VI2
      REAL*8  VP(3),A(9)
*
*
*       Save the velocity components and square velocity.
      VI2 = 0.0
      DO 1 K = 1,3
          VP(K) = XDOT(K,I)
          VI2 = VI2 + XDOT(K,I)**2
    1 CONTINUE
*
*       Include velocity kick in case of new neutron star or BH.
      IF(KW.NE.KSTAR(I).AND.(KW.EQ.13.OR.KW.EQ.14))THEN
*       Distinguish between single star and binary (called from EXPEL).
          IF (I.LE.N) THEN
              CALL KICK(I,1)
          END IF
      END IF
*
*       Define consistent c.m. variables for KS mass loss (exclude Roche).
      IF (I.GT.N.AND.KSTAR(I).LE.10) THEN
          I2 = 2*(I - N)
          I1 = I2 - 1
          IF (LIST(1,I1).EQ.0) THEN
              CALL KSRES2(I-N,J1,J2,0.0D0)
          END IF
          VF2 = 0.0
          DV2 = 0.0
          BODYI = BODY(I)
	  if (body(i).eq.0.0d0) bodyi = body(i1) + body(i2)
          DO 5 K = 1,3
              X(K,I) = (BODY(I1)*X(K,I1) + BODY(I2)*X(K,I2))/BODYI
              XDOT(K,I) = (BODY(I1)*XDOT(K,I1) + BODY(I2)*XDOT(K,I2))/
     &                                                           BODYI
              X0(K,I) = X(K,I)
              X0DOT(K,I) = XDOT(K,I)
              VF2 = VF2 + XDOT(K,I)**2
              DV2 = DV2 + (XDOT(K,I) - VP(K))**2
    5     CONTINUE
          VFAC = SQRT(VF2/VI2)
      END IF
*
*       Exclude self-interaction for KS component.
      IF (I.LT.IFIRST) THEN
          J = KVEC(I) + N
      ELSE
          J = I
      END IF
*
*       Form potential on host for any common-envelope mass as first case.
      DMS = DM*ZMBAR
      IF (IPHASE.EQ.0.AND.J.GT.N) THEN
          CALL POTI(J,POTJ)
*       Obtain potential for non-standard case (DMS > 0.1 & NS done in MDOT).
      ELSE IF (IPHASE.NE.-3) THEN
          CALL POTI(J,POTJ)
      ELSE IF (DMS.GT.0.1.OR.KW.GE.13) THEN
          POTJ = 0.0D0
      ELSE
*       Decide between GRAPE value with correction and direct summation.
          IF ((STEP(J).GT.0.004.AND.DMS.LT.0.01).OR.T0(J).EQ.TIME.OR.
     &        DMS.LT.0.001) THEN
              CALL PHIDOT(J,DPHI)
              POTJ = PHI(J) + DPHI
          ELSE
              CALL POTI(J,POTJ)
          END IF
      END IF
*
*       Perform explicit restart instead of corrections if DMSUN > 0.1.
      IF (ABS(DMS).GT.0.1)THEN
*        IF(I.LE.N) GO TO 50
         IF(I.GT.N.AND.KSTAR(I).LE.10) GO TO 50
      END IF
*
*       Define neigbour distances and copy membership from NBLIST.
      FACM = MIN(1.0 + DMS,3.0D0)
      RS2 = (FACM*RSCALE)**2/FLOAT(N)**0.66667
      RCR2 = 4.0*RS2
      NNB = ILIST(1)
*
*       Correct neighbour forces & first derivatives.
      DO 40 L = 2,NNB+1
          J = ILIST(L)
          RIJ2 = 0.0D0
*
          DO 10 K = 1,3
              A(K) = X(K,I) - X(K,J)
              RIJ2 = RIJ2 + A(K)**2
   10     CONTINUE
*
*       Skip force corrections for distant particles.
*         IF (RIJ2.GT.RCR2) GO TO 40
*
          RIJDOT = 0.0D0
          RDVDOT = 0.0D0
*
          DO 12 K = 1,3
              A(K+3) = VP(K) - XDOT(K,J)
              RIJDOT = RIJDOT + A(K)*A(K+3)
              RDVDOT = RDVDOT + A(K)*(XDOT(K,I) - VP(K))
   12     CONTINUE
*
          A3 = 1.0/(RIJ2*SQRT(RIJ2))
          A4 = BODY(I)*A3
          A5 = DM*A3
          A6 = 3.0*RIJDOT/RIJ2
          A7 = 3.0*RDVDOT/RIJ2
*
          DO 15 K = 1,3
              A(K+3) = (A(K+3) - A6*A(K))*A5
              IF(KW.EQ.13.OR.KW.EQ.14)THEN
*       Include FDOT corrections due to increased velocity.
                  A(K+3) = A(K+3) + (XDOT(K,I) - VP(K))*A4
                  A(K+3) = A(K+3) - A7*A(K)*A4
              ENDIF
   15     CONTINUE
*
*       Subtract contributions to potential, force & first derivative.
          PHI(J) = PHI(J) + DM/SQRT(RIJ2)
          DO 20 K = 1,3
              F(K,J) = F(K,J) - 0.5*A(K)*A5
              FDOT(K,J) = FDOT(K,J) - ONE6*A(K+3)
   20     CONTINUE
   40 CONTINUE
*
*       Update the potential and kinetic energy loss.
   50 ECDOT = ECDOT + DM*POTJ + 0.5*DM*VI2
*
*     IF (IPHASE.NE.-3.AND.KW.GE.10) THEN
*     WRITE (6,52) I, IPHASE, ILIST(1), DMS, DM*POTJ,DM*PHI(I)
*  52 FORMAT (' FCORR    I IPH NB DMS DM*POTJ DM*PHI ',
*    &                   2I6,I4,F8.3,1P,2E10.2)
*     END IF
*
*       Modify energy loss further for c.m. body (exclude Roche cases).
      IF (I.GT.N.AND.KSTAR(I).LE.10) THEN
          ECDOT = ECDOT - 0.5*BODY(I)*VI2*(VFAC**2 - 1.0)
*       Improve solution for c.m. system (binary or hierarchy).
          CALL FPOLYI(I)
      END IF
*
*       See whether tidal terms should be included (standard or scaled).
      IF (KZ(14).GT.0) THEN
          IF (KZ(14).LE.2) THEN
              ECDOT = ECDOT - 0.5*DM*(TIDAL(1)*X(1,I)**2 +
     &                                TIDAL(3)*X(3,I)**2)
          ELSE
              BODY(I) = BODY(I) + DM
              CALL XTRNLV(I,I)
              ECDOT = ECDOT + HT
              BODY(I) = BODY(I) - DM
              CALL XTRNLV(I,I)
              ECDOT = ECDOT - HT
              ECDOT = ECDOT - 0.5*TIDAL(4)*DM*(X(1,I)*XDOT(2,I) -
     &                                         X(2,I)*XDOT(1,I))
          END IF
      END IF
*
*       Accumulate energy loss for conservation check.
      E(12) = ECDOT
*
      RETURN
*
      END
