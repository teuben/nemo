      SUBROUTINE KSRECT(IPAIR)
*
*
*       Rectification of KS orbit.
*       --------------------------
*
      INCLUDE 'common4.h'
*
*
*       Skip procedure for circularized orbits.
      I = N + IPAIR
      SEMI = -0.5*BODY(I)/H(IPAIR)
      ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(SEMI*BODY(I))
      IF(ECC2.EQ.0.D0) GO TO 50
      IF(ECC2.LT.0.0)THEN
         WRITE(6,*)' KSRECT ERROR ECC < 0 ',ipair,ecc2
         WRITE(6,*)body(i),r(ipair),h(ipair),tdot2(ipair)
         CALL gpfree
         STOP
      ENDIF
      ECC = SQRT(ECC2)
      IF (ECC.LE.0.002) GO TO 50
*     IF (ECC.LE.0.01) GO TO 50
*
*       Include some diagnostic output.
      UPR2 = 0.0
      DO 5 K = 1,4
          UPR2 = UPR2 + UDOT(K,IPAIR)**2
    5 CONTINUE
*
      HI = (2.0*UPR2 - BODY(I))/R(IPAIR)
      EB = BODY(2*IPAIR-1)*BODY(2*IPAIR)*HI/BODY(I)
      ERR = (HI - H(IPAIR))/HI
      ZMU = BODY(2*IPAIR)*BODY(2*IPAIR-1)/BODY(I)
      DB = ZMU*(HI - H(IPAIR))
      IF (ABS(DB).GT.1.0D-08) THEN
      RA = R(IPAIR)/SEMI
      ERR2 = HI + 0.5*BODY(I)/SEMI
      ERR3 = UPR2 - 0.25*BODY(I)*(1.0 - ECC)
      IF (SEMI.LT.0.0) RA = R(IPAIR)
      WRITE (16,3)  TTOT, IPAIR, ECC, H(IPAIR),
     &              GAMMA(IPAIR), DB, ERR, ERR2,ERR3,TDOT2(IPAIR)
    3 FORMAT (' KSRECT:   T # E H G DB DH/H ',
     &                    F8.2,I5,F8.4,F7.1,F7.3,1P,5E10.1)
      CALL FLUSH(16)
      END IF
*
*       Initialize iteration counter for difficult case (SJA 10/97).
      ITER = 0
*
*       Form square regularized velocity for the explicit binding energy.
   10 UPR2 = 0.0
      DO 15 K = 1,4
          UPR2 = UPR2 + UDOT(K,IPAIR)**2
   15 CONTINUE
*
*       Form KS scaling factors from energy and angular momentum relation.
      A1 = 0.25D0*BODY(I)/UPR2
*       Solve for C1 from H = (2*U'*U'*C1**2 - M)/(U*U*C2**2) with C2 = 1/C1.
      A2 = A1**2 + 0.5D0*H(IPAIR)*R(IPAIR)/UPR2
*
*       Avoid negative round-off value on second try (NB! no change in CK).
      IF (ITER.EQ.2.AND.A2.LT.0.0) A2 = 0.0D0
*
*       Check for undefined case (circular orbit or eccentric anomaly = 90).
      IF (A2.GE.0.0D0) THEN
          IF (A1.LT.1.0) THEN
*       Choose square root sign from eccentric anomaly (e*cos(E) = 1 - R/a).
              C1 = SQRT(A1 + SQRT(A2))
          ELSE
              C1 = SQRT(A1 - SQRT(A2))
          END IF
          CK = 1.0
      ELSE
*       Adopt C1*C2 = CK for difficult case (Seppo's suggestion of 1991).
          C1 = 1.0
          CK = BODY(I)/SQRT(-8.0D0*H(IPAIR)*R(IPAIR)*UPR2)
          WRITE (6,20)  NAME(2*IPAIR-1), IPHASE, KSTAR(I), ECC,
     &                  R(IPAIR), GAMMA(IPAIR), DB, A2, CK-1.0
   20     FORMAT (' WARNING!    KSRECT    NM IPH K* E R G DB A2 CK-1 ',
     &                                    I6,2I4,F7.3,1P,5E10.2)
          ITER = ITER + 1
      END IF
*
*       Specify KS coordinate scaling from angular momentum conservation.
      C2 = CK/C1
*
*       Transform KS variables to yield the prescribed elements.
      R(IPAIR) = 0.0D0
*     UPR2 = 0.0D0
      TD2 = 0.0D0
      DO 25 K = 1,4
          U(K,IPAIR) = C2*U(K,IPAIR)
          UDOT(K,IPAIR) = C1*UDOT(K,IPAIR)
          U0(K,IPAIR) = U(K,IPAIR)
          R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
          TD2 = TD2 + U(K,IPAIR)*UDOT(K,IPAIR)
*         UPR2 = UPR2 + UDOT(K,IPAIR)**2
   25 CONTINUE
      TDOT2(IPAIR) = 2.0*TD2
*
*       Include diagnostic output.
*     HI = (2.0*UPR2 - BODY(I))/R(IPAIR)
*     WRITE (16,30)  IPAIR, IPHASE, KSTAR(I), R(IPAIR),
*    &               GAMMA(IPAIR), (HI-H(IPAIR))/HI
*  30 FORMAT (' KSRECT:    KS IPH K* R G DH/H ',3I4,1P,3E10.2)
*     CALL FLUSH(16)
*
*       Improve solution by second iteration in case of CK procedure.
      ITER = ITER + 1
      IF (ITER.EQ.2) GO TO 10
*
   50 RETURN
*
      END
