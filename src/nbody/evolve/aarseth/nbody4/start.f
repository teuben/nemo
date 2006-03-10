***
      SUBROUTINE START
*
*
*       Initialization of data & polynomials.
*       -------------------------------------
*
      INCLUDE 'common4.h'
      PARAMETER  (NS=12)
      EXTERNAL SCALE
*
*
*       Initialize global scalars, counters & useful constants.
      CALL ZERO
*
*       Read input parameters.
      CALL INPUT
*
*       Set initial conditions: BODY(I), X(K,I), XDOT(K,I); I=1,N & K=1,3.
      CALL DATA
*
*       Scale initial conditions to new units.
      CALL SCALE
*
*       Set total mass in case routines DATA & SCALE are not used.
      ZMASS = 0.0D0
      DO 10 I = 1,N
          ZMASS = ZMASS + BODY(I)
   10 CONTINUE
*
*       Define mean mass in scaled units and solar mass conversion factor.
      BODYM = ZMASS/FLOAT(N)
      IF (KZ(5).NE.3) THEN
          ZMBAR = ZMBAR/BODYM
      END IF
*
*       Introduce scaling factors DAYS, YRS, SU, RAU, SMU, TSTAR & VSTAR.
      CALL UNITS
*
*       Check option for external force.
      IF (KZ(14).GT.0.OR.KZ(14).EQ.-1) THEN 
          CALL XTRNL0
      END IF 
*
*       Check optional scaling to hot system.
      IF (KZ(29).GT.0) THEN
          CALL HOTSYS
      END IF
*
*       Check option for initial binaries.
      IF (KZ(8).EQ.1.OR.KZ(8).GE.3) THEN
          CALL BINPOP
      END IF
*
*       Include stable primordial triples.
      IF (KZ(18).GT.1.AND.KZ(8).GT.0) THEN
          CALL HIPOP
      END IF
*
*       Check optional initialization for tidal two-body capture.
      IF (KZ(27).GT.0) THEN
          CALL INTIDE
      END IF
*
*       Set sequential name, maximum & minimum mass and primary velocity.
      BODY1 = 0.0
      BODYN = 1.0
      DO 30 I = 1,N
          NAME(I) = I
          BODY1 = MAX(BODY1,BODY(I))
          BODYN = MIN(BODYN,BODY(I))
          DO 25 K = 1,3
              X0DOT(K,I) = XDOT(K,I)
   25     CONTINUE
          IF(I.LT.2*NBIN0) LIST(1,I) = 0
   30 CONTINUE
*
*       Initialize fixed block steps (40 levels).
      CALL IBLOCK
*
*       Generate initial coefficients for Stumpff functions.
      DO 40 I = 1,NS
          SCOEFF(i) = 1.0/((I+1)*(I+2))
   40 CONTINUE
*
*       Set optional stellar evolution parameters and initialize KTYPE.
      IF (KZ(19).GT.2.OR.KZ(19).LT.0) THEN
          CALL INSTAR
      ELSE
*       Define first quantized step < 100 yrs (minimum interval for MDOT).
          DT = 1.0d-04/TSCALE
          CALL STEPK(DT,DTN)
          IF (DTN*TSCALE.LT.100.0) DTN = 2.0*DTN
          STEPX = DTN
      END IF
*
*       Initialize optional cloud parameters (NB! negative option).
      IF (KZ(12).LT.0) THEN
          CALL CLOUD0
      END IF
*
*        Initialize GRAPE, force and first derivative and time-steps. 
      CALL FPOLY0
*
*       Regularize any hard primordial binaries (assume sequential ordering).
      IF (NBIN0.GT.0) THEN
*       Change indicator temporarily for FPOLYI with GRAPE (then restore).
          IPHASE = 1
          DO 50 IPAIR = 1,NBIN0
              ICOMP = 2*IPAIR - 1
              JCOMP = 2*IPAIR
              RIJ2 = 0.0
*       Include standard distance criterion.
              DO 45 K = 1,3
                  RIJ2 = RIJ2 + (X(K,ICOMP) - X(K,JCOMP))**2
   45         CONTINUE
              IF (RIJ2.LT.RMIN**2) THEN
                  CALL KSREG
              END IF
   50     CONTINUE
*
*       Make new perturber list if any member is a KS component.
          DO 60 JPAIR = 1,NPAIRS
              J1 = 2*JPAIR - 1
              NNB1 = LIST(1,J1) + 1
              DO 55 L = 2,NNB1
                  IF (LIST(L,J1).LT.IFIRST) THEN
                      CALL KSLIST(JPAIR)
                      GO TO 60
                  END IF
   55         CONTINUE
   60     CONTINUE
*         IPHASE = -2
      END IF
*
*       Initialize the time-step list used to find next body (Hermite only).
      DTLIST = 100.0
      DO 70 I = IFIRST,NTOT
          DTLIST = MIN(DTLIST,STEP(I))
   70 CONTINUE
*
*       Set initial time-step list interval twice the smallest step.
      IF (N.GT.10) THEN
          DTLIST = 2.0*DTLIST
      ELSE
          DTLIST = TCRIT
      END IF
      NNB = 1
   80 TLIST = TLIST + DTLIST
*
*       Select all members due in the interval (0,TLIST).
      DO 90 J = IFIRST,NTOT
          IF (T0(J) + STEP(J).LT.TLIST) THEN
              NNB = NNB + 1
              NLIST(NNB) = J
              IF(NNB.GE.LMAX-3) GOTO 92
          END IF
   90 CONTINUE
*
*       Check whether membership range is acceptable. 
   92 IF (NNB.EQ.1) GO TO 80
*
*       Reduce new DTLIST to prevent early crowding and set membership.
      DTLIST = 0.2*DTLIST
      NLIST(1) = NNB - 1
*
      RETURN
      END
***
