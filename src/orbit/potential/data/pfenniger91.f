C
C      23-feb-2002   derived from pfennifer84
C                    - added 2nd Miyamoto disk
C                    
C
C

c
c-----------------------------------------------------------------------
c
      SUBROUTINE inipotential(npar,par,name)
      IMPLICIT NONE
      INTEGER npar
      DOUBLE PRECISION par(*)
      CHARACTER name*(*)
c
      DOUBLE PRECISION gmb,a,aob,aoc

      DOUBLE PRECISION gmd1,gmd2,ad1,ad2,bd,bd2
      COMMON /miya/gmd1,gmd2,ad1,ad2,bd,bd2



      gmb = 0.10D0
      gmd1= 0.10D0
      gmd2= 0.80D0
      a   = 6.00D0
      aob = 4.00D0
      aoc = 10.0D0
      ad1 = 0.14D0
      ad2 = 3.00D0
      bd  = 1.00D0
      IF (npar.GE.2) gmb = par(2)

      CALL inpo3 (a,a/aob,a/aoc,gmb)

      bd2 = bd*bd

c     return the pattern speed of the standard model (R_cr = a)
      par(1) = 0.05471

      RETURN
      END
c
c-----------------------------------------------------------------------
c
      SUBROUTINE potential(ndim, pos, acc, pot, time)
      IMPLICIT NONE
      INTEGER ndim
      DOUBLE PRECISION pos(*), acc(*), pot, time
c
      DOUBLE PRECISION r2,q2,s2,q,tmp
c
      DOUBLE PRECISION gmd1,gmd2,ad1,ad2,bd,bd2
      COMMON /miya/gmd1,gmd2,ad1,ad2,bd,bd2

c  bar
      CALL pbar3(pos(1),pos(2),pos(3),pot)
      CALL fbar3(pos(1),pos(2),pos(3),acc(1),acc(2),acc(3))
c  shared miyamoto disks variables
      r2 = pos(1)*pos(1)+pos(2)*pos(2)
      q2 = pos(3)*pos(3) + bd2
      q = sqrt(q2)

c  miyamoto disk 1 ('bulge')
      s2 = r2 + (ad1+q)*(ad1+q)
      pot = pot - gmd1/sqrt(s2)
      tmp = pot/s2
      acc(1) = acc(1) + tmp*pos(1)
      acc(2) = acc(2) + tmp*pos(2)
      acc(3) = acc(3) + tmp*pos(3)*(ad1+q)/q
c  miyamoto disk 2 ('disk')
      s2 = r2 + (ad2+q)*(ad2+q)
      pot = pot - gmd2/sqrt(s2)
      tmp = pot/s2
      acc(1) = acc(1) + tmp*pos(1)
      acc(2) = acc(2) + tmp*pos(2)
      acc(3) = acc(3) + tmp*pos(3)*(ad2+q)/q

      END

!-----------------------------------------------------------------------

      SUBROUTINE INPO3 (A,B,C,GMB)

!     Initialisation of the Ferrers n=2 potential
!     Semi-axes A,B,C and mass GMB
!     Restriction : A > B > C >= 0

!     D. Pfenniger, Observatoire de Geneve, 1984, cleaned 2002
!     Reference: Pfenniger D., Astron. Astrophys. 134, 373 (1984) (Appendix A)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /BARR3/
     &   A2,B2,C2,UA2,UB2,UC2,CTE,RHOC,CT6,XK,UA2B2,UA2C2,UB2C2,SUA2C2,
     &   V000,V100,V010,V001,V110,V101,V011,V200,V020,V002,
     &   V111,V210,V201,V120,V021,V102,V012,V300,V030,V003
!
      PARAMETER (US3 = 1D0/3D0, PI=3.141592653589793D0)

!     Constant terms

      A2 = A*A
      B2 = B*B
      C2 = C*C

      UA2 = 1D0/A2
      UB2 = 1D0/B2
      UC2 = 1D0/C2

      CT2 = GMB*105D0/32D0
      CTE = -2D0*CT2
      CT6 = CTE/6D0

      RHOC = CT2/(PI*A*B*C)

      UA2B2 = 1D0/(A2 - B2)
      UA2C2 = 1D0/(A2 - C2)
      UB2C2 = 1D0/(B2 - C2)

      SUA2C2 = DSQRT(UA2C2)

      PHI = ASIN(SQRT(1D0 - C2*UA2))
      XK = SQRT(UA2C2*(A2-B2))
      CALL ELINT(PHI,XK,F,E)
      D2 = 2D0*SQRT(UA2*UB2*UC2)

!     Wijk inside the ellipsoid are constant!

      V000 = 2D0*F*SUA2C2
      V100 = 2D0*(F-E)*UA2B2*SUA2C2
      V001 = (D2*B2 - 2D0*E*SUA2C2)*UB2C2
      V010 = D2 - V100 - V001

      V110 = (V010 - V100)*UA2B2
      V101 = (V001 - V100)*UA2C2
      V011 = (V001 - V010)*UB2C2

      V200 = (D2*UA2 - V110 - V101)*US3
      V020 = (D2*UB2 - V110 - V011)*US3
      V002 = (D2*UC2 - V011 - V101)*US3

      V111 = (V011 - V110)*UA2C2
      V210 = (V110 - V200)*UA2B2
      V201 = (V101 - V200)*UA2C2
      V120 = (V020 - V110)*UA2B2
      V021 = (V011 - V020)*UB2C2
      V102 = (V002 - V101)*UA2C2
      V012 = (V002 - V011)*UB2C2

      V300 = (D2*UA2*UA2 - V210 - V201)*.2D0
      V030 = (D2*UB2*UB2 - V120 - V021)*.2D0
      V003 = (D2*UC2*UC2 - V102 - V012)*.2D0

      RETURN
      END

!-----------------------------------------------------------------------

      SUBROUTINE RBAR3 (X,Y,Z,RHO)

!     Mass density RHO of the Ferrers n=2 potential at X, Y, Z
!     Semi-axes A,B,C and mass GMB

!     D. Pfenniger, Observatoire de Geneve, 1984, cleaned 2002
!     Reference: Pfenniger D., Astron. Astrophys. 134, 373 (1984) (Appendix A)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!***********************************************************************
!     The following COMMON parameters must be initialized by INPO3 
!     before the first use of this routine
!***********************************************************************

      COMMON /BARR3/
     &   A2,B2,C2,UA2,UB2,UC2,CTE,RHOC,CT6,XK,UA2B2,UA2C2,UB2C2,SUA2C2,
     &   V000,V100,V010,V001,V110,V101,V011,V200,V020,V002,
     &   V111,V210,V201,V120,V021,V102,V012,V300,V030,V003

      XM2 = X*X*UA2+Y*Y*UB2+Z*Z*UC2

      IF (XM2 .LT. 1D0) THEN

         RHO = RHOC * (1D0 - XM2)**2

      ELSE

         RHO = 0D0

      ENDIF

      RETURN
      END

!-----------------------------------------------------------------------

      SUBROUTINE FBAR3 (X,Y,Z,FX,FY,FZ)

!     Forces FX, FY, FZ of the Ferrers n=2 potential at X, Y, Z
!     Semi-axes A,B,C and mass GMB
!     Restriction : A > B > C >= 0
!     Warning : Ferrers' polynomial formulation produces 
!               large numerical errors at distances >> A,B,C
!               or when 2 semi-axes are very close

!     D. Pfenniger, Observatoire de Geneve, 1984, cleaned 2002
!     Reference: Pfenniger D., Astron. Astrophys. 134, 373 (1984) (Appendix A)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!***********************************************************************
!     The following COMMON parameters must be initialized by INPO3 
!     before the first use of this routine
!***********************************************************************

      COMMON /BARR3/
     &   A2,B2,C2,UA2,UB2,UC2,CTE,RHOC,CT6,XK,UA2B2,UA2C2,UB2C2,SUA2C2,
     &   V000,V100,V010,V001,V110,V101,V011,V200,V020,V002,
     &   V111,V210,V201,V120,V021,V102,V012,V300,V030,V003

      PARAMETER (US3 = 1D0/3D0)

      X2 = X*X
      Y2 = Y*Y
      Z2 = Z*Z

      IF (X2*UA2+Y2*UB2+Z2*UC2 .LE. 1D0) THEN

         FX = CTE * X * (V100+X2*(X2*V300+2D0*(Y2*V210-V200))
     &                       +Y2*(Y2*V120+2D0*(Z2*V111-V110))
     &                       +Z2*(Z2*V102+2D0*(X2*V201-V101)))
         FY = CTE * Y * (V010+X2*(X2*V210+2D0*(Y2*V120-V110))
     &                       +Y2*(Y2*V030+2D0*(Z2*V021-V020))
     &                       +Z2*(Z2*V012+2D0*(X2*V111-V011)))
         FZ = CTE * Z * (V001+X2*(X2*V201+2D0*(Y2*V111-V101))
     &                       +Y2*(Y2*V021+2D0*(Z2*V012-V011))
     &                       +Z2*(Z2*V003+2D0*(X2*V102-V002)))

      ELSE

         XL = XLMBD(X2,Y2,Z2,A2,B2,C2)
         UA3 = 1D0/(A2+XL)
         B3 = B2+XL
         UB3 = 1D0/B3
         UC3 = 1D0/(C2+XL)
         PHI= ASIN(SQRT(UA3/UA2C2))
         CALL ELINT(PHI,XK,F,E)
         D2 = 2D0*SQRT(UA3*UB3*UC3)

         W100 = 2D0*(F-E)*UA2B2*SUA2C2
         W001 = (D2*B3 - 2D0*E*SUA2C2)*UB2C2
         W010 = D2 - W100 - W001
     
         W110 = (W010 - W100)*UA2B2
         W101 = (W001 - W100)*UA2C2
         W011 = (W001 - W010)*UB2C2
     
         W200 = (D2*UA3 - W110 - W101)*US3
         W020 = (D2*UB3 - W110 - W011)*US3
         W002 = (D2*UC3 - W011 - W101)*US3

         W111 = (W011 - W110)*UA2C2
         W210 = (W110 - W200)*UA2B2
         W201 = (W101 - W200)*UA2C2
         W120 = (W020 - W110)*UA2B2
         W021 = (W011 - W020)*UB2C2
         W102 = (W002 - W101)*UA2C2
         W012 = (W002 - W011)*UB2C2

         W300 = (D2*UA3*UA3 - W210 - W201)*.2D0
         W030 = (D2*UB3*UB3 - W120 - W021)*.2D0
         W003 = (D2*UC3*UC3 - W102 - W012)*.2D0

         FX = CTE * X * (W100+X2*(X2*W300+2D0*(Y2*W210-W200))
     &                       +Y2*(Y2*W120+2D0*(Z2*W111-W110))
     &                       +Z2*(Z2*W102+2D0*(X2*W201-W101)))
         FY = CTE * Y * (W010+X2*(X2*W210+2D0*(Y2*W120-W110))
     &                       +Y2*(Y2*W030+2D0*(Z2*W021-W020))
     &                       +Z2*(Z2*W012+2D0*(X2*W111-W011)))
         FZ = CTE * Z * (W001+X2*(X2*W201+2D0*(Y2*W111-W101))
     &                       +Y2*(Y2*W021+2D0*(Z2*W012-W011))
     &                       +Z2*(Z2*W003+2D0*(X2*W102-W002)))

      END IF

      RETURN
      END

!-----------------------------------------------------------------------

      SUBROUTINE PBAR3 (X,Y,Z,PHI)

!     Potential PHI of the Ferrers n=2 potential at X, Y, Z
!     Semi-axes A,B,C and mass GMB
!     Restriction : A > B > C >= 0
!     Warning : Ferrers' polynomial formulation produces 
!               large numerical errors at large distances >> A,B,C
!               or when 2 semi-axes are very close

!     D. Pfenniger, Observatoire de Geneve, 1984, cleaned 2002
!     Reference: Pfenniger D., Astron. Astrophys. 134, 373 (1984) (Appendix A)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!***********************************************************************
!     The following COMMON parameters must be initialized by INPO3 
!     before the first use of this routine
!***********************************************************************

      COMMON /BARR3/
     &   A2,B2,C2,UA2,UB2,UC2,CTE,RHOC,CT6,XK,UA2B2,UA2C2,UB2C2,SUA2C2,
     &   V000,V100,V010,V001,V110,V101,V011,V200,V020,V002,
     &   V111,V210,V201,V120,V021,V102,V012,V300,V030,V003

      PARAMETER (US3=1D0/3D0)

      X2 = X*X
      Y2 = Y*Y
      Z2 = Z*Z

      IF (X2*UA2+Y2*UB2+Z2*UC2 .LE. 1D0) THEN

        PHI = CT6 * 
     &          (V000 - 6D0*X2*Y2*Z2*V111 +
     &             X2*(X2*(3D0*V200-X2*V300) + 
     &                     3D0*(Y2*(V110+V110-Y2*V120-X2*V210)-V100)) +
     &             Y2*(Y2*(3D0*V020-Y2*V030) + 
     &                     3D0*(Z2*(V011+V011-Z2*V012-Y2*V021)-V010)) +
     &             Z2*(Z2*(3D0*V002-Z2*V003) + 
     &                     3D0*(X2*(V101+V101-X2*V201-Z2*V102)-V001)))

      ELSE

         XL = XLMBD(X2,Y2,Z2,A2,B2,C2)
         UA3 = 1D0/(A2+XL)
         B3 = B2+XL
         UB3 = 1D0/B3
         UC3 = 1D0/(C2+XL)
         PHI= ASIN(SQRT(UA3/UA2C2))
         CALL ELINT(PHI,XK,F,E)
         D2 = 2D0*SQRT(UA3*UB3*UC3)

         W100 = 2D0*(F-E)*UA2B2*SUA2C2
         W001 = (D2*B3 - 2D0*E*SUA2C2)*UB2C2
         W010 = D2 - W100 - W001

         W110 = (W010 - W100)*UA2B2
         W101 = (W001 - W100)*UA2C2
         W011 = (W001 - W010)*UB2C2

         W200 = (D2*UA3 - W110 - W101)*US3
         W020 = (D2*UB3 - W110 - W011)*US3
         W002 = (D2*UC3 - W011 - W101)*US3

         W111 = (W011 - W110)*UA2C2
         W210 = (W110 - W200)*UA2B2
         W201 = (W101 - W200)*UA2C2
         W120 = (W020 - W110)*UA2B2
         W021 = (W011 - W020)*UB2C2
         W102 = (W002 - W101)*UA2C2
         W012 = (W002 - W011)*UB2C2

         W300 = (D2*UA3*UA3 - W210 - W201)*.2D0
         W030 = (D2*UB3*UB3 - W120 - W021)*.2D0
         W003 = (D2*UC3*UC3 - W102 - W012)*.2D0

         PHI = CT6 * 
     &          (2D0*F*SUA2C2 - 6D0*X2*Y2*Z2*W111 +
     &             X2*(X2*(3D0*W200-X2*W300) + 
     &                     3D0*(Y2*(W110+W110-Y2*W120-X2*W210)-W100)) +
     &             Y2*(Y2*(3D0*W020-Y2*W030) + 
     &                     3D0*(Z2*(W011+W011-Z2*W012-Y2*W021)-W010)) +
     &             Z2*(Z2*(3D0*W002-Z2*W003) + 
     &                     3D0*(X2*(W101+W101-X2*W201-Z2*W102)-W001)))
      END IF

      RETURN
      END

!-----------------------------------------------------------------------

      FUNCTION XLMBD (X,Y,Z,A,B,C)

!     Resolution of the equation for L:

!     X/(A+L) + Y/(B+L) + Z/(C+L) = 1

!     Restriction :     X/A + Y/B + Z/C   >   1       ,

!***********************************************************************
!     The input parameters are assumed to be positive
!***********************************************************************

!     D. Pfenniger, Observatoire de Geneve, 1984, cleaned 2002

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER (US3 = 1D0/3D0)

      C2 = (A+B+C-X-Y-Z)*US3
      C1 = A*(B-Y-Z) + B*(C-X-Z) + C*(A-X-Y)
      C0 = A*B*(C-Z) - (A*Y+X*B)*C
      P = C1*US3-C2*C2
      Q = C2**3 + (C0-C1*C2)*.5D0
      DE = P**3+Q*Q
      IF (DE .LT. 0D0) THEN
         R = SQRT(-P)
         XLMBD = 2D0*R*COS( ACOS(-Q/R**3)*US3 ) - C2
      ELSE
         R = -Q+SQRT(DE)
         IF (R .LT. 0D0) THEN
            U = -(-R)**US3
         ELSE
            U = R**US3
         END IF
         R = -Q-SQRT(DE)
         IF (R .LT. 0D0) THEN
            V = -(-R)**US3
         ELSE
            V = R**US3
         END IF
         XLMBD = U + V - C2
      END IF

      RETURN
      END

!-----------------------------------------------------------------------

      SUBROUTINE ELINT(PHI,K,F,E)

!     Incomplete elliptic integrals F & E

!     Method : descending Landen transformation
!              Ref.: Abramowitz & Stegun  17.5 p. 597

!     Input parameters : PHI, K : amplitude angle and module
!                        0 <= PHI < PI/2   ;   0 <= K < 1
!     Output values    : F, E    : incomplete elliptic integrals of the
!                                  first and second kind

!     D. Pfenniger, Observatoire de Geneve, 1984, cleaned 2002

      IMPLICIT DOUBLE PRECISION (A-Z)

      PARAMETER (PI=3.141592653589793D0, TPI=2D0*PI)

      A = 1D0
      C = K
      S = C * C
      B = SQRT(1D0 - S)
      P = PHI
      P2 = P
      T = 0D0
      I = 1D0

 999  CONTINUE
         I = I + I
         P = P + INT(P/PI + 0.5D0)*PI + ATAN(TAN(P2) * B/A)
         P2 = MOD(P,TPI)
         C = 0.5D0 * (A - B)
         A1= 0.5D0 * (A + B)
         B = SQRT(A * B)
         A = A1
         S = S + C * C * I  
         T = T + C * SIN(P2)
      IF (C .GE. 1D-15) GOTO 999

      F = P / (I * A)
      E = T + F * (1D0 - 0.5D0 * S)

      RETURN
      END
