C
C   History:
C
C   Taken from Lia:
C       Date: Thu, 21 Oct 1993 08:22:15 EDT
C       From: lia@OBMARA.CNRS-MRS.FR
C       Subject: RE: Bar Potentials
C
C   21-oct-93   Put into potential(5NEMO) format                PJT
C		Special case for (0,0,0) -- doesn't work because of LOG()
c   23-apr-96   removed fortran specific things for loadobj()       pjt
c   20-jun-01   geez, SGI doesn't like ctrl-L in the file....       pjt
C
C   ToDO:
C	code DACSH()
C=======================================================================
      SUBROUTINE inipotential(npar, par, name)
C
      IMPLICIT NONE
      INTEGER npar
      REAL*8 par(*)
      CHARACTER name*(*)
C-----------------------------------------------------------------------
C INITIALIZE THE POTENTIAL AND FORCE CALCULATIONS (after Lia's INPOT)
C
      COMMON/ AREAD/FCMASS,AXIRAT,QUAD,CENDEN,RADLAG, INDEX
      REAL*8 FCMASS, AXIRAT, QUAD, CENDEN, RADLAG
      INTEGER INDEX
C 
      COMMON/AXISY/C1BUL,CONST, CORAD,VT2,VTDR2, RTI2
      REAL*8 C1BUL, CONST, CORAD, VT2, VTDR2, RTI2
C
      COMMON /COM1/ OMEGS, A2, B2, C2, E, NEQ
      INTEGER NEQ
      REAL*8 OMEGS, A2, B2, C2, E
C
      COMMON /COM2/ H, EPS, IOUT, IORB, ISIGN, ICHOOSE
      INTEGER IOUT, IORB, ISIGN, ICHOOSE
      REAL*8 H, EPS
C
      COMMON/PROLA/ELIPM,A,B,C,AAIN,CCIN,AAI,CCO,ACONST
      REAL*8 ELIPM, A, B, C, ACONST, CCO, CCIN, AAI, AAIN
C
      COMMON/PROLA2/COEFF,COEFP,EP,EP2,SHO,CHO,THO,PSIO,CSHO,CTHO
      REAL*8 COEFF, COEFP, EP, EP2, SHO, CHO, THO, PSIO, CSHO, CTHO
C
      INTEGER NITER
      REAL*8 FUN1, ASINH, Z, DACSH
      REAL*8 RTMAX, VTMAX
      REAL*8 AAXIS, BAXIS, BAXIS2, QUADM, BARMAS, BARDEN, BULMAS, BULDEN
      REAL*8 CENMAS
      REAL*8 RQ, FMAS, AMAS, DMAS, DDMAS, QQ, AT, AX, AY, AZ
      REAL*8 PIGROC
C
      REAL*8 PI, FOURPI, GRAVC
      PARAMETER ( PI = 3.141592654, FOURPI = 4. * PI, GRAVC = 4.29569)
C
      FUN1 (Z) = SQRT( Z * Z + 1.)
      ASINH (Z) = LOG( Z + SQRT( 1.0 + Z * Z))
c
c set constants from the initpotential par-list
c   1 = pattern speed (returned)
c   2 = ignored in this routine  (should be 1 though)
c   3 = index (0,1,2)
c   4 = radlag
c
c      COMMON/ AREAD/FCMASS,AXIRAT,QUAD,CENDEN,RADLAG, INDEX
c      REAL*8 FCMASS, AXIRAT, QUAD, CENDEN, RADLAG
c      INTEGER INDEX

c first, set defaults:
      index = 1
      axirat = 2.5
      radlag = 6.0
      quad = 45000.0
      cenden = 24000.0

c override the first npar:
c but omega (par(1)) and amode (par(2) not used in this version)
      IF (npar.GE.3) index = par(3)
      IF (npar.GE.4) axirat = par(4)
      IF (npar.GE.5) radlag = par(5)
      IF (npar.GE.6) quad = par(6)
      IF (npar.GE.7) cenden = par(7)
C
C
C  GRAVC IN UNITS OF KPC**3 GYR**-2 1.E6MSUN**-1
C
      RTMAX = 20.
      RTI2 = 1. / (RTMAX * RTMAX)
      VTMAX = 164.204
      VT2 = VTMAX * VTMAX
      VTDR2= VT2 * RTI2
C
      fcmass = 1.
      AAXIS = 5.0
      CENMAS = 4.87333E4
      CENMAS = CENMAS * FCMASS
      BAXIS = AAXIS / AXIRAT
      BAXIS2 = BAXIS ** 2
      IF( INDEX .EQ. 0) QUADM = CENMAS * (AAXIS ** 2 - BAXIS2) / 5.0
      IF( INDEX .EQ. 1) QUADM = CENMAS * (AAXIS ** 2 - BAXIS2) / 7.0
      QUAD = MIN(QUAD, QUADM)
      BARMAS = CENMAS * QUAD / QUADM
      IF(INDEX .EQ. 0) BARDEN = BARMAS * 0.2387324    / (AAXIS * BAXIS2)
      IF(INDEX .EQ. 1) BARDEN = BARMAS * 0.5968310366 / (AAXIS * BAXIS2)
      BULMAS = CENMAS - BARMAS
      BULDEN = 3.0E-3 * BULMAS / FOURPI
      CENDEN = MAX(CENDEN, BARDEN + BULDEN)
      BULDEN = CENDEN - BARDEN
      CORAD = 1.0
      BULMAS = CENMAS - BARMAS
      NITER = 0
   10 RQ = 10.0 / CORAD
      AMAS = FOURPI * BULDEN * CORAD * CORAD * CORAD *
     1       (ASINH(RQ) - RQ / FUN1(RQ))
      DMAS = BULMAS - AMAS
      FMAS = DMAS / BULMAS
      IF (ABS(FMAS) .LE. 1.0E-6) GO TO 11
      IF (NITER .GE. 20) GO TO 13
      NITER = NITER + 1
      DDMAS = FOURPI * BULDEN * CORAD * CORAD * (3.0 * ASINH(RQ) -
     1        RQ * (3.0 + 4.0 * RQ * RQ) / (FUN1(RQ) ** 3))
      CORAD = MAX(0.01 * CORAD, CORAD + DMAS / DDMAS)
      GO TO 10
   11 CONTINUE
      C1BUL = GRAVC * BULDEN * FOURPI * CORAD * CORAD
      CONST=12.56637061*GRAVC*CORAD*CORAD*CORAD*BULDEN
C
C     BAR
C
      A=AAXIS
      B=BAXIS
      AAIN = A * A
      AAI = 1. / AAIN
      C = BAXIS
      CCIN = C * C
      CCO = 1. / CCIN
      ELIPM = BARMAS
      ACONST = ELIPM * 4.027209375 / (A * C * C)
C
      IF( INDEX .EQ. 0) THEN
         PIGROC = PI * GRAVC * BARDEN
         EP = DSQRT( AAXIS**2 - BAXIS2)
         EP2 = EP * EP
c oops, need to define DACSH here for INDEX=0 case
cpjt         PSIO = DACSH( AAXIS / BAXIS)
c
c	 STOP
c
         SHO = DSINH (PSIO) 
         CHO = AAXIS / BAXIS
         THO = SHO / CHO 
         CSHO = SHO ** 2 / EP ** 2
         CTHO = THO ** 2 / EP ** 2
         COEFF = 2. * PIGROC / (THO * SHO * SHO)
         COEFP = PIGROC *BAXIS2 / THO
      END IF
C
C
C     CALCULATE PATTERN SPEED
C
      RQ = RADLAG / CORAD
      AMAS = FOURPI * BULDEN * CORAD * CORAD * CORAD *
     1       (ASINH(RQ) - RQ / FUN1(RQ))
      QQ = GRAVC * AMAS / RADLAG
C  BAR NOW LIES ON THE Y AXIS. SO INVERT X AND Y TO HAVE THE SAME LAGRAN
C  POINT AS IN THE HYDRO CODE
      IF( INDEX .EQ. 0) CALL FORBA0( 0.D0, RADLAG, 0.D0,AY,AX,AZ)
      IF( INDEX .EQ. 1) CALL FORBA ( 0.D0, RADLAG, 0.D0,AY,AX,AZ)
      QQ = QQ - AX * RADLAG
      AT = RADLAG * VT2 *
     1     (3. / (1. + 2. * RADLAG * RADLAG * RTI2)) ** 1.5
      AT = AT * RTI2
      QQ = (QQ + AT * RADLAG) / (RADLAG * RADLAG)
      OMEGS = SQRT(QQ)
      par(1) = omegs
      RETURN
13    CONTINUE
c13    PRINT 613, CORAD,FMAS,NITER
c613   FORMAT(' NO CONVERGENCE FOR CORE RADIUS; RADIUS=',1PE12.5,
c     &       ' FRACTION MASS DIFFERENCE =',E12.5,' NITER =',I5)
c      STOP
      END
C-----------------------------------------------------------------------      
      SUBROUTINE potential(ndim, pos, acc, pot, time)
      INTEGER ndim
      REAL*8 pos(*), acc(*), pot, time
   
      acc(1) = 0.0
      acc(2) = 0.0
      acc(3) = 0.0
      pot = 0.0
c			note exchanged X and Y, since their bar
c			is along the Y axis, ours along the X
      CALL FISOB(pos(2), pos(1), pos(3), acc(2), acc(1), acc(3))
      CALL PISOB(pos(2), pos(1), pos(3), pot)
      END
C-----------------------------------------------------------------------
      SUBROUTINE FISOB( X, Y, Z, FX, FY, FZ)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     CALCUL DES FORCES DANS LE POTENTIEL (PISOB)
C     ROUTINE APPELLEE PAR DISOB (DERIVEES POUR RK78)
C
      IMPLICIT NONE
      REAL*8 X, Y, Z, FX, FY, FZ
C
      COMMON/ AREAD/FCMASS,AXIRAT,QUAD,CENDEN,RADLAG, INDEX
      REAL*8 FCMASS, AXIRAT, QUAD, CENDEN, RADLAG
      INTEGER INDEX
C 
      REAL*8 FAX, FAY, FAZ, FBX, FBY, FBZ
      CALL FORAX( X, Y, Z, FAX, FAY, FAZ) 
      IF( QUAD .GT. 1.E-9) THEN
         IF( INDEX .EQ. 0) CALL FORBA0( X, Y, Z, FBX, FBY, FBZ) 
         IF( INDEX .EQ. 1) CALL FORBA ( X, Y, Z, FBX, FBY, FBZ) 
         FX = FAX + FBX
         FY = FAY + FBY
         FZ = 0.
      ELSE
         FX = FAX
         FY = FAY
         FZ = 0.
      END IF
C 
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FORAX( X, Y, Z, FX, FY, FZ)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C FORCES DUES A UN POTENTIEL AXISYMETRIQUE
C
      IMPLICIT NONE
      REAL*8 X, Y, Z, FX, FY, FZ
C
      COMMON/AXISY/C1BUL,CONST,CORAD,VT2,VTDR2,RTI2
      REAL*8 C1BUL, CONST, CORAD, VT2, VTDR2, RTI2
C
      REAL*8 FUN1, ASINH, ZZ
      REAL*8 R2, R, RQ, QQ
      FUN1 (ZZ) = SQRT( ZZ * ZZ + 1.)
      ASINH (ZZ) = LOG( ZZ + FUN1 (ZZ))
C
      R2 = X * X + Y * Y
      IF (R2.GT.0.0D0) THEN
         R = SQRT (R2)
         RQ = R / CORAD
         QQ = CONST * (ASINH(RQ) - RQ / FUN1(RQ))/(R * R * R) +
     1        VTDR2 * (3./ (1. + 2. * R * R * RTI2)) ** 1.5
         FX = - X * QQ
         FY = - Y * QQ
      ELSE
         FX = 0.
         FY = 0.
      ENDIF
      FZ = 0.
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FORBA( XINPUT, YINPUT, Z, AXOUT, AYOUT, AZ)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C FORCES FROM A FERRER'S BAR WITH INDEX = 1
C     SUBROUTINE PROL6 (ELIPM,A,B,C,X,Y,Z,AX,AY,AZ,EPT)
C      IMPLICIT real*8 (A-Z)
      IMPLICIT NONE
      REAL*8 XINPUT, YINPUT, Z, AXOUT, AYOUT, AZ
C
      COMMON/PROLA/ELIPM,A,B,C,AAIN,CCIN,AAI,CCO,ACONST
      REAL*8 ELIPM, A, B, C, ACONST, CCO, CCIN, AAI, AAIN
C
      REAL*8 X, Y, AA, CC, XX, YY, RR, GI, SIDE, BK, CK
      REAL*8 KAPPA, EM, EE, E, LNE, A3, A1, I, A33, A13, A11
      REAL*8 X4, Y4, AX, AY
C     CAREFUL : IN HYDRO PROGRAM THE BAR IS ALONG THE X AXIS
C     IN THE ORBIT CALCULATIONS ALONG THE Y AXIS
C
      X = YINPUT
      Y = XINPUT
      AA = AAIN
      CC = CCIN
      XX = X*X
      YY = Y*Y
      RR = XX+YY
      GI = 1.
      SIDE = XX/AA + (YY)/CC - 1.
      IF (SIDE .LE. 0.) GO TO 10
         BK = RR-AA-CC
         CK = AA*CC*SIDE
         KAPPA = (BK+SQRT(BK*BK+4.*CK))/2.
         AA = AA + KAPPA
         CC = CC + KAPPA
         GI = A*C*C/SQRT(AA)/CC
   10 EM = CC/AA
      EE = 1. - EM
      E = SQRT(EE)
      LNE = LOG((1.+E)/(1.-E))
C
C     CALCULATE INDEX SYMBOLS IN EASIEST WAY
C
      A3 = GI*(1.-EM*LNE/(E + E))/EE
      A1 = 2.*(GI - A3)
      I = AA*A1 + 2.*CC*A3
 
      A33 = GI*(4.*EE/EM - 6. + 3.*EM*LNE/E)/(8.*AA*EE*EE)
      A13 = 2.*(GI/CC - 2.*A33)
      A11 = 2.*(GI/AA - A13)/3.
 
      X4 = A11*XX + A13*(YY)
      Y4 = A13*XX + A33*(YY)
C     Z4 = Y4
 
C     EPT = I - 2.*(A1*XX + A3*(YY)) +(XX*X4+YY*Y4)
 
      AX = 4. * X * (X4 - A1)
 
      AY = 4. * Y * (Y4 - A3)
 
      ACONST = ELIPM*4.027209375/(A*C*C)
C     EPT = -EPT * ACONST
      AX = AX * ACONST
      AY = AY * ACONST
      AZ = 0.
      AXOUT=AY
      AYOUT=AX
      RETURN
      END
C-------------------------------------------------------------------------------
      SUBROUTINE FORBA0( X, Y, Z, FX, FY, FZ)
C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X, Y, Z, FX, FY, FZ
C
      COMMON/PROLA2/COEFF,COEFP,EP,EP2,SHO,CHO,THO,PSIO,CSHO,CTHO
      REAL*8 COEFF, COEFP, EP, EP2, SHO, CHO, THO, PSIO, CSHO, CTHO
C
      REAL*8 R2, Y2, XM, PSI, SH, CH, TH, AA, W10, W11, W20
      R2 = X * X
      Y2 = Y * Y
      XM = SQRT( R2 * CSHO + Y2 * CTHO)
      IF( XM .LE. 1.) THEN
         PSI = PSIO
         SH = SHO
         CH = CHO
         TH = THO
      ELSE
         AA = Y2 + R2 - EP2 
         IF( R2 * EP2 .LE. 1.D-10 * AA) THEN
            SH = EP / SQRT( Y2 - EP2)
         ELSE
            SH = SQRT(( SQRT( AA * AA + 4. * R2 * EP2) - AA) / (R2+R2))
         END IF
         PSI = LOG( SH + SQRT( SH ** 2 + 1.)) 
         CH = DCOSH (PSI) 
         TH = SH / CH
      END IF
      W10 = PSI + PSI 
      W20 = .5 * W10 - CH * SH
      W11 = TH + TH - W10 
      FX = COEFF * X * W20
      FY = COEFF * Y * W11
      FZ = 0.
      RETURN
      END 
C-----------------------------------------------------------------------
      SUBROUTINE PISOB( X, Y, Z, PHI)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     POTENTIEL AXISYMETRIQUE + BARRE PROLATE
      IMPLICIT NONE
      REAL*8 X, Y, Z, PHI
C
      COMMON/ AREAD/FCMASS,AXIRAT,QUAD,CENDEN,RADLAG, INDEX
      REAL*8 FCMASS, AXIRAT, QUAD, CENDEN, RADLAG
      INTEGER INDEX
C
      REAL*8 PAX, PBA
C
      CALL POTAX( X, Y, Z, PAX) 
      IF ( INDEX .EQ. 0) CALL POTBA0 ( X, Y, Z, PBA)
      IF ( INDEX .EQ. 1) CALL POTBA  ( X, Y, Z, PBA)
      PHI = PAX + PBA 
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE POTAX( X, Y, Z, POT)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C POTENTIEL AXISYMETRIQUE DU BULBE ET DU DISQUE
      IMPLICIT NONE
      REAL*8 X, Y, Z, POT
C
      COMMON/AXISY/C1BUL,CONST,CORAD,VT2,VTDR2,RTI2
      REAL*8 C1BUL, CONST, CORAD, VT2, VTDR2, RTI2
C
      REAL*8 R2, R, RREL, RREL2, POTBU, POTDI
      R2 = X * X + Y * Y
      R = SQRT (R2)
C     BULGE
      RREL = R / CORAD
      RREL2 = RREL * RREL
      POTBU = - C1BUL * LOG( RREL + SQRT( RREL2 + 1.)) / RREL
C     DISK
      POTDI = - 1.5 * VT2 * SQRT( 1.5 / ( 0.5 + R2 * RTI2))
C     TOTAL
      POT = POTBU + POTDI
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE POTBA( XINPUT, YINPUT, Z, EPT)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C POTENTIAL OF A FERRERS BAR WITH INDEX 1
C     SUBROUTINE PROL6 (ELIPM,A,B,C,X,Y,Z,AX,AY,AZ,EPT)
C      IMPLICIT real*8 (A-Z)
      IMPLICIT NONE
      real*8 XINPUT, YINPUT, Z, EPT
C
      COMMON/PROLA/ELIPM,A,B,C,AAIN,CCIN,AAI,CCO,ACONST
      REAL*8 ELIPM, A, B, C, ACONST, CCO, CCIN, AAI, AAIN
C
      REAL*8 X, Y, AA, CC, XX, YY, RR, GI, SIDE, BK, CK
      REAL*8 KAPPA, EM, EE, E, LNE, A3, A1, I, A33, A13, A11
      REAL*8 X4, Y4
C
C     CAREFUL : IN HYDRO PROGRAM THE BAR IS ALONG THE X AXIS
C     IN THE ORBIT CALCULATIONS ALONG THE Y
C
      X = YINPUT
      Y = XINPUT
      AA = AAIN
      CC = CCIN
      XX = X*X
      YY = Y*Y
      RR = XX+YY
      GI = 1.
      SIDE = XX/AA + (YY)/CC - 1.
      IF (SIDE .LE. 0.) GO TO 10
         BK = RR-AA-CC
         CK = AA*CC*SIDE
         KAPPA = (BK+SQRT(BK*BK+4.*CK))/2.
         AA = AA + KAPPA
         CC = CC + KAPPA
         GI = A*C*C/SQRT(AA)/CC
   10 EM = CC/AA
      EE = 1. - EM
      E = SQRT(EE)
      LNE = LOG((1.+E)/(1.-E))
C
C     CALCULATE INDEX SYMBOLS IN EASIEST WAY
C
      A3 = GI*(1.-EM*LNE/(E + E))/EE
      A1 = 2.*(GI - A3)
      I = AA*A1 + 2.*CC*A3
 
      A33 = GI*(4.*EE/EM - 6. + 3.*EM*LNE/E)/(8.*AA*EE*EE)
      A13 = 2.*(GI/CC - 2.*A33)
      A11 = 2.*(GI/AA - A13)/3.
 
      X4 = A11*XX + A13*(YY)
      Y4 = A13*XX + A33*(YY)
 
      EPT = I - 2.*(A1*XX + A3*(YY)) +(XX*X4+YY*Y4)
 
C     ACONST = ELIPM*4.027209375/(A*C*C)
      EPT = -EPT * ACONST
      RETURN
      END
C-------------------------------------------------------------------------------
      SUBROUTINE POTBA0( X, Y, Z, POT) 
C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------
C     POTENTIEL D'UN SPHEROIDE HOMOGENE PROLAT 
C     PARAMETRES
C     A : PETIT AXE 
C     C : GRAND AXE 
C     EP : EXCENTRICITE   EP = SQRT(C * C - A * A)
C     SHO = SINH(PSIO) .........
      IMPLICIT NONE
      REAL*8 X, Y, Z, POT
C
      COMMON/PROLA2/COEFF,COEFP,EP,EP2,SHO,CHO,THO,PSIO,CSHO,CTHO
      REAL*8 COEFF, COEFP, EP, EP2, SHO, CHO, THO, PSIO, CSHO, CTHO
C
      REAL*8 R2, Y2, XM, PSI, SH, CH, TH, AA, W10, W11, W20
      R2 = X * X
      Y2 = Y * Y
      XM = SQRT( R2 * CSHO + Y2 * CTHO)
      IF( XM .LE. 1.) THEN
         PSI = PSIO
         SH = SHO
         CH = CHO
         TH = THO
      ELSE
         AA = Y2 + R2 - EP2 
         IF( R2 * EP2 .LE. 1.D-30*AA) THEN
            SH = EP / SQRT( Y2 - EP2)
         ELSE
            SH = SQRT(( SQRT( AA * AA + 4. * R2 * EP2) - AA) / (R2+R2))
         END IF
         PSI = LOG( SH + SQRT( SH ** 2 + 1.)) 
         CH = DCOSH (PSI) 
         TH = SH / CH
      END IF
      W10 = PSI + PSI 
      W20 = .5 * W10 - CH * SH
      W11 = TH + TH - W10 
      POT= ((W20 * R2 + W11 * Y2) / EP2 + W10) * COEFP
      POT = - POT
      RETURN
      END 
