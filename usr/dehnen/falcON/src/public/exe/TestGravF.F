C--------------------------------------------------------------------------+
C                                                                          |
C TestGravC.cc                                                             |
C                                                                          |
C Copyright (C) 2000-2005 Walter Dehnen                                    |
C                                                                          |
C This program is free software; you can redistribute it and/or modify     |
C it under the terms of the GNU General Public License as published by     |
C the Free Software Foundation; either version 2 of the License, or (at    |
C your option) any later version.                                          |
C                                                                          |
C This program is distributed in the hope that it will be useful, but      |
C WITHOUT ANY WARRANTY; without even the implied warranty of               |
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU        |
C General Public License for more details.                                 |
C                                                                          |
C You should have received a copy of the GNU General Public License        |
C along with this program; if not, write to the Free Software              |
C Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                |
C                                                                          |
C--------------------------------------------------------------------------+
C                                                                          |
C     Here, we use C-style macros to obtain correct input data type.       |
C     For this to work, we need first to process this file with the        |
C     C preprocessor cpp, see makefile make.public                         |
C                                                                          |
C--------------------------------------------------------------------------+
#ifdef falcON_DOUBLE
#  define INPUT_TYPE DOUBLE PRECISION
#else
#  define INPUT_TYPE REAL
#endif
C---------------------------------------------------------------------------
      PROGRAM TESTGRAVF

      IMPLICIT NONE

      INCLUDE 'forces.f'

      INTEGER NMAX
      INPUT_TYPE RMAX,PI,TPI,GRAV
      PARAMETER(NMAX=1000000)
      PARAMETER(RMAX=1.E3, PI=3.14159265358979323846264338328, TPI=2*PI)
      PARAMETER(GRAV=1)

      INTEGER    MODEL,K,N,NCRIT,I,J,S,CPU0,CPU1
      INPUT_TYPE EPS,TH,GAM,TIME,DAMX
      INPUT_TYPE G1,G2,IG2,G3,IG3,TM,RAD,RCY
      INPUT_TYPE PHI,FR,P,CTH,MI,ASE,DACC,SF(3)
      INTEGER    FL(NMAX)
      INPUT_TYPE M(NMAX),PH(NMAX),RH(NMAX)
      INPUT_TYPE X(3,NMAX),A(3,NMAX),AT(3,NMAX)
      INPUT_TYPE MSFORCE

      CPU0=CLOCK()

      WRITE(*,'("model? (0/1/2/3/4 = hom disk /Plummer /gamma-model")')
      WRITE(*,'("                    /Kuzmin disk /hom disk)")')
      READ (*,*) MODEL
      IF(MODEL.EQ.2) THEN
         WRITE(*,'("gamma for gamma-model ? (1->Hernquist model)")')
         READ (*,*) GAM
      ELSE
         GAM = 0
      ENDIF
      WRITE(*,'("N ? (max "I7")")') NMAX
      READ(*,*) N
      WRITE(*,'("seed ?")')
      READ(*,*) S
      WRITE(*,'("epsilon ?")')
      READ(*,*) EPS
      WRITE(*,'("theta? (if <0, use theta=const)")')
      READ(*,*) TH
      WRITE(*,'("kernel ? (n) = P_n)")')
      READ(*,*) K
      WRITE(*,'("Ncrit ? ")')
      READ(*,*) NCRIT
      
      CALL SRAND48(S)

      G1  = 1.-GAM
      G2  = 2.-GAM
      IG2 = 1./G2
      G3  = 3.-GAM
      IG3 = 1./G3
      MI  = 1./REAL(N)

      DO I=1,N

 10      IF(MODEL.EQ.0) THEN

            RAD = DRAND48()**(1./3.)
            P   = 0.5*(3*RAD*RAD-1)
            FR  =-1.
         ELSE IF(MODEL.EQ.1) THEN
            TM  = DRAND48()**(2./3.)
            TM  = TM/(1.-TM)
            RAD = SQRT(TM)
            P   =-1./SQRT(TM+1)
            FR  = P*P*P
         ELSE IF(MODEL.EQ.2) THEN
            TM  = DRAND48()**IG3
            RAD = TM/(1.-TM)
            IF(GAM.EQ.2.) THEN
               P= LOG(TM)
            ELSE
               P= -IG2*(1-TM**G2)   
            ENDIF
            FR  =-TM**G1 * (1-TM)**2 / RAD
         ELSE IF(MODEL.EQ.3) THEN
            TM  = DRAND48()
            RAD = SQRT((2-TM)*TM)/(1.-TM)
            P   =-1./SQRT(RAD**2+1)
            FR  = P*P*P
         ELSE
            TM  = DRAND48()
            RAD = SQRT(TM)
            P   = 2*(RAD*RAD/3.-2)
            FR  =-4*RAD/3.
         ENDIF
         IF(RAD.GT.RMAX) GOTO 10
         IF(MODEL.GE.3) THEN
            CTH    = 0.
         ELSE
            CTH    = 2.*DRAND48()-1
         ENDIF
         RCY    = RAD*SQRT(1-CTH*CTH)
         PHI    = TPI*DRAND48()
         FL(I)  = 1
         X(1,I) = RCY * SIN(PHI)
         X(2,I) = RCY * COS(PHI)
         X(3,I) = RAD * CTH
         M(I)   = MI
         AT(1,I)= FR * X(1,I)
         AT(2,I)= FR * X(2,I)
         AT(3,I)= FR * X(3,I)
C TEST
C         WRITE(*,'("X="3E25.13)') X(I),Y(I),Z(I)
C TSET
      ENDDO
      CPU1 = CLOCK()
      TIME = (CPU1-CPU0) / 1.E6
      WRITE(*,'(" time needed for set up of X_i:  "F13.6)') TIME
      CPU0 = CPU1
      CALL FALCON_INITIALIZE(FL,M,X,
#ifdef falcON_INDI
     C%VAL(0),
#endif
     CA,PH,RH, N,0,EPS,TH,K,GRAV)

      DO I=1,2
         CALL FALCON_GROW(NCRIT)
         CPU1 = CLOCK()
         TIME = (CPU1-CPU0) / CLOCKS_PER_SECOND()
         WRITE(*,'(" time needed by growing a tree:  "F13.6)') TIME
         CPU0 = CPU1
      ENDDO

      CALL FALCON_APPROX_GRAV()
      CPU1 = CLOCK()
      TIME = (CPU1-CPU0) / CLOCKS_PER_SECOND()
      WRITE(*,'(" time needed by force algorithm: "F13.6)') TIME
      CPU0 = CPU1

      CALL FALCON_STATS()

      ASE   = 0
      DAMX  = 0
      SF(1) = 0
      SF(2) = 0
      SF(3) = 0
      DO I=1,N
         DACC = 0
         DO J=1,3
            SF(J) = SF(J) + M(I)*A(J,I) 
            DACC  = DACC + (A(J,I) - AT(J,I))**2
         ENDDO
         ASE   = ASE + M(I)*DACC
         IF(DACC.GT.DAMX) THEN
            DAMX = DACC
         ENDIF
C           WRITE(*,'(" "6F13.8)') X(1,I),X(2,I),X(3,I),A(1,I),A(2,I),A(3,I)
      ENDDO
      IF(MODEL.EQ.1) THEN
         ASE = ASE * 13.125
      ELSE IF(MODEL.EQ.3) THEN
         ASE = ASE * 17.5
      ELSE IF(MODEL.EQ.4) THEN
         ASE = ASE * 1.125
      ELSE IF(MODEL.EQ.0) THEN
         ASE = ASE / 0.6
      ELSE IF(MODEL.EQ.2 .AND. GAM.LT.1.66) THEN
         ASE = ASE / MSFORCE(GAM)
      ENDIF
      IF(MODEL.EQ.2 .AND. GAM.GE.1.66) THEN
         WRITE(*,'(" ASE(F)        = "F13.10)') ASE
      ELSE
         WRITE(*,'(" ASE(F)/<F^2>  = "F13.10)') ASE
      ENDIF
      WRITE(*,'(" max (dF)^2    = "F13.10)') DAMX
      WRITE(*,'(" Sum m_i acc_i = "3E15.6)') SF(1),SF(2),SF(3)

      END

      INPUT_TYPE FUNCTION MSFORCE(GAM)
      INPUT_TYPE GAM
      INPUT_TYPE GAMMA,X,F
      GAMMA = GAM
      X = 5-3*GAMMA
      IF(X .LT. 0) THEN
         MSFORCE = 0
      ELSE
         F =     1/X
         X = X + 1
         F = F - 4/X
         X = X + 1
         F = F + 6/X
         X = X + 1
         F = F - 4/X
         X = X + 1
         F = F + 1/X
         MSFORCE = (3-GAMMA) * F
      ENDIF
      RETURN
      END

