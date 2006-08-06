C--------------------------------------------------------------------------+
C     Here, we use C-style macros to obtain correct input data type.       |
C     For this to work, we need first to process this file with the        |
C     C preprocessor cpp, see makefile make.public                         |
C--------------------------------------------------------------------------+
#if defined(falcON_SINGLE_DOUBLE) || defined(falcON_DOUBLE_DOUBLE)
#define INPUT_TYPE DOUBLE PRECISION
#else
#define INPUT_TYPE REAL
#endif
C---------------------------------------------------------------------------
      SUBROUTINE DMINTMIN(R,V,TAU,DMIN,TMIN)
      INPUT_TYPE R(3),V(3)
      INPUT_TYPE TAU,TMIN,DMIN,RV,VQ,TMP
      INTEGER    I
      IF(TAU .LE. 0.) THEN
         TMIN = 0.
         DMIN = 0.
         DO I=1,3
            DMIN = DMIN + R(I)*R(I)
         ENDDO
         DMIN = SQRT(DMIN)
      ELSE
         RV = 0
         VQ = 0
         DO I=1,3
            RV = RV + R(I)*V(I)
            VQ = VQ + V(I)*V(I)
         ENDDO
         TMIN = -RV/VQ
         IF(TMIN .LT. 0. ) TMIN = 0.
         IF(TMIN .GT. TAU) TMIN = TAU
         DMIN = 0.
         DO I=1,3
            TMP  = R(I) + TMIN * V(I)
            DMIN = DMIN + TMP*TMP
         ENDDO
         DMIN = SQRT(DMIN)
      ENDIF
      RETURN
      END
C---------------------------------------------------------------------------
      PROGRAM TESTPAIRF

      IMPLICIT NONE

      INCLUDE 'forces.f'

      INTEGER MAX,NMAX
      INPUT_TYPE RMAX,PI,TPI
      PARAMETER(NMAX=100000)
      PARAMETER(RMAX=1.E3, PI=3.14159265358979323846264338328, TPI=2*PI)

      INTEGER    N,NCRIT,NS,NA,IN,I,S,CPU0,CPU1,I1(NMAX),I2(NMAX)
      INPUT_TYPE TAU,SIZE,TIME,DIST,TMP,DX(3),DV(3)
      INPUT_TYPE DMIN,TMIN,EPS,TH
      INTEGER    FL(NMAX),ZERO(3),D,K
      INPUT_TYPE M(NMAX),X(NMAX),Y(NMAX),Z(NMAX)
      INPUT_TYPE VX(NMAX),VY(NMAX),VZ(NMAX),SZ(NMAX)
      INPUT_TYPE MR,RQ,RAD,POT,VEQ,P,VQ,V,F0,FS,CTH,PHI,RCY,MI

      CPU0=CLOCK()

      DO i=1,3
         ZERO(i)=0
      END DO
      D   = 100
      K   = 1
      EPS = 0.1
      TH  = 0.4

      WRITE(*,'(" Testing the pair finder using a Plummer model")')
      WRITE(*,'(" sticky particle for tau>=0, SPH for tau<0")')
      WRITE(*,'("N(<="I6")   = ")') NMAX
      READ(*,*) N
      WRITE(*,'("SEED            = ")')
      READ(*,*) S
      WRITE(*,'("NS <= N         = ")')
      READ(*,*) NS
      WRITE(*,'("body size       = ")')
      READ(*,*) SIZE
      WRITE(*,'("TAU             = ")')
      READ(*,*) TAU
      WRITE(*,'("NCRIT           = ")')
      READ(*,*) NCRIT
      WRITE(*,'(" |r| < max(h_i,h_j) (1) OR")')
      WRITE(*,'(" |r| < h_i+h_j      (0)  ? :")')
      READ(*,*) MAX
      
      CALL SRAND48(S)
      MI  =1./N

      DO I=1,N
         MR  = DRAND48()**(2./3.)
         RQ  = MR/(1-MR)
         POT =-1/SQRT(1+RQ)
         RAD = SQRT(RQ)
         VEQ =-2*POT
 10         P  = VEQ*DRAND48()**(1/4.5)
            VQ = VEQ-P
            V  = SQRT(VQ)
            F0 = V*P**3.5
            FS = VEQ*F0
            F0 = F0*V
         IF(FS*DRAND48() .GT. F0) GOTO 10
         CTH   = 2*DRAND48()-1
         PHI   = TPI*DRAND48()
         RCY   = RAD * SQRT(1-CTH**2)
         X(I)  = RCY * SIN(PHI)
         Y(I)  = RCY * COS(PHI)
         Z(I)  = RAD * CTH
         CTH   = 2*DRAND48()-1
         PHI   = TPI*DRAND48()
         RCY   = V * SQRT(1-CTH**2)
         VX(I) = RCY * SIN(PHI)
         VY(I) = RCY * COS(PHI)
         VZ(I) = V   * CTH
         M(I)  = MI
         SZ(I) = SIZE
         IF(I .LE. NS) THEN 
            FL(I) = 13
         ELSE
            FL(I) = 1
         ENDIF
      ENDDO
      CPU1 = CLOCK()
      TIME = (CPU1-CPU0) / 1.E6
      CPU0 = CPU1
      WRITE(*,'(" time needed for set up of (x,v): "F13.6)') TIME
      CALL FALCON_INITIALIZE(FL,M,X,Y,Z,
#ifdef falcON_INDI
     C%VAL(0),
#endif
     C%VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0), N,0.1,0.55,1)

      CALL FALCON_GROW(NCRIT)

      CPU1  = CLOCK()
      TIME = (CPU1-CPU0) / CLOCKS_PER_SECOND()
      CPU0  = CPU1
      WRITE(*,'(" time needed by tree growth:      "F13.6)') TIME

      IF(TAU .LT. 0) THEN
         CALL FALCON_SPH(I1,I2,NMAX,NA,MAX,SZ)
      ELSE
         CALL FALCON_STICKY(TAU,I1,I2,NMAX,NA,VX,VY,VZ,SZ)
      ENDIF

      CPU1  = CLOCK()
      TIME = (CPU1-CPU0) / CLOCKS_PER_SECOND()
      CPU0  = CPU1
      WRITE(*,'(" time needed by pair search:      "F13.6)') TIME

      IF(TAU .LT. 0) THEN
         OPEN(5,FILE='pairs.for')
         IN = NA
         IF(NA .GT. NMAX) IN = NMAX
         DO I = 1,IN
            DIST = 0.
            TMP  = X(I1(I)) - X(I2(I))
            DIST = DIST + TMP*TMP
            TMP  = Y(I1(I)) - Y(I2(I))
            DIST = DIST + TMP*TMP
            TMP  = Z(I1(I)) - Z(I2(I))
            DIST = DIST + TMP*TMP
            DIST = SQRT(DIST)
            WRITE(5,'(I8" "I8"  "F13.6)') I1(I), I2(I), DIST
         ENDDO
         CLOSE(5)
      ELSE
         OPEN(5,FILE='pairs.for')
         IN = NA
         IF(NA .GT. NMAX) IN = NMAX
         DO I = 1,IN
            DX(1) =  X(I1(I)) -  X(I2(I))
            DX(2) =  Y(I1(I)) -  Y(I2(I))
            DX(3) =  Z(I1(I)) -  Z(I2(I))
            DV(1) = VX(I1(I)) - VX(I2(I))
            DV(2) = VY(I1(I)) - VY(I2(I))
            DV(3) = VZ(I1(I)) - VZ(I2(I))
            CALL DMINTMIN(DX,DV,TAU,DMIN,TMIN)
            WRITE(5,'(I8" "I8"  "F13.6" "F13.6)') I1(I),I2(I),TMIN,DMIN
         ENDDO
         CLOSE(5)
      ENDIF
      IF(NA .LT. NMAX) THEN
         WRITE(*,'(I8" pairs found and written to file pairs.for")')
     C         NA
      ELSE 
         WRITE(*,'(I8" pairs found; "I8" written to file pairs.for")')
     C         NA,NMAX
      ENDIF

      END
