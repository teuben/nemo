CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  CGS - Collisionless Galactic Simulator - Source Files  C  
C                                                         C
C  Created by M. Trenti & T.van Albada in Fortran77 2003  C
C                                                         C
C  Version 1.0 Alpha                                      C
C                                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C_________________________________________________________
C SUBROUTINE FROM THE OLD VERSION OF THE CODE
      SUBROUTINE PLUMMER (SEED,RRC,RRT)
C
C       **DETERMINES COORDINATES AND VELOCITIES OF N
C         PARTICLES IN A POLYTROPE N=5 (PLUMMER MODEL)
C         AND STORES THE RESULTS IN THE COMMON BLOCK FPCORD
C
C         PARAMETERS :
C     NP   : NUMBER OF PARTICLES         !via common.blk
C     MM          : MASS OF POLYTROPE    !via common.blk
C         RRC         : CORE RADIUS
C         RRT         : TRUNCATION RADIUS
C         X,Y,Z      : PARTICLE COORDINATES
C         VX,VY,VZ   : PARTICLE VELOCITIES
C
C         UNITS :
C         MASS       : 10**6 MSUN
C         LENGTH     : 10 PC
C         TIME       : 10**6 YEAR
C         G          : 4.497
C         VELOCITY   : 9.78 KM/S
C       **
C

C_________________________________________________________
      include 'common.blk'

C_________________________________________________________
      include 'com_part.blk'

      real mm,rrc,rrt  
        

      REAL PSI
C     
      REAL  U,U1,U2,U3,UU,FU,R,A,SS
C
      INTEGER SEED
C
C
       
C
      REAL GM,RC,RC2,RT
C 

c random number common 
        integer idum
        common /ised/ idum

        real ran2
        external ran2

C___________________________________________________
C MM equivalent to Mtot
        MM=Mtot

c units adjustment
        GM=MM*GRAVC             !G*Mass
        RC=RRC                  !core radius
        RC2=RRC*RRC             !square rc
        RT=RRT                  !truncation radius


c random number init
        idum=seed


c     open storage file
c       open(10,file='initconPL.dat')
      

C
C       **DETERMINE COORDINATES**
C
             DO 100 L= 1,NP

   10           U= RAN2()
                R= RC*U**(1.0/3.0)/SQRT(1.0-U**(2.0/3.0))
                IF (R.GT.RT) GO TO 10
C
C       **RADIUS R ACCEPTED, DETERMINE X,Y,Z**
C
   20           U1= 2.0*RAN2()-1.0
                U2= 2.0*RAN2()-1.0
                U3= 2.0*RAN2()-1.0
                SS= SQRT(U1*U1+U2*U2+U3*U3)
                IF ((SS.GT.1.0).OR.(SS.LT.0.01))GO TO 20
C
C       **DIRECTION COSINES ACCEPTED**
C
                X(L)= R*U1/SS
                Y(L)= R*U2/SS
                Z(L)= R*U3/SS
                S(L)=R          !to be checked 
                TH(L)=ACOS(Z(L)/S(L))
                PH(L)=ATAN2(Y(L),X(L))
                IF (PH(L).LT.0.0) PH(L)=PH(L)+2.0*PI
C
C       **DETERMINE VELOCITIES**
C
                PSI= -GM/SQRT(RC2+R*R)
                A= -2.0*PSI


   30           U= RAN2()
                FU= 4.5*(9.0/7.0*(1.0-U*U))**3.5*U*U
                UU= RAN2()
                IF (UU.GT.FU) GO TO 30
C
C       **VELOCITY ACCEPTED**
C
                V= U*SQRT(A)
C
C       **DETERMINE VX,VY,VZ**
C
   40           U1= 2.0*RAN2()-1
                U2= 2.0*RAN2()-1
                U3= 2.0*RAN2()-1
                SS= SQRT(U1*U1+U2*U2+U3*U3)
                IF ((SS.GT.1.0).OR.(SS.LT.0.01)) GO TO 40
C
C       **DIRECTION COSINES ACCEPTED**
C
                VX(L)= V*U1/SS
                VY(L)= V*U2/SS
                VZ(L)= V*U3/SS
 
C            DUMMY-VALUES IP,JP,KP AND FORCES                           
      IP(L)=1                                                          
      FX0(L)=0.0                                                        
      FY0(L)=0.0                                                        
      FZ0(L)=0.0                  


C______________________________________________________________
C SAVES PARTICLES POSITION & VELOCITY
c            WRITE(10,222) X(L),Y(L),Z(L),VX(L),VY(L),VZ(L)
c222         FORMAT(6(F15.9))
              

  100           CONTINUE


c            close(10)
C

        END
C
C*********************************************************


c
c 
c uniform random between 0 and 1
c from Press et al. Numerical recipes in Fortran.
c
c modified for init (common block) 
c
c NOTE initialize IDUM < 0
c 
      FUNCTION ran1()
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      DOUBLE PRECISION ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      
      common /ised/ idum
      
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END



ccccccccccccccccccccccccccccccccc

      FUNCTION ran2()

      integer idum
      common /ised/ idum

      INTEGER IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
