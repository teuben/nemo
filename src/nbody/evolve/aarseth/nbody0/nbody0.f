C***********************************************************************
       SUBROUTINE nbody0
C
C	S.J. Aarseth's Standard N-Body Program (nbody0)
C
C	Adapted from the source code in Binney & Tremaine's (1987)
C       book 'Galactic Dynamics',
C         by Peter Teuben - June '89:
C       - all variables to be declared and a choice of 'real' or 
C         'double precision'
C	- all I/O is done through subroutines in order to allow easy
C	  interface with development different systems, e.g. NEMO89.
C	- PARAMETERS for NMAX and NDIM via an include file (see Makefile)
C	  Note that NDIM should NEVER be reset from 3, the code does
C	  not handle 2D stuff (yet).
C
C          jun-89   Created
C	23-jan-90   Back to double precision  - PJT
C	 5-apr-90   started to declare all variables - 	PJT
C       14-nov-91   finished(!) declaring all variables; 
C                   tested using IMPLICIT NONE and compile with -u
C       15-nov-91   Connection Machine testing in nbody0.fcm file
C       11-feb-98   Fixed array index for f2dot(k)  a(16)->a(15)
C       22-jan-00   Fixed reset time > tnext after datadumps
C       21-feb-04   stop if d2test is 0
C***********************************************************************
      INCLUDE 'fdefs.inc'
      INCLUDE 'nmax.inc'
CSED	The next statement can be modified with SED to toggle type
      DOUBLE PRECISION
     -    time,tnext,dt,s,deltat,tcrit,e,eta,eps2,
     -    t1pr,t2pr,t3pr,dt1,dt2,dt3,d2test,
     -    x(NDIM,NMAX),x0(NDIM,NMAX),x0dot(NDIM,NMAX),t0(NMAX),
     -    body(NMAX),step(NMAX),
     -    f(NDIM,NMAX),fdot(NDIM,NMAX),
     -    d1(NDIM,NMAX),d2(NDIM,NMAX),d3(NDIM,NMAX),
     -    t1(NMAX),t2(NMAX),t3(NMAX),a(17),
     -    f1(NDIM),f1dot(NDIM),f2dot(NDIM),f3dot(NDIM)
      INTEGER   nsteps, n, i, j, k, makesure,reset,use3dot

      DATA  time,tnext,nsteps /0.0,0.0,0/
C-----------------------------------------------------------------------
      CALL INPARS(NMAX,n,eta,deltat,tcrit,eps2,reset,use3dot)
      CALL INBODS (n, body, x0, x0dot)
   
C           obtain total forces and first derivative for each body
      DO 20 i=1,n
         DO 2 k=1,NDIM
            f(k,i)=0.0
            fdot(k,i)=0.0
            d2(k,i)=0.0
            d3(k,i)=0.0
    2    CONTINUE
         DO 10 j=1,n
            IF (j.EQ.i) GO TO 10
            DO 5 k=1,NDIM
               a(k) = x0(k,j) - x0(k,i)
               a(k+3) = x0dot(k,j) - x0dot(k,i)
    5       CONTINUE
            a(7) = 1.0/(a(1)**2 + a(2)**2 + a(3)**2 + eps2)
            a(8) = body(j)*a(7)*sqrt(a(7))
            a(9) = 3.0*(a(1)*a(4)+a(2)*a(5)+a(3)*a(6))*a(7)
            DO 8 k=1,NDIM
               f(k,i) = f(k,i) + a(k)*a(8)
               fdot(k,i) = fdot(k,i) + (a(k+3) - a(k)*a(9))*a(8)
    8       CONTINUE
   10    CONTINUE
   20 CONTINUE
C              form second and third derivative
      DO 40 i=1,n
         DO 30 j=1,n
            IF (i.EQ.j) GO TO 30
            DO 25 k=1,NDIM
               a(k) = x0(k,j) - x0(k,i)
               a(k+3) = x0dot(k,j) - x0dot(k,i)
               a(k+6) = f(k,j) - f(k,i)
               a(k+9) = fdot(k,j) - fdot(k,i)
   25       CONTINUE
            a(13) = 1.0/(a(1)**2 + a(2)**2 + a(3)**2 + eps2)
            a(14) = body(j)*a(13)*sqrt (a(13))
            a(15) = (a(1)*a(4) + a(2)*a(5) + a(3)*a(6))*a(13)
            a(16) = (a(4)**2 + a(5)**2 + a(6)**2 +
     -             a(1)*a(7) + a(2)*a(8) + a(3)*a(9))*a(13) + a(15)**2
            a(17) = (9.0*(a(4)*a(7) + a(5)*a(8) + a(6)*a(9)) + 
     -             3.0*(a(1)*a(10) + a(2)*a(11) + a(3)*a(12)))*a(13) +
     -              a(15)*(9.0*a(16) - 12.0*a(15)**2)
            DO 28 k=1,NDIM
               f1dot(k) = a(k+3) - 3.0*a(15)*a(k)
               f2dot(k) = (a(k+6)-6*a(15)*f1dot(k)-3*a(16)*a(k))*a(14)
               f3dot(k) = (a(k+9)-9*a(16)*f1dot(k)-a(17)*a(k))*a(14)
               d2(k,i) = d2(k,i) + f2dot(k)
               d3(k,i) = d3(k,i) + f3dot(k) - 9.0*a(15)*f2dot(k)
   28       CONTINUE
   30    CONTINUE
   40 CONTINUE               
     
C              initialize integration steps and convert to force differences

      DO 50 i=1,n
         d2test = d2(1,i)**2 + d2(2,i)**2 + d2(3,i)**2
         if (d2test.eq.0) RETURN
         step(i) = sqrt(eta*sqrt ((f(1,i)**2 + f(2,i)**2 + f(3,i)**2)/
     -              (d2(1,i)**2 + d2(2,i)**2 + d2(3,i)**2)))
         t0(i) = time
         t1(i) = time - step(i)
         t2(i) = time - 2.0*step(i)
         t3(i) = time - 3.0*step(i)
         DO 45 k=1,NDIM
            d1(k,i) = (d3(k,i)*step(i)/6-0.5*d2(k,i))*step(i)+fdot(k,i)
            d2(k,i) = 0.5*d2(k,i) - 0.5*d3(k,i)*step(i)
            d3(k,i) = d3(k,i)/6
            f(k,i) = 0.5*f(k,i)
            fdot(k,i) = fdot(k,i)/6
   45    CONTINUE
   50 CONTINUE
C-----------------------------------------------------------------------
C==============> entry point when major output + diagnostics to be done
C   
C              energy check and ouput
  100 CONTINUE
      e = 0.0
      DO 110 i=1,n
         dt = tnext - t0(i)
         DO 101 k=1,NDIM
            f2dot(k) = d3(k,i)*((t0(i)-t1(i))+(t0(i)-t2(i))) + d2(k,i)
            x(k,i) = ((((0.05*d3(k,i)*dt + f2dot(k)/12.0)*dt +
     -                  fdot(k,i))*dt + f(k,i))*dt + x0dot(k,i))*dt
     -               + x0(k,i)
            a(k) = (((0.25*d3(k,i)*dt + f2dot(k)/3.0)*dt + 
     -                3.0*fdot(k,i))*dt + 2.0*f(k,i))*dt + x0dot(k,i)
  101    CONTINUE
         CALL OUTBODS (body(i),x(1,i),a,step(i),i)
      e = e + 0.5*body(i)*(a(1)**2 + a(2)**2 + a(3)**2)
  110 CONTINUE
      DO 130 i=1,n
         DO 120 j=1,n
            IF (j.EQ.i) GO TO 120
            e = e - 0.5*body(i)*body(j)/sqrt ((x(1,i) - x(1,j))**2 +
     -           (x(2,i) - x(2,j))**2 + (x(3,i) - x(3,j))**2 + eps2)
  120    CONTINUE
  130 CONTINUE
      CALL OUTENE(tnext,nsteps,e)
      IF (time.GT.tcrit) RETURN
      tnext = tnext + deltat
      makesure = 1

C-----------------------------------------------------------------------
C   ===============>  Normal entry point for next timstep
C
C    
C              find next body to be advanced and set new time
   
  200 CONTINUE
      time = 1.0e+10
      DO 210 j=1,n
         IF (time.GT.t0(j)+step(j)) i=j
         IF (time.GT.t0(j)+step(j)) time = t0(j) + step(j)
  210 CONTINUE
c	write(*,*) 'time=',time,' i=',i
      IF(makesure.GT.0) THEN
c         write (*,*) 'TIME?: ',time,tnext
         IF(reset.GT.0 .AND. time.GT.tnext) THEN
c           write (*,*) 'reset done!'
	   time = tnext
         ENDIF
         makesure = 0
      ENDIF
 
C              predict all coordinates to first order in force derivative

      DO 220 j=1,n
         s = time - t0(j)
         x(1,j) = ((fdot(1,j)*s + f(1,j))*s + x0dot(1,j))*s + x0(1,j)
         x(2,j) = ((fdot(2,j)*s + f(2,j))*s + x0dot(2,j))*s + x0(2,j)
         x(3,j) = ((fdot(3,j)*s + f(3,j))*s + x0dot(3,j))*s + x0(3,j)
  220 CONTINUE
           
C              include 2nd and 3rd order and obtain the velocity

      dt = time - t0(i)
      DO 230 k=1,NDIM
         f2dot(k) = d3(k,i)*((t0(i)-t1(i)) + (t0(i)-t2(i))) + d2(k,i)
         x(k,i) = (0.05*d3(k,i)*dt + f2dot(k)/12.0)*dt**4 + x(k,i)
         x0dot(k,i) = (((0.25*d3(k,i)*dt + f2dot(k)/3.0)*dt +
     -                3.0*fdot(k,i))*dt + 2.0*f(k,i))*dt + x0dot(k,i)
         f1(k) = 0.0
  230 CONTINUE         

C              obtain current force on i-th body

      DO 240 j=1,n
         IF (j.EQ.i) GO TO 240
         a(1) = x(1,j) - x(1,i)
         a(2) = x(2,j) - x(2,i)
         a(3) = x(3,j) - x(3,i)
         a(4) = 1.0/(a(1)**2 + a(2)**2 + a(3)**2 + eps2)
         a(5) = body(j)*a(4)*sqrt(a(4))
         f1(1) = f1(1) + a(1)*a(5)
         f1(2) = f1(2) + a(2)*a(5)
         f1(3) = f1(3) + a(3)*a(5)
  240 CONTINUE

C              set time intervals for new difference and update the times

      dt1 = time - t1(i)
      dt2 = time - t2(i)
      dt3 = time - t3(i)
      t1pr = t0(i) - t1(i)
      t2pr = t0(i) - t2(i)
      t3pr = t0(i) - t3(i)
      t3(i) = t2(i)
      t2(i) = t1(i)
      t1(i) = t0(i)
      t0(i) = time

C              form new differences and include 4th order semi-iterative

      DO 250 k=1,NDIM
         a(k) = (f1(k) - 2.0*f(k,i))/dt
         a(k+3) = (a(k) - d1(k,i))/dt1
	 a(k+6) = (a(k+3) - d2(k,i))/dt2
	 a(k+9) = (a(k+6) - d3(k,i))/dt3
	 d1(k,i) = a(k)
	 d2(k,i) = a(k+3)
	 d3(k,i) = a(k+6)
	 f1dot(k) = t1pr*t2pr*t3pr*a(k+9)
	 f2dot(k) = (t1pr*t2pr + t3pr*(t1pr+t2pr))*a(k+9)
	 f3dot(k) = (t1pr + t2pr + t3pr)*a(k+9)
	 x0(k,i) = (((a(k+9)*dt/30.0 + 0.05*f3dot(k))*dt +
     -              f2dot(k)/12.0)*dt + f1dot(k)/6.0)*dt**3 + x(k,i)
         x0dot(k,i) = (((0.2*a(k+9)*dt + 0.25*f3dot(k))*dt +
     -              f2dot(k)/3.0)*dt + 0.5*f1dot(k))*dt**2 + x0dot(k,i)
  250 CONTINUE
    
C           scale F and FDOT by factorials and set new integration step

      DO 260 k=1,NDIM
         f(k,i) = 0.5*f1(k)
         fdot(k,i) = ((d3(k,i)*dt1 + d2(k,i))*dt + d1(k,i))/6.0
         f2dot(k)= 2.0*(d3(k,i)*(dt+dt1) + d2(k,i))
  260 CONTINUE
      step(i) = sqrt(eta*sqrt((f1(1)**2 + f1(2)**2 + f1(3)**2)/
     -         (f2dot(1)**2 + f2dot(2)**2 + f2dot(3)**2)))
      nsteps = nsteps + 1
      IF (time - tnext) 200,100,100
      END
