C            S.J.Aarseth's Standard N-Body Program
C
      REAL*8 x,x0,x0dot,t0,time
      DIMENSION x(3,50),x0(3,50),x0dot(3,50),t0(50),body(50),step(50),
     -     f(3,50),fdot(3,50),d1(3,50),d2(3,50),d3(3,50),t1(50),t2(50),
     -     t3(50),a(17),f1(3),f1dot(3),f2dot(3),f3dot(3)
      DATA time,tnext,nsteps /0.0,0.0,0/
      READ(5,*) n,eta,deltat,tcrit,eps2
      DO 1 i=1,n
 1    READ(5,*) body(i),(x0(k,i),k=1,3),(x0dot(k,i),k=1,3)
C           obtain total forces and first derivative for each body
      DO 20 i=1,n
      DO 2 k=1,3
        f(k,i)=0.0
        fdot(k,i)=0.0
        d2(k,i)=0.0
 2      d3(k,i)=0.0
      DO 10 j=1,n
      IF (j.EQ.i) GO TO 10
      DO 5 k=1,3
        a(k) = x0(k,j) - x0(k,i)
 5      a(k+3) = x0dot(k,j) - x0dot(k,i)
        a(7) = 1.0/(a(1)**2 + a(2)**2 + a(3)**2 + eps2)
        a(8) = body(j)*a(7)*sqrt(a(7))
        a(9) = 3.0*(a(1)*a(4)+a(2)*a(5)+a(3)*a(6))*a(7)
        DO 8 k=1,3
          f(k,i) = f(k,i) + a(k)*a(8)
 8        fdot(k,i) = fdot(k,i) + (a(k+3) - a(k)*a(9))*a(8)
 10   CONTINUE
 20   CONTINUE
C              form second and third derivative
      DO 40 i=1,n
        DO 30 j=1,n
          IF (i.EQ.j) GO TO 30
          DO 25 k=1,3
            a(k) = x0(k,j) - x0(k,i)
            a(k+3) = x0dot(k,j) - x0dot(k,i)
            a(k+6) = f(k,j) - f(k,i)
 25         a(k+9) = fdot(k,j) - fdot(k,i)
          a(13) = 1.0/(a(1)**2 + a(2)**2 + a(3)**2 + eps2)
          a(14) = body(j)*a(13)*sqrt (a(13))
          a(15) = (a(1)*a(4) + a(2)*a(5) + a(3)*a(6))*a(13)
          a(16) = (a(4)**2 + a(5)**2 + a(6)**2 + a(1)*a(7) + a(2)*a(8)
     -                                    + a(3)*a(9))*a(13) + a(15)**2
      a(17) = (9.0*(a(4)*a(7) + a(5)*a(8) + a(6)*a(9)) + 3.0*(a(1)*a(10)
     -+ a(2)*a(11) + a(3)*a(12)))*a(13) + a(15)*(9.*a(16) -12.*a(15)**2)
      DO 28 k=1,3
        f1dot(k) = a(k+3) - 3.0*a(15)*a(k)
        f2dot(k) = (a(k+6)-6*a(15)*f1dot(k)-3*a(16)*a(k))*a(14)
        f3dot(k) = (a(k+9)-9*a(16)*f1dot(k)-a(17)*a(k))*a(14)
        d2(k,i) = d2(k,i) + f2dot(k)
 28     d3(k,i) = d3(k,i) + f3dot(k) - 9.0*a(15)*f2dot(k)
   30 CONTINUE
   40 CONTINUE               
C              initialize integration steps and convert to force differences
      DO 50 i=1,n
        step(i) = sqrt(eta*sqrt ((f(1,i)**2 + f(2,i)**2 + f(3,i)**2)/
     -            (d2(1,i)**2 + d2(2,i)**2 + d2(3,i)**2)))
        t0(i) = time
        t1(i) = time - step(i)
        t2(i) = time - 2.0*step(i)
        t3(i) = time - 3.0*step(i)
        DO 45 k=1,3
          d1(k,i) = (d3(k,i)*step(i)/6-0.5*d2(k,i))*step(i)+fdot(k,i)
          d2(k,i) = 0.5*d2(k,i) - 0.5*d3(k,i)*step(i)
          d3(k,i) = d3(k,i)/6
          f(k,i) = 0.5*f(k,i)
 45       fdot(k,i) = fdot(k,i)/6
 50   CONTINUE
C              energy check and ouput
 100  e = 0.0
      DO 110 i=1,n
      dt = tnext - t0(i)
      DO 101 k=1,3
      f2dot(k) = d3(k,i)*((t0(i)-t1(i))+(t0(i)-t2(i))) + d2(k,i)
      x(k,i) = ((((0.05*d3(k,i)*dt + f2dot(k)/12.0)*dt + fdot(k,i))*dt
     -              + f(k,i))*dt + x0dot(k,i))*dt + x0(k,i)
 101        a(k) = (((0.25*d3(k,i)*dt + f2dot(k)/3.0)*dt + 
     -                3.0*fdot(k,i))*dt + 2.0*f(k,i))*dt + x0dot(k,i)
      WRITE(6,105) i,body(i),step(i),(x(k,i),k=1,3),(a(k),k=1,3)
 105  FORMAT(1h,i10,f10.2,f12.4,3x,3f10.2,3x,3f10.2)
 110  e = e + 0.5*body(i)*(a(1)**2 + a(2)**2 + a(3)**2)
      DO 130 i=1,n
        DO 120 j=1,n
          IF (j.EQ.i) GO TO 120
          e = e - 0.5*body(i)*body(j)/sqrt ((x(1,i) - x(1,j))**2 +
     -           (x(2,i) - x(2,j))**2 + (x(3,i) - x(3,j))**2 + eps2)
 120    CONTINUE
 130  CONTINUE
      WRITE(6,140) tnext,nsteps,e
  140 FORMAT(1h0,5x,'time =',f7.2,'  steps =',i6,' energy =',f10.4,/)
      IF (time.GT.tcrit) STOP
      tnext = tnext + deltat
C              find next body to be advanced and set new time   
 200  time = 1.0e+10
      DO 210 j=1,n
        IF (time.GT.t0(j)+step(j)) i=j
        IF (time.GT.t0(j)+step(j)) time = t0(j) + step(j)
  210 CONTINUE
C              predict all coordinates to first order in force derivative
      DO 220 j=1,n
        s = time - t0(j)
        x(1,j) = ((fdot(1,j)*s + f(1,j))*s + x0dot(1,j))*s + x0(1,j)
        x(2,j) = ((fdot(2,j)*s + f(2,j))*s + x0dot(2,j))*s + x0(2,j)
 220    x(3,j) = ((fdot(3,j)*s + f(3,j))*s + x0dot(3,j))*s + x0(3,j)        
C              include 2nd and 3rd order and obtain the velocity
      dt = time - t0(i)
      DO 230 k=1,3
        f2dot(k) = d3(k,i)*((t0(i)-t1(i)) + (t0(i)-t2(i))) + d2(k,i)
        x(k,i) = (0.05*d3(k,i)*dt + f2dot(k)/12.0)*dt**4 + x(k,i)
        x0dot(k,i) = (((0.25*d3(k,i)*dt + f2dot(k)/3.0)*dt +
     -                3.0*fdot(k,i))*dt + 2.0*f(k,i))*dt + x0dot(k,i)
 230    f1(k) = 0.0      
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
      DO 250 k=1,3
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
 250    x0dot(k,i) = (((0.2*a(k+9)*dt + 0.25*f3dot(k))*dt +
     -              f2dot(k)/3.0)*dt + 0.5*f1dot(k))*dt**2 + x0dot(k,i)
C           scale F and FDOT by factorials and set new integration step
      DO 260 k=1,3
        f(k,i) = 0.5*f1(k)
        fdot(k,i) = ((d3(k,i)*dt1 + d2(k,i))*dt + d1(k,i))/6.0
 260    f2dot(k)= 2.0*(d3(k,i)*(dt+dt1) + d2(k,i))
      step(i) = sqrt(eta*sqrt((f1(1)**2 + f1(2)**2 + f1(3)**2)/
     -         (f2dot(1)**2 + f2dot(2)**2 + f2dot(3)**2)))
      nsteps = nsteps + 1
      IF (time - tnext) 200,100,100
      END
