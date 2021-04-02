      program anisot
c
c        This program is supposed to calculate the distribution function
c     for a given density distribution.  The density distribution will be
c     given numerically, in the form of a vector.  The potential will also
c     have to be computed numerically.  The models are supposed to be regular
c     in the center, and finite.
c
c        Once density and potential are known, the distribution
c     function is be computed with the aid of the formula
c
c           f(E) = {\sqrt{2} \over 4\pi} {d \over dE} \int_E^0
c        {d \rho1 \over d U} { d U \over \sqrt{U-E} } \> ,
c
c     \rho1  being the `corrected' (in Merritt's sense) density.
c
c        It may be possible to integrate by parts.  In this case, one
c     can avoid either numerical differentiation or integration over a
c     singularity, possibly both.  We will indeed use the alternative method
c     of integrating once by parts to calculate  DISTR2, as a check.
c     The first tests seem to indicate that DISTR2 is LESS accurate than DISTR.
c
c        The vectors RADIUS, POTEN and DENS will contain the potential and
c     density as a function of radius.  They will be filled in a separate
c     subroutine, which is also expected to fill in  EMINT,  the mass interior
c     to a given radius, and  THDIS,  the theoretical distribution
c     function (when known at least for the isotropic problem).
c
      implicit double precision (a-h, o-z)
      parameter (npoint = 10001)
      dimension poten(npoint), radius(npoint), distr(npoint)
      dimension dens(npoint), dens1(npoint), emint(npoint)
      dimension distr2(npoint), value(npoint), thdis(npoint)
      dimension aux1(npoint), aux2(npoint), aux3(npoint)
      dimension cdens(npoint, 3), cdens1(npoint, 3), cpoten(npoint, 3)
      dimension caux1(npoint, 3), caux2(npoint, 3), caux3(npoint, 3)
      dimension cvalue(npoint, 3)
      dimension qplus(npoint), index(npoint)
      parameter (ngauto = 200)
      dimension abscis(ngauto), weight(ngauto)
      data idens, idens1, ipoten, ivalue, iaux1,  iaux2,  iaux3
     >  / npoint, npoint, npoint, npoint, npoint, npoint, npoint /
      data pi, g / 3.1415926535897932d0, 1.d0 /
c
      write (6, 6400)
 6400 format ('$ nrad, njump, ngauss, plummer/king (0/1) ... ')
      read (5, *) nrad, njump, ngauss, iking
      write (16, *) ' nrad, njump, ngauss = ', nrad, njump, ngauss
c
      call gauleg (abscis, weight, ngauss)
c
c        The subroutine  GAULEG  will calculate abscissae X_i and
c     weights W_i for the Gauss-Legendre integration with N = 2*NGAUSS
c     (from Press et al, {\it Numerical Recipes}, p. 110).
c     They have been tested - for NGAUSS = 48 - against the values
c     in AS, Table 25.4, pp. 916-919.  GAULEG returns the values:
c
c     ABSCIS(i) = {X_i}^2
c     WEIGHT(i) = 2 * W_i
c
c     which are appropriate for integrals of the form:
c
c     I = \int_a^b {f(t) \, dt \over \sqrt(t-a)} \, ;
c
c     in fact we have
c
c     I \approx \sqrt{b-a} \cdot \sum_{i=1)^n WEIGHT(i) \, f(t_i) \, ,
c
c     with  t_i = a + ABSCIS(i) * (b-a)
c
c     (derived from AS 25.4.36, with a simple change of variable).
c
c--------------------------------------------------------------------------
c
      if (iking .ne. 1) call plummer (radius, dens, poten, emint,
     >                                               thdis, nrad, b)
      if (iking .eq. 1) call king (radius, dens, poten, emint,
     >                                               thdis, nrad, b)
      write (6, *) ' model complete'
c
c        Define the `corrected' density  DENS1  by
c
c               DENS1 = DENS * (1 + (B*R)**2)
c
c     B  is the anisotropy parameter = 1/RA (in Merritt's notation).
c     B=0 for isotropic models
c
      do 300 i = 1, nrad
            dens1(i) = dens(i) * (1.d0 + (b * radius(i))**2)
  300 continue
c
c        Now prepare the interpolation subroutine for the density as
c     a function of the potential
c
      call icsccu (poten, dens1, nrad, cdens1, idens1, ier)
c
c        In order to calculate the distribution function we will need
c     first the value of the integral
c
c     \int_Q^0 {dU \over \sqrt{U-Q)} {d \rho1 \over dU}
c
c     This we will calculate for values of  Q  coinciding
c     with tabulated values of the potential.  The integral will
c     be performed by Gauss's method (see beginning).
c
      idistr = 0
      do 500 i = 1, nrad, njump
            idistr = idistr + 1
            index(idistr) = i
            qplus(idistr) = poten(i)
            width = - qplus(idistr)
c
c        lower limit of integration and width of the interval
c
            value(idistr) = 0.d0
            distr2(idistr) = 0.d0
            if (i .eq. nrad) go to 500
            do 410 j = 1, ngauss
                  aux1(j) = qplus(idistr) + abscis(j) * width
  410       continue
c
c     AUX1  now contains the abscissae of the points where the function
c     has to be evaluated;  AUX2  will contain the values (of d\rho/dU);
c     AUX3  those of the second derivative.
c
            npri = ngauss
            nsec = ngauss
            call dcsevu (poten, dens1, nrad, cdens1, idens1, aux1,
     >                               aux2, npri, aux3, nsec, ier)
c
c        values available; sum up for the integral
c
            do 450 j = 1, ngauss
                  value(idistr) = value(idistr) + aux2(j) * weight(j)
                  distr2(idistr) = distr2(idistr) + aux3(j) * weight(j)
  450       continue
c
c        correct for the scale and for the factor  {\sqrt{2} \over 4 \pi}
c
            value(idistr) = value(idistr) * sqrt(2.d0 * width) /
     >                                              (4.d0 * pi**2)
            distr2(idistr) = distr2(idistr) * sqrt(2.d0 * width) /
     >                                              (4.d0 * pi**2)
c
c        fill in the vectors reporting density, mass and theoretical
c     distr funct at selected points
c
  500 continue
      ndistr = idistr
      write (6, *) ' distribution function computed'
c
c        The distribution function  f  is given by
c
c     f = {\sqrt 2 \over 4 \pi} {d \over d Q_+} {the previous integral}
c
c     Again, splines will do the job for the differentiation.
c
      call icsccu (qplus, value, ndistr, cvalue, ivalue, ier)
      do 600 idistr = 1, ndistr
            aux1(idistr) = qplus(idistr)
  600 continue
c
      npri = ndistr
      nsec = 0
      call dcsevu (qplus, value, ndistr, cvalue, ivalue, aux1,
     >                            distr, npri, aux2, nsec, ier)
c
      do 700 idistr = 1, ndistr
            i = index(idistr)
            write (16, 6500) radius(i), dens(i),
     >        emint(i), qplus(idistr), value(idistr),
     >        distr(idistr), distr2(idistr), thdis(i)
  700 continue
 6500 format (1x, 8f15.9)
      stop
      end
c
c-----------------------------------------------------------------------
c
      subroutine gauleg (abscis, weight, ngauss)
c
c        Finds the abscissae and weights for the Gauss-Legendre
c     integration scheme.  To be used in DISTR.
c
      implicit real*8 (a-h, o-z)
      dimension x(2000), w(2000), abscis(ngauss), weight(ngauss)
      n = 2 * ngauss
c
      m = (n+1) / 2
c
      do 12 i = 1, m
         z = cos (3.141592654d0 * (i-.25d0) / (n+.5d0))
    1    continue
            p1 = 1.d0
            p2 = 0.d0
            do 11 j = 1, n
                  p3 = p2
                  p2 = p1
                  p1 = ((2.d0*j - 1.d0) * z * p2 - (j-1.d0)*p3) / j
   11       continue
            pp = n * (z*p1 - p2) / (z*z - 1.d0)
            z1 = z
            z = z1 - p1 / pp
         if (abs(z-z1) .gt. 5.d-16) go to 1
         x(i) = -z
         x(n+1-i) = z
         w(i) = 2.d0 / ((1.d0 - z*z) * pp*pp)
         w(n+1-i) = w(i)
         abscis(ngauss+1-i) = z*z
         weight(ngauss+1-i) = 2 * w(i)
   12 continue
c     do 15 i = 1, ngauss
c     write (16, 6100) i, abscis(i), weight(i)
c  15 continue
c6100 format (1x, i6, 2f25.15)
      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine plummer (radius, dens, poten, emint, thdis, nrad, b)
      implicit double precision (a-h, o-z)
      dimension radius(nrad), dens(nrad), poten(nrad), emint(nrad)
      dimension thdis(nrad)
      parameter (npoint = 10001)
      dimension cdens(npoint, 3)
      data idens / npoint /
      data pi, g / 3.1415926535897932d0, 1.d0 /
      data nexp / 2 /
      write (6, 1000)
 1000 format ('$ rmax, sigma, r0, b (=1/ra) ... ')
      read (5, *) rmax, sigma, r0, b
      rstep = rmax / (nrad - 1)**nexp
      rhocen = 9.d0 * sigma**2 / (2.d0 * pi * g * r0**2)
      do 10 i = 1, nrad
            radius(i) = (i-1)**nexp * rstep
            dens(i) = rhocen / (sqrt(1.d0 + radius(i)**2 / r0**2))**5
   10 continue
c
c        calculate gravitational potential (in CALCPOT)
c
      call calcpot (radius, dens, poten, emint, nrad)
c
c        fill in the theoretical distribution function
c
      thdis0 = sqrt(2.d0) / (378 * pi**3 * g * r0**2 * sigma)
c
      do 300 i = 1, nrad
            qtrue = - 6.d0*sigma**2 / sqrt(1.d0 + (radius(i)/r0)**2)
            qqq = - qtrue / sigma**2
            thdis(i) = thdis0 * ( (sqrt(qqq))**7 * (1.d0-(b*r0)**2)
     >                     + 63.d0/4.d0 * (b*r0)**2 * (sqrt(qqq))**3 )
  300 continue
c
      write (16, 1100) sigma, r0, rhocen, rmax, emint(nrad), b
 1100 format (' Plummer model;  sigma, r0, rho0 = ', 3f12.6,
     > '  rmax, totmas = ', 2f12.6, '  b ='/ 5x, f10.5)
c
      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine calcpot (radius, dens, poten, emint, nrad)
      implicit double precision (a-h, o-z)
      dimension radius(nrad), dens(nrad), poten(nrad), emint(nrad)
      parameter (npoint = 10001)
      dimension cdens(npoint, 3)
      data pi, g / 3.1415926535897932d0, 1.d0 /
      data idens / npoint /
c
c        This subroutine will compute the potential and the mass inside
c     r given the density.  It will assume that the density, as given
c     in the vector  DENS  and spline-interpolated, is exact.
c     We will define two auxiliary functions,  EM  and  TI, being the
c     integrals of  RHO*R*R  and of  RHO*R:
c
c        EM(r) = 4 \pi \int_0^r \rho(s) s^2 ds
c
c        TI(r) = 4 \pi \int_0^r \rho(s) s ds
c
c        U(r) = U(0) + G * (TI(r) - EM(r)/r)
c
c     (the last is obtained from the definition of the potential, after
c     exchanging the order of integration)
c
c     The gravitational potential is defined to be zero at the outer
c     boundary (supposed finite = RADIUS(NRAD)), and negative inside
c     (a physicist's definition).
c
      call icsccu (radius, dens, nrad, cdens, idens, ier)
c
c        Now we know the coefficients of the spline interpolation
c     routine.  The density is described as a polynomial inside each
c     subinterval.  Therefore the integration -- to calculate the
c     quantities  EM  and  TI  -- is straightforward, implying only
c     recalculating the relevant coefficients.
c     In particular, we will define the coefficients  a0-a4 and b0-b5,
c     of the polynomial representation of  \rho \cdot r and of
c     \rho \cdot r^2  respectively.
c
      em = 0.d0
      ti = 0.d0
      poten(1) = 0.d0
      emint(1) = em
c
c        assumes  RADIUS(1) = 0  -- necessary for the interpolation
c     routines to work correctly
c
      do 200 i = 2, nrad
            s = radius(i) - radius(i-1)
            c0 = dens(i-1)
            c1 = cdens(i-1, 1)
            c2 = cdens(i-1, 2)
            c3 = cdens(i-1, 3)
            r = radius(i-1)
            a0 = c0 * r
            a1 = c1 * r + c0
            a2 = c2 * r + c1
            a3 = c3 * r + c2
            a4 =          c3
            b0 = (c0 * r            ) * r
            b1 = (c1 * r + 2.d0 * c0) * r
            b2 = (c2 * r + 2.d0 * c1) * r + c0
            b3 = (c3 * r + 2.d0 * c2) * r + c1
            b4 = (         2.d0 * c3) * r + c2
            b5 =                            c3
            ti = ti + 4 * pi * s * (a0 + s * (a1/2.d0 + s * (a2/3.d0 +
     >              s * (a3/4.d0 + s * a4/5.d0))))
            em = em + 4 * pi * s * (b0 + s * (b1/2.d0 + s * (b2/3.d0 +
     >              s * (b3/4.d0 + s * (b4/5.d0 + s * b5/6.d0)))))
            poten(i) = g * (ti - em / radius(i))
            emint(i) = em
  200 continue
      rout = radius(nrad)
      emtot = emint(nrad)
      poten0 = - poten(nrad)
c
c        Redefine the zero of the potential
c
      do 220 i = 1, nrad
            poten(i) = poten(i) + poten0
c           write (6, *) radius(i), poten(i)
  220 continue
c
      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine king (radius, dens, poten, emint, thdis, nrad, b)
      implicit double precision (a-h, o-z)
      dimension radius(nrad), dens(nrad), poten(nrad), emint(nrad)
      dimension thdis(nrad)
      dimension yy(2), wk(2, 9), comm(24)
      common / KINGPAR / w0, rho0
      external func
      data pi, g / 3.1415926535897932d0, 1.d0 /
      write (6, 1000)
 1000 format ('$ w0, emtot, rc, b (=1/ra) ... ')
      read (5, *) w0, emtot, rc, b
      if (w0 .gt. 0) w0 = -w0
      wstep = abs(w0) / (nrad - 1)
      rho0 = rho(w0)
      do 10 i = 1, nrad
            poten(i) = w0 + wstep * (i-1)
   10 continue
c
c        first solve for the unscaled model (depending on  W0  only)
c
      neqs = 2
      yy(1) = 0.d0
      yy(2) = 2.d0/3.d0
      tol = 1.d-6
      ind = 1
      nwk = 2
      w = w0
      radius(1) = 0.d0
      dens(1) = 1.d0
      emint(1) = 0.d0
      do 100 i = 2, nrad
            wend = poten(i)
            call dverk (neqs, func, w, yy, wend, tol, ind, comm, nwk,
     >                                                         wk, ie)
c
            radius(i) = sqrt(yy(1))
            dens(i) = rho(w) / rho0
            emint(i) = 2.d0 * yy(1) * sqrt(yy(1)) / yy(2)
  100 continue
c
c     now find the scaling factors (from the core radius, given, and the
c     total mass) and calculate the theoretical ISOTROPIC distr func
c
      sigma = sqrt (g * emtot / rc / emint(nrad))
      rhocen = (9.d0 * sigma**2) / (4.d0 * pi * g * rc**2)
      thdis0 = 9.d0 / (4.d0*pi*g*rc**2*sigma * 4.d0*pi*sqrt(2.d0)
     >                                                        * rho0)
c
      write (16, 1100) w0, emtot, rc, sigma, rhocen, b
 1100 format (' King model; w0 = ', f10.5, '  M, rc, sigma, rhocen = ',
     > 4f15.9, 5x, 'b ='/ 5x, f10.5 )
c
      do 200 i = 1, nrad
            radius(i) = radius(i) * rc
            emint(i) = emint(i) * sigma**2 * rc / g
            dens(i) = dens(i) * rhocen 
            thdis(i) = thdis0 * (exp(-poten(i)) - 1)
            poten(i) = poten(i) * sigma**2
c     write (16, 2000) radius(i), emint(i), poten(i), dens(i)
c2000 format (1x, 4f16.10)
  200 continue
      return
      end
c
      double precision function func(n, w, y, yprime)
      implicit double precision (a-h, o-z)
      dimension y(n), yprime(n)
      common / KINGPAR / w0, rho0
c
c        The differential equation to be solved is
c
c     x" - (3/2) (x')^2 + (9/4) (rho(w)/rho(w0) (x')^3 = 0
c
c     where             x = radius^2
c
c        The two components of the vector  Y  are  X, X' respectively
c
      yprime(1) = y(2)
c     write (6, 900) n, w, y(1), y(2)
c 900 format (i10, 3f15.9)
      if (y(1) .gt. 1.d-8) yprime(2) = 1.5d0 *  y(2)*y(2) * (1.d0 -
     >                               1.5 * y(2) * rho(w)/rho0) / y(1)
      if (y(1) .le. 1.d-8) yprime(2) = 0.4d0 * (1.d0 +
     >                            (sqrt(-w0))**3 / rho0)
c
c        Despite the apparent singularity in the first form, the
c     differential equation - with the right initial conditions - has a
c     regular solution.  The second statement amounts to saying that
c                      X" = - (2/5) rho'(w0)/rho(w0)
c
      return
      end
c
c
c
      double precision function rho(w)
      implicit double precision (a-h, o-z)
      data coeff / 0.8862269255d0 /
c
c     COEFF  =  \sqrt{\pi} / 2  -- to convert  ERF  into the integral
c
c
c        The function  RHO  calculates the integral
c
c     rho = [\exp(y^2) \int_0^y \exp(-v^2) dv] - y - 2 y^3 / 3 
c                                                      (y = \sqrt{-w})
c
c     In terms of this, the density is defined as
c
c     dens = 4 \sqrt(2) \pi k \sigma^3 \exp(w0)
c
c     where  k  is the constants in the definition of the distribution
c     function, which can be expressed in terms of the other quantities
c     by using the condition
c
c        dens(w0) = {9 \sigma^2 \over 4 \pi g rc^2},
c
c     so that
c
c        k = {9 \over {(4 \pi)}^2 g rc^2 \sigma \sqrt{2} rho(w0) \exp(w0)}
c
      rho = 0.d0
      if (w .ge. 0.d0) return
      y = sqrt(-w)
      rho = exp(-w) * coeff * erf(y) - y - 2.d0 * y**3 / 3.d0
      return
      end
	subroutine dcsevu(a)
	return
	end

	subroutine icsccu(a)
	return
	end

	subroutine dverk(a)
	return
	end


