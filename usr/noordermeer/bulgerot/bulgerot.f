c     #########################################################################
c                           WELCOME TO BULGEROT  
c
c     This program calculates the rotation curve for an oblate spheroidal 
c     bulge. It assumes a projected surface density that follows a Sersic 
c     profile. The central surface magnitude M_0, characteristic radius R_0 and
c     Sersic index n_sers must be given. Also, the inclination and intrinsic 
c     axis ratio of the bulge must be given. Finally, the distance and the 
c     Galactic foreground extinction of the galaxy is needed.

c     Copyright (C) 2007, Edo Noordermeer
c     E-mail: edo.noordermeer@gmail.com
c
c     Updates:
c     sept-2013:  Minor modification to read from stdin for NEMO's runbulgerot 
c
c     If you use this program for a scientific publication, I would appreciate 
c     a reference to the paper 'The rotation curves of flattened Sersic bulges
c     (Noordermeer 2008)'.
c
c     #########################################################################

      real*8       pi,G,pc,Masssun,MagsunR
      character*20 galaxy
      character*99 logfile
      logical      onscreen 
      real*8       Mu_0_app,R_0,n_sers
      real*8       incdeg,axisratio,distpc,A_R
      real*8       Rstart, dR
      integer      N
      real*8       inclination,distance,arcsectom,arcsectopc
      real*8       Mu_0_abs,I_0,I0,R0,nsers
      real*8       asympmass, mass1, mass2, asympmag_abs
      real*8       asympmag_app
      integer      i
      real*8       radius(1000),intensity(1000),density(1000)
      real*8       velocity(1000), mass(1000)
      real*8       r
      real*8       rho, intens, vc, cummass
      external     rho, intens, vc, cummass, qromb, qromo, midinf
      external     integrand3

      common /sersic/ I0,R0,nsers
      common /geometry/ inclination,axisratio
      common /univconst/ pi,G,pc,Masssun,MagsunR



c     #########################################################################
c     ## DEFINE CONSTANTS: CHANGE FOR EACH GALAXY #############################
c     #########################################################################

      if (.true.) then
         write (*,*) 'Interactive input: (5 lines expected)'
         read (*,*) logfile
         read (*,*) galaxy
         read (*,*) Mu_0_app, R_0, n_sers
         read (*,*) incdeg,axisratio,distpc,A_R
         read (*,*) Rstart,dR,N
         onscreen = .false.
      else
c     Define logfile
      logfile = 'bulgerot.dat'
      onscreen = .true.
c     Define galaxy name
      galaxy = 'thisgalaxy'

c     Define Sersic parameters: 
c     apparent central R-band surface magnitude M_0, 
c     characteristic radius R_0 in arcseconds, 
c     Sersic index n_sers
      Mu_0_app = 15.0d0
      R_0      = 1.0d0
      n_sers   = 1.0d0 

c     Define geometric parameters: 
c     inclination in degrees, 
c     intrinsic axis ratio, 
c     distance in pc, 
c     galactic foreground extinction
      incdeg    = 0.0d0
      axisratio = 1.0d0
      distpc    = 10.0d6
      A_R       = 0.0d0

c     Define radii to calculate quantities
c     Starting radius, increment, and number of points (all in kpc)
      Rstart = 0.0d0
      dR     = 0.01d0
      N      = 250
      endif

c     open outputfile for writing
      open (11,file=logfile)
     
c     #########################################################################



c     #########################################################################
c     ## NO NEED TO CHANGE ANYTHING BELOW #####################################
c     #########################################################################

c     #########################################################################
c     Define universal constants:
c     pi; G; 1 parsec; 1 solar mass; absolute R-band solar magnitude
      pi      = 2.0d0 * asin(1.0d0)
      G       = 6.672d-11
      pc      = 3.08567802d16
      Masssun = 1.989d30
      MagsunR = 4.28d0                     
c     #########################################################################


c     #########################################################################
c     Convert values to physical values
c     inclination in radians
c     distance in meters
c     1 arcsec on sky in meters in galaxy
      inclination = incdeg/180.0d0*pi
      distance    = distpc*pc
      arcsectom   = 4.847*pc*(distpc/1.0d6)
      arcsectopc  = 4.847*(distpc/1.0d6)
 
c     absolute central R-band surface magnitude
c     central surface density in Msun/"^2
      Mu_0_abs = Mu_0_app - 5.0d0*(log(distpc)/log(10.0d0)) + 5.0d0 
     $     - A_R
      I_0 = 10.0d0**( (Mu_0_abs - MagsunR)/(-2.5d0) )

c     central surface density in kg/m^2.0
c     characteristic radius in m
c     Sersic index
      I0 = I_0 * Masssun/(arcsectom**2.0d0)
      R0 = R_0 * arcsectom
      nsers = n_sers

c     Total mass, integrated to infinity
      call QROMB(integrand3, 0.0d0, 1.0d20, mass1)
      call QROMO(integrand3, 1.0d20,1.0d60, mass2, midinf)
      asympmass = (mass1+mass2)/Masssun
      asympmag_abs = -2.5d0* (log(asympmass)/log(10.0d0)) + MagsunR
      asympmag_app = asympmag_abs + 5.0d0*(log(distpc)/log(10.0d0)) - 
     $     5.0d0 + A_R
c     #########################################################################


c     #########################################################################
c     Write header to screen
      if (onscreen) then
      write (*,fmt=9990) '-----------------------------------',
     $        '---------------------------------------------'
      write (*,*) 'Calculating rotation curve for Sersic bulge of ',
     $     galaxy
      write (*,fmt=9990) '-----------------------------------',
     $        '---------------------------------------------'

      write (*,*)
      write (*,fmt=9995) 'Central R-band surface magnitude: ', 
     $     'Mu_0 = ',Mu_0_app
      write (*,fmt=9996) 'Characteristic radius:            ',
     $     'R_0 = ', R_0, '"'
      write (*,fmt=9995) 'Sersic index:                     ',
     $     'n = ', n_sers
      write (*,*)

      write (*,fmt=9994) 'Assumed inclination:              ',
     $     'i = ', incdeg, 'deg'
      write (*,fmt=9995) 'Assumed intrinsic axis-ratio:     ',
     $     'B/A = ', axisratio
      write (*,fmt=9994) 'Assumed distance:                 ',
     $     'D = ', distpc/1.0d6, 'Mpc'
      write (*,*)

      write (*,fmt=9993) 'Total integrated mass:            ',
     $     'M_inf = ', asympmass, 'Msun'
      write (*,fmt=9995) 'Total apparent R-band magnitude:  ',
     $     'Mu_inf = ', asympmag_app
      write (*,*)

      write (*,fmt=9992) '1 arcsec corresponds to           ',
     $     arcsectopc, 'pc'
      write (*,fmt=9992) '1 kpc corresponds to              ',
     $     1.0D3/arcsectopc, '" '
      write (*,*)
      write (*,fmt=9990) '-----------------------------------',
     $        '---------------------------------------------' 
      write (*,fmt=9989)         'radius', 'intensity', 'density ', 
     $        'mass  ', 'velocity', 'progress'
      write (*,fmt=9989)         '(kpc) ', '(kg/m^2) ', '(kg/m^3)',
     $        '(Msun) ', ' (km/s) ', '  (%)  '
      write (*,fmt=9990) '-----------------------------------',
     $        '---------------------------------------------'
      write (*,fmt=9991) 0.0,0.0,0.0,0.0,0.0,0.0,'%'
      endif

c     Write header to file
      write (*,*) 'Writing to ',logfile

      write (11,9997) '----------------------------------',
     $        '------------------------------------'
      write (11,'(A,A)') '# Rotation curve for Sersic bulge of ',galaxy
      write (11,9997) '----------------------------------',
     $        '------------------------------------'     

      write (11,*)
      write (11,fmt=9995) 'Central R-band surface magnitude: ', 
     $     'Mu_0 = ',Mu_0_app
      write (11,fmt=9996) 'Characteristic radius:            ',
     $     'R_0 = ', R_0, '"'
      write (11,fmt=9995) 'Sersic index:                     ',
     $     'n = ', n_sers
      write (11,*)

      write (11,fmt=9994) 'Assumed inclination:              ',
     $     'i = ', incdeg, 'deg'
      write (11,fmt=9995) 'Assumed intrinsic axis-ratio:     ',
     $     'B/A = ', axisratio
      write (11,fmt=9994) 'Assumed distance:                 ',
     $     'D = ', distpc/1.0d6, 'Mpc'
      write (11,*)

      write (11,fmt=9993) 'Total integrated mass:            ',
     $     'M_inf = ', asympmass, 'Msun'
      write (11,fmt=9995) 'Total apparent R-band magnitude:  ',
     $     'Mu_inf = ', asympmag_app
      write (11,*)

      write (11,fmt=9992) '1 arcsec corresponds to           ',
     $     arcsectopc, 'pc'
      write (11,fmt=9992) '1 kpc corresponds to           ',
     $     1.0D3/arcsectopc, '" '
      write (11,*)
      write (11,9997) '----------------------------------',
     $        '------------------------------------'          
      write (11,fmt=9998)         'radius', 'intensity', 'density ', 
     $        'mass  ', 'velocity'
      write (11,fmt=9998)         '(kpc) ', '(kg/m^2) ', '(kg/m^3)',
     $        '(Msun) ', ' (km/s) '
      write (11,fmt=9997) '----------------------------------',
     $        '------------------------------------'
      write (11,fmt=9999) 0.0,0.0,0.0,0.0,0.0
c     #########################################################################


c     #########################################################################
c     Calculate for a range of radii: 
c     the surface density in kg/m^2
c     the 3D mass density in kg/m^3
c     the cumulative mass inside radius R in solar masses
c     the rotation velocity in km/s
c     Write output to screen and to file.
      do i=1,N
         radius(i) = Rstart + i*dR
         r = radius(i)*1.0d3*pc

         intensity(i) = intens(r)
         density(i) = rho(r)
         mass(i) = cummass(r)
         velocity(i) = vc(r)         
         if (onscreen) then
         write (*,9991)      radius(i),intensity(i),density(i),
     $        mass(i),velocity(i),1.0d2*i/N, '%'
         endif
         write (11,fmt=9999) radius(i),intensity(i),density(i),
     $        mass(i),velocity(i)
      enddo
c     #########################################################################


c     #########################################################################
c     Wrap up
      close (11)

 9989 format ('#',A10,4A14,A12)
 9990 format ('#',A35,A45)
 9991 format (F12.6,4D14.6,F8.2,A1)
 9992 format ('#',A48,F12.6,A3)
 9993 format ('#',A40,A9,D12.6,A5)
 9994 format ('#',A40,A9,F12.6,A4)
 9995 format ('#',A40,A9,F12.6)
 9996 format ('#',A40,A9,F12.6,A2)
 9997 format ('#',A34,A36)
 9998 format ('#',A10,4A14)
 9999 format (F12.6,4D14.6)
      end
c     #########################################################################

      
c     -------------------------------------------------------------------------
      real*8   function vc(R)
      real*8   Rprime, inclination, axisratio
      real*8   constant, integral,R
      real*8   pi, G, pc, Masssun, MagsunR
      external QROMO, integrand2, midpnt

      common /Rradius/ Rprime
      common /geometry/ inclination,axisratio
      common /univconst/ pi,G,pc,Masssun,MagsunR

      Rprime = R

      constant =  4.0d0*pi*G*axisratio
      
      call QROMO(integrand2,0.0d0,R,integral,midpnt)

      vc = sqrt(constant*integral)/1.0d3


      return
      end
c     -------------------------------------------------------------------------



c     -------------------------------------------------------------------------
      real*8   function integrand2(m)
      real*8   inclination, axisratio, eccentricity2, Rprime
      real*8   m, teller, noemer, rho
      real*8   pi, G, pc, Masssun, MagsunR
      external rho

c     This function calculates the integrand 
c     rho*m^2/sqrt(R^2 - m^2*e^2)


      common /geometry/ inclination,axisratio
      common /Rradius/ Rprime
      common /univconst/ pi,G,pc,Masssun,MagsunR

      eccentricity2 = 1.0d0 - axisratio**2.0d0
c     This is the square of the eccentricity

      teller = rho(m)*(m**2.0d0)
      noemer = sqrt( Rprime**2.0d0 - (m**2.0d0)*eccentricity2 )

      integrand2 = teller/noemer
     
      return
      end
c     -------------------------------------------------------------------------



c     -------------------------------------------------------------------------
      real*8   function cummass(R)
      real*8   R, mass
      real*8   pi, G, pc, Masssun, MagsunR
      external qromb, integrand3, midpnt

      common /univconst/ pi,G,pc,Masssun,MagsunR

      call QROMB(integrand3, 0.0d0, R, mass)

      cummass = mass/Masssun

      return
      end
c     -------------------------------------------------------------------------



c     -------------------------------------------------------------------------
      real*8   function integrand3(q)
      real*8   q, intens
      real*8   inclination, axisratio, f2, BA
      real*8   pi, G, pc, Masssun, MagsunR
      external intens

      common /geometry/ inclination,axisratio
      common /univconst/ pi,G,pc,Masssun,MagsunR

      f2 = axisratio**2.0d0
      BA = sqrt(f2 + (1.0d0 - f2)*(cos(inclination))**2.0d0)

      integrand3 = 2.0d0*pi*BA*intens(q)*q

      return
      end
c     -------------------------------------------------------------------------



c     -------------------------------------------------------------------------
      real*8  function intens(q)
      real*8  I0,R0,nsers
      real*8  q, qrel
      real*8  pi, G, pc, Masssun, MagsunR

      common /sersic/ I0,R0,nsers
      common /univconst/ pi,G,pc,Masssun,MagsunR

      qrel = q/R0
      intens = I0*exp( -(qrel**(1.0d0/nsers)) )
      
      return
      end
c     -------------------------------------------------------------------------



c     -------------------------------------------------------------------------
      real*8    function rho(m)
      real*8    mprime, inclination, axisratio
      real*8    m, constant, integral, integrala, integralb
      external  QROMO, integrand1,midexp, midinf
      real*8    pi, G, pc, Masssun, MagsunR

      common  /mradius/ mprime
      common  /geometry/ inclination,axisratio
      common /univconst/ pi,G,pc,Masssun,MagsunR

      mprime = m

      constant = sqrt( (sin(inclination))**2.0 + 
     $               ( (cos(inclination))/axisratio )**2.0) / PI

      call QROMO(integrand1, m, 1.0d60, integral, midinf)

      rho = -1.0d0*constant*integral

      return
      end
c     -------------------------------------------------------------------------



c     -------------------------------------------------------------------------
      real*8  function integrand1(q)
      real*8  I0,R0,nsers,mprime
      real*8  q, qrel, term1,term2,term3,term4
      real*8    pi, G, pc, Masssun, MagsunR

c     This function calculates the integrand didq/sqrt(q^2 - m^2)

      common /sersic/ I0,R0,nsers
      common /mradius/ mprime
      common /univconst/ pi,G,pc,Masssun,MagsunR

      qrel = q/R0
      term1 = -I0/(nsers*R0)
      term2 = exp( -(qrel**(1.0d0/nsers)) )
      term3 = qrel**(1.0d0/nsers - 1.0d0)
      term4 = sqrt(q**2.0d0 - mprime**2.0d0)
      
      integrand1 = term1*term2*term3/term4

      return
      end
c     -------------------------------------------------------------------------



c     -------------------------------------------------------------------------
      SUBROUTINE qromo(func,a,b,ss,choose)
      INTEGER    JMAX,JMAXP,K,KM
      REAL*8     a,b,func,ss,EPS
      EXTERNAL   func,choose
      PARAMETER  (EPS=1.d-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint
      INTEGER    j
      REAL*8     dss,h(JMAXP),s(JMAXP)

      h(1)=1.
      do 11 j=1,JMAX
        call choose(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=h(j)/9.
11    continue
      pause 'too many steps in qromo'
      END
c     -------------------------------------------------------------------------


c     -------------------------------------------------------------------------
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER    n,NMAX
      REAL*8     dy,x,y,xa(n),ya(n)
      PARAMETER  (NMAX=10)
      INTEGER    i,m,ns
      REAL*8     den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
c     -------------------------------------------------------------------------


c     -------------------------------------------------------------------------
      SUBROUTINE midexp(funk,aa,bb,s,n)
      INTEGER    n
      REAL*8     aa,bb,s,funk
      EXTERNAL   funk
      INTEGER    it,j
      REAL*8     ddel,del,sum,tnm,x,func,a,b

      func(x)=funk(-log(x))/x
      b=exp(-aa)
      a=0.
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END
c     -------------------------------------------------------------------------


c     -------------------------------------------------------------------------
      SUBROUTINE midinf(funk,aa,bb,s,n)
      INTEGER n
      REAL*8 aa,bb,s,funk
      EXTERNAL funk
      INTEGER it,j
      REAL*8 a,b,ddel,del,sum,tnm,func,x
      func(x)=funk(1./x)/x**2
      b=1./aa
      a=1./bb
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END
c     -------------------------------------------------------------------------


c     -------------------------------------------------------------------------
      SUBROUTINE qromb(func,a,b,ss)
      INTEGER JMAX,JMAXP,K,KM
      REAL*8 a,b,func,ss,EPS
      EXTERNAL func
      PARAMETER (EPS=1.d-8, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint,trapzd
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call trapzd(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25*h(j)
11    continue
      pause 'too many steps in qromb'
      END
c     -------------------------------------------------------------------------


c     -------------------------------------------------------------------------
      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      REAL*8 a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*8 del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END
c     -------------------------------------------------------------------------


c     -------------------------------------------------------------------------
      SUBROUTINE midpnt(func,a,b,s,n)
      INTEGER n
      REAL*8 a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*8 ddel,del,sum,tnm,x
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END
