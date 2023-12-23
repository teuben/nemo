c***********************************************************************
c+      Keywords for NEMO retrieved with the ftoc(1NEMO) program.
c   x=1\n                   x coordinate
c   y=0\n                   y coordinate
c   VERSION=1.0\n           24-sep-91 PJT
c-      
c-----------------------------------------------------------------------
      SUBROUTINE nemomain
      IMPLICIT REAL*8(a-z)
      include 'henyey.h'
      INTEGER n,n1,n2,l
      LOGICAL out
c                   (LUN1 is for input, LUN2 for output) 
      lun1=5
      lun2=6
c                   (Get stuff from the commandline via NEMO)
      CALL userinp(1)
c                   (some constants, kept in the /CONSTS/ common in model.h)
      pi=3.14159265358979
      gravc=0.4298
      ome2=ome*ome
c                   (MASS, A, B are those of the bar)
      mass = 2.29
      mass=mass*1.33333333/2.0
      a=0.5
      b=0.49
c                   MDi mass of a disk, ADi length scale of a disk (i=1,2)
c                   Note MDi is actually M*GRAVC
      md1=9.9
      ad1=1.414
      md2=1.23
      ad2=0.106
c--BEGIN TEST - reset some variables for known analytical answers
      IF(.TRUE.) THEN
         mass=1.0
         md1=0.0
         ad1=2.0
         md2=0.0
         ad2=1.0
      ENDIF
c--END TEST
c
      CALL bar1(x0,y0,pot,fx,fy,fxx,fxy,fyy)
      den = fxx + 2.0*fyy
      WRITE(lun2,*) x0,y0,pot,fx,fy,fxx,fxy,fyy,den

      END
      SUBROUTINE bar1(x,y,pot,fx,fy,fxx,fxy,fyy)
      IMPLICIT REAL*8(a-z)
      include 'model.h'
c
c       Homogeneous prolate spheroid with axes A,B and mass MASS
c       Potential is defined > 0
c       BUG: can't handle spherical (a.EQ.b) cases
c
c                           If zero mass, quick exit
      IF(mass.EQ.0.0) THEN
         pot=0.0
         fx=0.0
         fy=0.0
         fxx=0.0
         fxy=0.0
         fyy=0.0
         RETURN
      ENDIF
      gm=mass*gravc
      ecc=SQRT(a*a-b*b)/a
      ae=a*ecc
      ae3=ae*ae*ae
      ale=LOG((1.0-ecc)/(1.0+ecc))
      aint1=(-2.0*ecc-ale)/ae3
c      omcr=SQRT(1.5*gm*aint1)         ! redundant for now ...
      r1=SQRT((x+ae)*(x+ae)+y*y)      
      r2=SQRT((x-ae)*(x-ae)+y*y)      
      IF(r1+r2.LE.2.0*a) THEN
c                                         (x,y) interior point
         al=0.0
         ale=LOG((1.0-ecc)/(1.0+ecc))
         aint1=(-2.0*ecc-ale)/ae3
         aint2=(ecc/(1.0-ecc*ecc)+0.5*ale)/ae3
         aint3=-ale/ae
         aux1=0.0
         aux2=0.0
         aux3=0.0
c                                         (x,y) exterior point
      ELSE
         bb=a*a+b*b-x*x-y*y
         cc=a*a*b*b-b*b*x*x-a*a*y*y
         al=-0.5*bb+0.5*SQRT(bb*bb-4.0*cc)
         rtal=SQRT(a*a+al)
         ale=LOG((rtal-ae)/(rtal+ae))
         aint1=(-2.0*ae/rtal-ale)/ae3
c         aint2=(ae*rtal/(b*b+al)+0.5*ale)/ae3
         aint2=(ae*rtal/(bb*bb+al)+0.5*ale)/ae3
         aint3=-ale/ae
         daldx=2.0*x*(al+b*b)/(2.0*al+bb)
         daldy=2.0*y*(al+a*a)/(2.0*al+bb)
         aux1=1.5*gm*x*daldx/(rtal*(a*a+al)*(b*b+al))
         aux2=1.5*gm*y*daldy/(rtal*(b*b+al)*(b*b+al))
         aux3=1.5*gm*x*daldy/(rtal*(a*a+al)*(b*b+al))
      ENDIF
      pot=-x*x*aint1-y*y*aint2+aint3
      pot=0.75*gm*pot
      fx=-1.5*gm*x*aint1
      fy=-1.5*gm*y*aint2
      fxx=-1.5*gm*aint1+aux1
      fyy=-1.5*gm*aint2+aux2
      fxy=              aux3
      
      RETURN
      END
c***********************************************************************           
      SUBROUTINE disk2(gm,a,x,y,psi,psix,psiy,psixx,psixy,psiyy)
      IMPLICIT REAL*8(a-z)
c
c Toomre disk no. 2
c gm = gravc*mass     a = scale length
c potential psi > 0   fx = dspidx = psix
c   BUG: can't handle the origin
c
c                           If zero mass, quick exit
      IF(gm.EQ.0.0) THEN
         psi=0.0
         psix=0.0
         psiy=0.0
         psixx=0.0
         psixy=0.0
         psiyy=0.0
         RETURN
      ENDIF
      r2=x*x+y*y
      r=SQRT(r2)
      r3=r*r2
      a2=a*a
      a2r2=a2+r2
      wa2r2=SQRT(a2r2)
      psi=gm/wa2r2
      psir=-gm*r/(a2r2*wa2r2)
      psirr=-gm*(a2-2.0*r2)/(a2r2*a2r2*wa2r2)
      psix=psir*x/r
      psiy=psir*y/r
      psixx= psir*y*y/r3+psirr*x*x/r2
      psixy=-psir*x*y/r3+psirr*x*y/r2
      psiyy= psir*x*x/r3+psirr*y*y/r2
      RETURN
      END
c***********************************************************************            
      SUBROUTINE userinp(n)
      IMPLICIT REAL*8(a-h,o-z)
      INTEGER n
c
c  This routine obtains parameters from the command line and
c  places them in the common blocks for later access 
c
      include 'henyey.h'
      DOUBLE PRECISION getdparam
      INTEGER getiparam

      WRITE(lun2,'(''REAL*8 version'')')
      IF (n.EQ.1) THEN
c                                NEMO user interface
         x0 = getdparam('x')
         y0 = getdparam('y')
         ome = 0.0
      ELSE
c                               do your own READ and WRITE here...
         WRITE(lun2,'(''Enter: x0,y0,u0,v0,per,type,ome,norbit,step'')')
         READ(lun1,*) x0,y0,u0,v0,per,type,ome,norbits,step
      ENDIF

      RETURN
      END      
      
