c CHANGES 25-FEB-1994: (KK)
c    Core radius of halo is now calculated from halo density only.
c    Before, the King radius was taken to depend on the total density.
c    Makes no difference if coreparam is set to zero.
c
c    R1 (the radius which appears for dimensional reasons in the halo
c    DF) has been reinstated (before it was set to 1.). 
c    It is set as follows:
c    The radius ra is prompted for: it is such that the halo rotation
c    curve near r=0 is (r/ra)*v0 
c    This radius is then used to calculate the central density of the
c    halo, which then gives r1. (accounting only for the C-term with
c    q=1 in the halo DF, i.e. taking the isothermal sphere halo, and 
c    ignoring truncation effects).
c
c    Convergence has been increased greatly by using a weighted mean of
c    old and new potential as start for the next step, rather than the
c    new potential only. .75*new+.25*old seems to be very effective
c    in damping waves which sometimes took a long time to settle down
c    before. The first time a particular harmonic is introduced its full
c    value is used, though.
c
c    Now also write out the second R-derivative of the potential, 
c    to enable a good calculation of the epicycle frequency. This is read
c    by getfreqs3.
c
c  INCORPORATED WITH THE CHANGES OF JD 21-MAR-1994
c
c CHANGES 21 MARCH 1994 (KK)
c  
c  disk cutoff radius implemented as a sharp cutoff in angular momentum.
c  This corresponds to multiplying the disk density by a complementary 
c  error function.
c  At every iteration the disk density now changes slightly, though, since
c  the density near the cutoff depends on omega and kappa.
c
c  The output file now also includes a table of disk surface mass density vs. 
c  radius. The disk parameters as input are also to be found there. 
c  Getfreqs5 reads them OK.
c
c  CHANGES 11 MAY 1994 (KK): 
c        
c  Changes the disk component to a dynamic one, which depends more on 
c  the potential. Rather than enforcing a sech**2 disk, which
c  is difficult to realize selfconsistently in the presence of other
c  gravitating components (the disk is then no longer isothermal in 
c  the vertical direction), the disk is now isothermal, with density
c  sur0/(2*zdisk)*exp(-r/rdisk)*exp(-(psi(r,z)-psi(r,0))/sigz2(r)). 
c  This makes the vertical kinematics closer to isothermal.
c  The vertical dispersion is chosen to give the correct half-density 
c  scale height: sigz2=[psi(r,0)-psi(r,zdisk)]/ln[sech**2(1)]
c  Instead of writing the disk surface density, now write the disk
c  midplane density to the dbh.dat file.
c
c  Also write out the potential at every iteration on a series of 10
c  grid points, in the file 'refpoints.dat'. This can be used 
c  afterwards to see how the different harmonics contribute.
c  The grid lies in the disk equator, at half a scale height above it,
c  and at 2 scale height above it.
c
c  CHANGES FROM DBH6 TO DBH7: (KK 13 MAY 1994)
c
c  Try to reduce the load on the higher harmonics by subtracting off
c  a potential obtained by integrating a sech**2 law vertically.
c  Might give problems near the cutoff radius if the truncation
c  width is very narrow (on the order of the scale height).
c
c  Have also simplified the disk cutoff ---- simply use an error 
c  function in radius. Makes life a lot easier for calculating the
c  density corresponding to the approximate disk potential...
c  
c  The parameter sigrtrunc has been turned into drtrunc, the gaussian
c  interval over which the truncation takes place.
c
c  FUNDAMENTAL SNAG: THE POTENTIAL HARMONICS NOW NO LONGER HAVE 
c  THE CORRECT ASYMPTOTIC BEHAVIOUR AT INFINITY. HENCE THE SERIES
c  EXPANSION BECOMES KINDA BOGUS....
c
c  CHANGES 19 MAY 1994 (KK)
c
c  Can be remedied by writing the potential in spherical coordinates.
c  So           psi=f(r) ln cosh r cos(theta)/zdisk
c  which at large radii does tend to zero.
c  f(r) again is such that the midplane density is correct.
c
c  ALSO CHANGES THE FILE NAMES FOR THE OUTPUT, TO DBH.DAT, H.DAT AND B.DAT
c
c  Output file dbh.dat now contains all input parameters, i.e. all
c  halo, disk and bulge parameters.

      program dbh

      common /potconstants/ apot(20,0:20000), fr(20,0:20000), 
     +     dr, nr, lmax, potcor
      common /legendre/ plcon(0:40)
      common /gparameters/  a, b, c, v0, q, psi00, 
     +                      psiout, rho1, sigbulge2, 
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /moreconstants/ v02, v03, rdisk2, diskconst, bulgea
      common /flags/ idiskflag, ibulgeflag, ihaloflag

      parameter(pi=3.1415926535)
      real adens(20,0:20000),s1(0:20000),s2(0:20000)
      real rref(10),zref(10),pref(10),oldpref(10)
      real fr2(20,0:20000)
      character ans*1
      
c  constants for legendre polynomials
      do l=0,40
         plcon(l) = sqrt((2*l+1)/(4.0*pi))
      enddo
c  set flags
      idiskflag = 0
      ibulgeflag = 0
      ihaloflag  = 0
c set default parameters
c disk: 
      rmdisk = 1.0
      rdisk  = 1.0
      zdisk  = .15
      outdisk = 5.0
      drtrunc = .3
c bulge:
      rho1 = 5
      psiout = -2.0
      sigbulge = .4
c halo:
      psi00 = -3
      q = 1.0
      v0 = 1.0
      rking = 1
c  
c  the core parameter is the size of the core radius in the DF (Evans's R_c)
c  as a multiple of the King radius, all squared. In the spherical case, 
c  a King model corresponds to coreparam=0.
c  The halo Ra radius is the radius at which the halo rotation curve, at its
c  initial slope ignoring cutoffs and the other components, reaches v0.
c  
c
c Enter halo parameters
      write(*,*) 'Include halo?'
      read(*,'(a)') ans

      if( ans .eq. 'y' ) then
         write(*,*) 'central potential, v0, q, 
     +               coreparam (=Rc/rK)^2, halo Ra?'
         read(*,*) psi00,v0,q,coreparam,ra
         ihaloflag = 1
         v02 = v0**2
         v03 = v0**3
         rhochalo=3*v02/(4*pi*ra**2)
         write(*,*) 'rhochalo = ', rhochalo
         r1=v0/sqrt(4*pi*rhochalo)*exp(-psi00/v02)
c  fix Evans's B constant such that the core radius R_c is 
c  prop. to the King radius of the model.
         a=(2/pi)**2.5*(1-q**2)/q**2/v03/r1**4
         b=0
         c=(2*q**2-1)/(4*pi**2.5*q**2*v0)/r1**2
      endif
      
c Enter disk parameters
      write(*,*) 'Include disk?'
      read(*,'(a)') ans
      if( ans .eq. 'y' ) then
         write(*,*) 'Disk mass, scale length, radius, 
     +	             scale height and trunc. width?'
         read(*,*) rmdisk, rdisk, outdisk, zdisk, drtrunc
         idiskflag = 1

         rdisk2 = rdisk*rdisk
         diskconst = rmdisk/(4.0*pi*rdisk2*zdisk)
         do iref=1,4
           rref(iref)=iref*rdisk*2
           zref(iref)=0
           enddo
         do iref=5,7
           rref(iref)=(iref-5)*rdisk*4
           zref(iref)=zdisk/2.
           enddo
         do iref=8,10
           rref(iref)=rref(iref-3)
           zref(iref)=zdisk*2.
           enddo
         open(18,file='refpoints.dat',status='unknown')
      endif


c Enter bulge parameters
      write(*,*) 'Include bulge?'
      read(*,'(a)') ans
      if( ans .eq. 'y' ) then
         write(*,*) 'Bulge central density (ish), cutoff potential, 
     +               velocity dispersion?'
         read(*,*) rho1,psiout,sigbulge
         sigbulge2=sigbulge**2
         bulgea=rho1*exp((psi00-psiout)/sigbulge2)
         ibulgeflag = 1
      endif

      if( ihaloflag .eq. 1 ) then 
         dr=1
         apot(1,0)=psi00*sqrt(4.*pi)
         lmax=0
         nr=2
         dens00=totdens(0.,0.)
         rking=sqrt(9*v0**2/(8*pi*dens00))
c  find Rc which gives desired Rc/RK by bisection
         if (coreparam.le.0.) then
            b=0
         else
            rclo=0
            rchi=rking*sqrt(coreparam)
            rcmid=rchi/2.
            do while (rchi-rclo.gt.1.e-5*rcmid)
               b=sqrt(2/pi**5)*rcmid**2/q**2/v0/r1**4
               dens00=totdens(0.,0.)
               rcnew=sqrt(coreparam)*sqrt(9*v02/(8*pi*dens00))
               if (rcnew.gt.rcmid) then
                  rclo=rcmid
               else
                  rchi=rcmid
               endif
               rcmid=(rclo+rchi)/2.
            enddo
            b=sqrt(2/pi**5)*rcmid**2/q**2/v0/r1**4
            dens00=totdens(0.,0.)
            rking=sqrt(9*v02/(8*pi*dens00))
         endif
c  report final King radius and central density         
         write(*,*) 'central density=',dens00
         write(*,*) 'King radius=',rking
      endif
      
c  cold start
      write(*,*) 'radial step size, no. of bins?'
      read(*,*) dr,nr
      write(*,*) 'max. azimuthal harmonic l?'
      read(*,*) lmaxx
c  initial potential is a spherical approximation to the unlowered Evans model.
      if( ihaloflag .eq. 0 ) then
         q = 1
         psi00 = -3.0
         rcmid = 1.0 
         v02 = 1.0
      endif
      z = 0.0
      do ir=0,nr
         r=ir*dr
c         apot(1,ir)=(psi00+0.5*v02*log((r*r+(z/q)**2+rcmid**2)/rcmid**2))
c     +           *sqrt(4.*pi)
         apot(1,ir)=psi00*sqrt(4.*pi)*exp(-100*(r*r+z*z)/(nr*dr)**2)
         do l=2,lmaxx,2
            apot(l/2+1,ir)=0
         enddo
      enddo
      niter=(2+lmaxx/2)*10
      lmaxstep=2
      lmax=0
      ntheta=lmax*10+2
      ntheta=max(10,ntheta)

c  now iterate. number of iterations and change in lmax depends on initial conditions.
c  iterate on first cycle until tidal radius is stable, then once for every harmonic added,
c  and until convergence once l has reached the desired maximum lmax.
      
c  if only a disk potential
      if( ihaloflag .eq. 0 .and. ibulgeflag .eq. 0 ) lmax = lmaxx-2
      drtidal = 2*dr
      rtidalold = 1e10
      lmaxold=-2
      iteroutside = 0

      call pgbegin(0,'?',5,4)
      call pgask(.false.)
      call pgsch(2.)
      write(*,*)

      do iter=0,299
c         do i=0,nr
c            write(17,*) iter*(nr*dr)+i*dr,apot(1,i),adens(1,i)
c         enddo
c  determine number of harmonics, and number of azimuthal bins, to use this iteration
         if( lmax .eq. 0 .or. lmax .eq. lmaxx ) then
            if( drtidal .lt. dr .and. iter.gt.10) then
               lmax=lmax+lmaxstep
               ntheta=lmax*4+2
            endif
         else
            lmax=lmax+lmaxstep
            ntheta=lmax*4+2
         endif
         if( lmax .eq. lmaxx+lmaxstep ) goto 199

c  Now get the harmonics of the density in this potential --> adens
c  NB that dens only gives the density without an approximate sech**2
c  component --- that part of the potential is not represented by the 
c  harmonics.
c  The function dens(r,z) returns the density - high-frequency cpt
c  The full density is returned by totdens
         adens(1,0)=dens(0.,0.)*sqrt(4.*pi)
         do l=2,lmax,2
            adens(l/2+1,0)=0
         enddo
         do ir=1,nr
            adens(1,ir)=0
         enddo
c  nrmax will mark the outermost radial bin with non-zero density.
         nrmax=nr
         do l=0,lmax,2
c  integrate density * spherical harmonic function over quadrant
c  use cos(theta) as independent variable.
            do ir=1,nrmax
               r=ir*dr
               s=0
               dctheta=1.0/ntheta
               s=s+polardens(r,1.0,l)+polardens(r,0.0,l)
               do is=1,ntheta-1,2
                  ctheta=is*dctheta
                  s=s+4*polardens(r,ctheta,l)
               enddo
               do is=2,ntheta-2,2
                  ctheta=is*dctheta
                  s=s+2*polardens(r,ctheta,l)
               enddo
               s=s*dctheta/3.
               s=s*4*pi
               adens(l/2+1,ir)=s
c  mark the first even radial bin on which the density has fallen to zero.
               if (l.eq.0 .and. s.eq.0.) then
                  nrmax=nr
                  nrmax=nrmax+mod(nrmax,2)
                  goto 77
               endif
c     write(*,*) 'a(',l,') =',s,'at r=',r
            enddo
 77      enddo
c  now get the potential harmonics of this new density. (BT 2-208)
c  Simpson's rule integration. 
         do l=0,lmax,2
            s1(0)=0
            r = 2*dr
            s1(2)=(r*dr/3.)*(4*adens(l/2+1,1)*(1.0-dr/r)**(l+2)+
     +             adens(l/2+1,2))
            rold = r 
c  doesn't matter but should be nonzero
c  to avoid round-off error
            do ir=4,nr,2
               r=ir*dr
               s1a = (r*dr/3.)*(adens(l/2+1,ir-2)*(1.0-2*dr/r)**(l+2)+
     &           4*adens(l/2+1,ir-1)*(1.0-dr/r)**(l+2)+adens(l/2+1,ir))
               s1(ir) = s1a + s1(ir-2)*(rold/r)**(l+1)
               rold = r
            enddo
            s2(nr)=0
            rold = nr*dr
            do ir=nr-2,2,-2
               r=ir*dr
               s2a = (r*dr/3.)*(adens(l/2+1,ir+2)*(1.0+2*dr/r)**(1-l)+
     &          4*adens(l/2+1,ir+1)*(1.0+dr/r)**(1-l)+adens(l/2+1,ir))
               s2(ir) = s2a + s2(ir+2)*(r/rold)**l
               rold = r
            enddo
c  replace the potential harmonics with a mean of the previous iteration (25%) 
c  and the current one (75%). This damps out oscillations that otherwise occur.
c  if this is the first time this harmonic is calculated, use the entire new 
c  value.
            do ir=2,nr,2
               if (l.le.lmaxold) then
                 apot(l/2+1,ir)=0.25*apot(l/2+1,ir)-0.75*4*pi/(2.*l+1.)*
     +                 (s1(ir)+s2(ir))
               else
                  apot(l/2+1,ir)=-4*pi/(2.*l+1.)*
     +                 (s1(ir)+s2(ir))
               endif
            enddo
c  
c  Calculate the 1st and 2nd-order radial gradients
c  
            do ir=2,nr,2
               r = ir*dr
               fr(l/2+1,ir)=-4*pi/(2.*l+1.)*(-(l+1)*s1(ir) + l*s2(ir))/r
               fr2(l/2+1,ir)=-4*pi/(2.*l+1.)*
     +              ((l+1)*(l+2)*s1(ir)/r**2+ 
     +              l*(l-1)*s2(ir)/r**2 -(2*l+1)*adens(l/2+1,ir))
            enddo
         enddo
c  now interpolate the gaps 
c  first quadratically interpolate the monopole back to the origin.
c  the remaining multipoles are zero there.
         apot(1,0)=(4*apot(1,2)-apot(1,4))/3.
         fr(1,0)=0.0
         fr2(1,0)=2*fr2(1,2)-fr2(1,4)
         do l=2,lmax,2
            apot(l/2+1,0)=0
            fr(l/2+1,0)=0
            fr2(l/2+1,0)=0
         enddo
c  then linearly interpolate other bins.
         do ir=1,nr-1,2
            do l=0,lmax,2
               apot(l/2+1,ir)=(apot(l/2+1,ir-1)+apot(l/2+1,ir+1))/2.
               fr(l/2+1,ir)=(fr(l/2+1,ir-1)+fr(l/2+1,ir+1))/2.
               fr2(l/2+1,ir)=(fr2(l/2+1,ir-1)+fr2(l/2+1,ir+1))/2.
            enddo
         enddo
c plot current harmonic functions

c if you are only a disk potential then no need to iterate exit the loop
c
         if( ihaloflag .eq. 0 .and. ibulgeflag .eq. 0 ) goto 199
c  finally reset the potential at the origin to psi00
c  Note that the fake disk potential is zero at the origin.
         a00=apot(1,0)
         do ir=0,nr
            apot(1,ir)=apot(1,ir)+psi00*sqrt(4.*pi)-a00
         enddo
         if (a00/sqrt(4.*pi).gt.psi00) then
            write(*,'(''Iter'',i4,'': lmax='',i4,
     +           '', tidal radius is infinite'')') iter,lmax
            iteroutside = iteroutside + 1
            if( iteroutside .gt. 20 ) then
                write(*,'(''nr='',i4,'' is too small'',
     +           '', try larger number of radial bins  
     +          - exiting program'')') nr
                goto 12345
            endif
            drtidal = 2.0*dr
         else
            potr=psi00
            do ir=1,nr
c  look for new tidal radius of the model
c  defined as the radius at the equator where the potential is zero
               potrm1=potr
               potr=pot(ir*dr,0.,apot,dr,nr,lmax)
               if (potrm1*potr.le.0.) then
                  dpot = potr - potrm1
                  if( dpot .eq. 0.0 ) then
                     rtidal = (ir-1)*dr
                  else
                     rtidal=(ir-1-potrm1/dpot)*dr
                  endif
                  write(*,'(''Iter'',i4,'': lmax='',i4,
     +                 '', tidal radius is '',g15.6)') iter,lmax,rtidal
                  drtidal = abs(rtidal - rtidalold)
                  rtidalold = rtidal
                  goto 9
               endif
            enddo
            write(*,'(''Iter'',i4,'': lmax='',i4,
     +           '', tidal radius is outside grid'')') iter,lmax
            drtidal = 2.0*dr
 9       endif
c write out the changes in the potential at the reference points 
c at this iteration.
         if( idiskflag .eq. 1 ) then 
            do iref=1,10
               oldpref(iref)=pref(iref)
               pref(iref)=pot(rref(iref),zref(iref),apot,dr,nr,lmax)
               enddo
            write(18,'(2i3,10g12.4)') iter,lmax,
     +                  (pref(iref)-oldpref(iref),iref=1,10)
          endif
c  now repeat with this new potential!
         lmaxold=lmax
         call dbhplot(apot,lmax,nr,dr)
      enddo
c  end of the main loop
 199  continue
      call pgend
      if( idiskflag .eq. 1 ) close(18)
c  write final density, potential and potential gradient harmonics
c
      totalmass = fr(1,nr)/sqrt(4*pi)*(dr*nr)**2
      write(*,*) 'psi00=',psi00,'Rt/RK=',rtidal/rking
      write(*,*) 'Total mass=', totalmass

c  
c  Calculate force and potential for halo only
c  
      lmax = lmaxx
      halomass = 0.0
      bulgemass = 0.0
      haloedge = 0
      bulgeedge = 0
      if( ihaloflag .eq. 1 ) call halopotential(halomass, haloedge)
      if( ibulgeflag .eq. 1 ) call bulgepotential(bulgemass, bulgeedge)

      diskmass = totalmass - halomass - bulgemass
      diskedge = outdisk + 2.0*drtrunc

      open(20,file='mr.dat',status='unknown')
      write(20,*) diskmass, diskedge
      write(20,*) bulgemass, bulgeedge
      write(20,*) halomass, haloedge
      close(20)

c      redge = nr*dr
c      do l=0,lmax,2
c          c0 = apot(l/2+1,nr) + fr(l/2+1,nr)*redge/(l+1)
c          do i=0,nr
c              apot(l/2+1,i) = apot(l/2+1,i) - c0
c          enddo
c      enddo

      open(11,file='dbh.dat',status='unknown')
      write(11,'('' # psi00,v0,q,a,b,c,dr,nr,lmax='')')
      write(11,'('' #'',7g15.5,i6,i4)') psi00,v0,q,a,b,c,dr,nr,lmaxx
      write(11,'('' # bulge psic, rho1, sig:'')')
      write(11,'('' #'',3g15.5)') psiout,rho1,sigbulge
      write(11,'('' # Mdisk, rdisk, zdisk, outdisk, drtrunc'')')
      write(11,'('' #'',5g15.5)') rmdisk, rdisk, zdisk, outdisk, drtrunc
      write(11,'('' #'',3i5)') idiskflag, ibulgeflag, ihaloflag
      write(11,'('' #  OUTPUT FROM DBH8. TOTAL POTENTIAL.'')')

      do ir=0,nr
         write(11,'(8g16.8)') ir*dr,(adens(l/2+1,ir),l=0,lmaxx,2)
      enddo
      write(11,*)
      write(11,*)
      do ir=0,nr
         write(11,'(8g16.8)') ir*dr,(apot(l/2+1,ir),l=0,lmaxx,2)
      enddo
      write(11,*)
      write(11,*)
      do ir=0,nr
         write(11,'(8g16.8)') ir*dr,(fr(l/2+1,ir),l=0,lmaxx,2)
      enddo
      write(11,*)
      write(11,*)
      do ir=0,nr
         write(11,'(8g16.8)') ir*dr,(fr2(l/2+1,ir),l=0,lmaxx,2)
      enddo
      write(11,*)
      do ir=0,nr
         psi=pot(ir*dr,0.,apot,dr,nr,lmax)
         write(11,'(2g16.8)') ir*dr,diskdens(ir*dr,0.,psi)
      enddo
      close(11)

      write(*,*) 'Final model written to file ''dbh.dat''.'
12345 continue
      
      end
