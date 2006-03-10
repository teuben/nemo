      SUBROUTINE TNEW
*
*
*       New astrophysical time scale.
*       -----------------------------
*
      INCLUDE 'common4.h'
      REAL*8 M1,TEV1,DLOGM,mrt,mrtfirst,trhstar
      integer itime,first
      real*8 tfirst
      data first/0/
*
*
*       Obtain logarithmic derivative dlog(m)/dlog(t) from CE relation.
      if (first.eq.0) then
         tfirst = time
         call dchesc(mrt,trhstar)
         mrtfirst = mrt
         tescstar = 1.0d+10
         myntry = 1
         first = 1
      else 
         if (time.ge.tfirst + myntry*dtadj) then
            call dchesc(mrt,trhstar)
c            WRITE (6,*) 'tnew: t*,M<rt ',time,mrt
            if (mrt.lt.mrtfirst) then
               tescstar = (time - tfirst)*mrtfirst/(mrtfirst - mrt)
               tesc = tescstar*dtstar
               mrtfirst = mrt
               tfirst = time
               myntry = 1
            else
               myntry = myntry + 1
            endif
         endif
      endif
      M1 = BODY1*ZMBAR
*       Check Chernoff-Weinberg updating.
      IF (KZ(19).GT.4) THEN
          tev1 = 0.d0
          CALL INTEV(M1,TEV1,DLOGM)
      END IF
*
*       Evaluate the characteristic time scale.
      IF (KZ(13).ne.0) THEN
         A2 = 2.d0 - ALPHA
         if (dlogm.gt.0.0) then
            TMl = (1.d0 - (BODYN/BODY1)**A2)*TPHYS/
     &                                 ((ALPHA - 2.d0)*DLOGM)
         else
            tml = 1.0d+11
         endif
      END IF
*
*       Restrict the time scale to the turnoff mass until first mass loss.
      IF (NMDOT.EQ.0) THEN
          TMl = 1.0d+06
      END IF
*
*       Update NGLOB and define reference time from relaxation time ratio.
      if (float(nzero).eq.0.0.or.float(n).eq.0.0.or.
     &     log(0.01*nglob).eq.0.0) then
         WRITE (6,*) 'tnew: float(nzero),float(n),log(0.01*nglob)',
     &        float(nzero),float(n),log(0.01*nglob)
         call gpfree
         stop
      endif
      nbin = (ifirst - 1)/2
cNGLOB0 has been corrected for the number of binaries in binpop
      NGLOB = NGLOB0*FLOAT(N - nbin)/FLOAT(NZERO - NBIN0)
      RATIO = FLOAT(NGLOB)/FLOAT(N - nbin)*LOG(1.d0+0.01d0*(N-nbin))
     &     /LOG(0.01d0*NGLOB)
      TREF = RATIO*TCR0*UT
*
*       Form differential time scale for the relevant case.
c      factor = 1.d0
c      factor = kz(13)
c      IF (tml.LT.factor*TCR0*UT.or.tescstar.lt.trhstar) THEN
      if (kz(13).eq.1) factor = 1.d0
      if (kz(13).eq.-2) then
         tml = 0.d0
         factor = 1.d0
      endif
      if (kz(13).eq.2) factor = (N-nbin)/
     &                          (11.d0*log(1.d0+0.01d0*(n-nbin)))
      if (kz(13).eq.-1) then
          DTSTAR = RATIO*UT
      else
         IF (tml.LT.factor*TCR0*UT) THEN
            DTSTAR = UT
         ELSE IF (tml.LT.TREF*factor) THEN
            DTSTAR = tml/(factor*TCR0)
         ELSE
            DTSTAR = RATIO*UT
         END IF
      endif
*
*       Update the astrophysical time and corresponding conversion factors.
c      WRITE (6,*) '1st tphys in tnew ',tphys
cCheck that no times for updating tstar have been missed.
      if (tstar.ne.0.0) then
         oldtime = tphys/tstar
c         i = (time - oldtime)/stepx + 0.5
c         if (i.ne.1) WRITE (37,*) 'i <> 1 in tnew - t,i:',t,i
      else
c         i = 1
         oldtime = 0.d0
      endif
c      
      TPHYS0 = TPHYS
      TPHYS = TPHYS + (time - oldtime)*DTSTAR
      TSTAR = TPHYS/TIME
      TSCALE = TSTAR
*
      RMSTAR = 1.d0/(nglob0/(nzero-nbin0))**(1.d0/3.d0)*(dtstar/ut)
*
*       Redetermine evolution time for first star if TPHYS > initial value.
      IF (NMDOT.EQ.0.AND.TPHYS.GT.TEV10) THEN
          TEV1 = 1.0d+10
          DO 10 I = 1,N
              IF (TEV(I).LT.TEV1) THEN
                  TEV1 = TEV(I)
                  IM = I
              END IF
   10     CONTINUE
*       Specify current mass loss time and update TMDOT for routine MDOT2.
          TEV(IM) = TIME
          TMDOT = TEV(IM) - STEPX
      END IF
*
*       Check updating of STEPX and TEV.
      IF (KZ(19).LE.4) THEN
          IF (STEPX*TSTAR.GT.2.0d-03) THEN
              STEPX = 0.5d0*STEPX
              WRITE (6,12)  TPHYS, TSTAR, STEPX, RMSTAR
   12         FORMAT (' TNEW:    TP T* STEPX RM* ',2F8.1,1P,3E10.2)
          END IF
          TFAC = TPHYS0/TPHYS
          TEV1 = 1.0d+10
          DO 15 I = 1,NTOT
              IF (TEV(I).LT.9.9E+09) THEN
                  TEV(I) = TFAC*TEV(I)
                  TEV0(I) = TFAC*TEV0(I)
                  TEV2 = TEV1
                  TEV1 = MIN(TEV(I),TEV2)
              END IF
   15     CONTINUE
          TMDOT = TEV1 - STEPX
      END IF
*
*     IF (KZ(13).GT.1) THEN
      itime = time
c      IF (TIME.LT.5*STEPX.OR.ABS(TIME-8.0).LT.2*STEPX.OR.
c     &    ABS(TIME-10.0).LT.2*STEPX.OR.ABS(TIME-6.0).LT.2*STEPX.OR.
c     7    ABS(TIME-4.0).LT.2*STEPX) THEN
      if (time - itime.lt.stepx) then
c      if (.true.) then
          WRITE (37,20)  TIME, M1, TMl, tesc, TREF, TPHYS, TSTAR,rmstar
   20     FORMAT (' TIME SCALE:    T M TMl tesc TREF TPH T* RM*',
     &                             F10.4,2F7.2,5F9.2)
c          WRITE (37,*) 'tcr0,ut,ratio,dtstar,stepx,tphys ',
c     &         tcr0,ut,ratio,dtstar,stepx,tphys
          CALL FLUSH(37)
      END IF
*
      RETURN
*
      END
c
c

      subroutine dchesc(mrt,trhstar)
      INCLUDE 'common4.h'
      REAL*8  R2,RHO,mrt,rt,mymrt,trhstar,mrth,rh,rlag(13),mfrac(13),
     &	mass,mc,rcore
      integer ncore
      COMMON/WORK1/ R2(NMAX),RHO(NMAX)
      REAL*8  C(3)
	data mfrac/0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.3,
     &     0.4,0.5,0.75,0.9/
*
*       Set square radii of single particles & c.m. bodies.
      NP = 0
      do i = 1,3
         c(i) = rdens(i)
      enddo
      DO 10 I = IFIRST,NTOT
          NP = NP + 1
          R2(NP) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2
     &                                + (X(3,I) - C(3))**2
          JLIST(NP) = I
   10 CONTINUE
*
*       Sort square distances with respect to the centre C.
      CALL SORT1(NP,R2,JLIST)
*
c Determine tidal radius and mass if kz(14) = 2
c
      if ((kz(14).eq.2.or.kz(14).eq.-1).and.tidal(1).ne.0.) then
c ZMASS may not be up to date:
         mrt = 0.d0
         do i = ifirst,ntot
            mrt = mrt + body(i)
         enddo
      I = np
c      mrt = zmass
   30 continue
      if (mrt.gt.0.) then
         rt = (mrt/tidal(1))**(1.d0/3.d0)
      else
         rt = 0.d0
      endif
      IF (sqrt(r2(i)).gt.rt.and.rt.gt.0.and.i.gt.1) then
         IM = JLIST(I)
         Mrt = Mrt - BODY(IM)
         i = i - 1
         GO TO 30
      endif
      myn = i
      mymrt = mrt
*
c	Now find Lagrangian radii within rt
c
	ilag = 1
	i = 1
        fmax = 0.d0
	im = jlist(i)
	mass = body(im)
35	continue
        ftemp = mass/r2(i)
        if (ftemp.gt.fmax) then
           fmax = ftemp
           mc = mass
           rcore = sqrt(r2(i))
           ncore = i
        endif
	if (mass.gt.mfrac(ilag)*mrt) then
		rlag(ilag) = sqrt(r2(i))
		ilag = ilag + 1
		if (ilag.le.13) go to 35
	endif
	if (ilag.le.13) then
		i = i + 1
		im = jlist(i)
		mass = mass + body(im)
		goto 35
	endif
	mrth = mrt
   40 continue
      IF (i.ge.1.and.mrth.gt.0.5*mrt) then
         IM = JLIST(I)
         Mrth = Mrth - BODY(IM)
         i = i - 1
         GO TO 40
      endif
	rh = sqrt(r2(i))
	if (mrt.gt.0.0) then
		trhstar = 0.138d0*myn*rh**1.5/
     &                    (mrt**0.5*log(1.d0+0.01d0*myn))
	else
		trhstar = 0.d0
	endif
      WRITE (6,'(''tnew:''f10.2,f10.5,f10.3,i6,f10.3,2f10.5,i6)')  
     &       time*tstar,mrt,rt,myn,trhstar,mc,rcore,ncore
	WRITE (6,'(''rlag: '',13f9.4)') (rlag(i),i=1,13)
        call flush(3)
c
cNow find the mass function. Save Sverre's rtide
c
        sjartide = rtide
        rtide = rt
        call global
c
cRestore Sverre's rtide
c
        rtide = sjartide
        if (myn.le.10) then
           WRITE (6,*) 'n(<rt) < 10; stopping'
           call gpfree
           stop
        endif
	else
		WRITE (6,*) 'dchesc entered with bad options'

      endif
      return
      end
