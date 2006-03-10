***
      SUBROUTINE INTGRT
*
*
*       N-body integrator flow control.
*       -------------------------------
*
      INCLUDE 'common4.h'
*
      INTEGER i,j,k,l,ii,kk,i10,j1,j2,kp1,kp2
      INTEGER nxtlen,next,ncpert,nlsum,ierr
      INTEGER nxtlst(nmax),icpert(kmax),icblk(kmax),listq(nmax),nl(20)
      INTEGER nn,np,np0,jnext,kcorr,jp,jpair,lp,ls
      INTEGER io,ireg,iq,icall,lq,isave,jsave,gpindx(48)
*
      REAL*8 xgp(3,48),vgp(3,48),gpacc(3,48),gpjerk(3,48),gppot(48)
*     REAL*8 fgp(3,nmax),f1gp(3,nmax),phgp(nmax)
      REAL*8 a(3,200),adot(3,200),potz(200)
      REAL*8 eps2,stepm,tcomp0,tcomp,share,ration
*
      DATA iq,icall,lq /0,0,12/
      DATA eps2,stepm,tcomp0 /0.d0,0.03125d0,0.d0/
      SAVE iq,icall,lq,nq,eps2,stepm,tcomp0
      PARAMETER (ration=5.d0)
      COMMON /CLOUDS/ xcl(3,mcl),xdotcl(3,mcl),bodycl(mcl),rcl2(mcl),
     &                clm(mcl),clmdot(mcl),vcl,sigma,rb2,rb3,pcl2,
     &                tcl,stepcl,tbig,dtbig,ncl,newcl
      EXTERNAL short
      SAVE ITER
      DATA ITER /0/

*
*       Initialize GRAPE and obtain current number of pipes.
 1    if(iphase.le.-2.and.gpstat.eq.0)then
*        CALL gpinit(gpid)
         if (NZERO.GE.8000) stepm = 1.0D0/64.0D0
         if (NZERO.GE.64000) stepm = 1.0D0/128.0D0
*        gpstat = 1
*        do k = 1,48
*           gph2(i) = 0.d0
*        enddo
      endif
*     npipe = g6_npipes()
      npipe = 48
*
*       Search for high velocities after escape or KS/chain termination.
 5    if(kz(37).gt.0.and.(iphase.lt.0.or.iphase.ge.2))then
          CALL hivel(0)
      endif
*
*       Reset regularization index and total block length on each return.
 6    ireg = 0
      if(iq.lt.0) icall = 0
      iphase = 0
      iq = 0
*       Enforce new block step search on significant changes.
      tlistq = time
*
*       Check whether to send all single particles and c.m. to GRAPE.
      if(isend.ne.0)then
         CALL gpsend
*        nn = ntot - ifirst + 1
*        CALL g6_setup_njdata(gpid,nn)
      endif
*
*       Form new list of active c.m. with zero mass on GRAPE.
      ncpert = 0
      j1 = -1
      do jpair = 1,npairs
         j1 = j1 + 2
         if(list(1,j1).gt.0)then
            ncpert = ncpert + 1
            icpert(ncpert) = n + jpair
         endif
      enddo
*
* Find next block to be advanced and set new time (restore nxtlen).
 10   icall = icall + 1
*     nn = ntot - ifirst + 1
* Reset TMIN second & third time after change to catch new chain step.
      if(time.ge.tlistq.or.icall.le.3)then
* Update interval by optimization at major times (square root of N).
         if(DMOD(tlistq,1.d0).eq.0.d0)then
            do l = 1,20
               nl(l) = 0
            enddo
            do i = ifirst,ntot
* Count steps at different levels for small values.
               do l = 10,16
                  if(step(i).lt.dtk(l)) nl(l) = nl(l) + 1
               enddo
            enddo
            nlsum = 0
* Determine interval by summing smallest steps until near sqrt(N).
            nsq = SQRT(float(n - npairs))
            do l = 16,10,-1
               nlsum = nlsum + nl(l)
               if(nlsum.le.nsq) lq = l
            enddo
         endif
*
* Increase interval by optimized value.
         nq = 0
         tmin = 1.0d+10
 18      tlistq = tlistq + dtk(lq)
         do i = ifirst,ntot
            if(tnext(i).le.tlistq)then
               nq = nq + 1
               listq(nq) = i
               tmin = MIN(tnext(i),tmin)
            endif
         enddo
* Increase interval in rare case of zero membership.
         if(nq.eq.0) goto 18
      endif
*
* Select members on new time-step block (t0 + step = tmin).
      CALL inext(nq,listq,tmin,nxtlen,nxtlst)
      i = nxtlst(1)
      time = t0(i) + step(i)
*
*     IF (TIME.GT.196.7921) THEN
*     WRITE (6,17) I, NXTLEN, NSTEPU, NSTEPI, TIME, STEP(I)
*  17 FORMAT (' INT   I LEN # T DT ',I6,I4,2I9,F12.6,1P,E10.2)
*     CALL FLUSH(6)
*     END IF
*
* Re-determine list if current time exceeds the boundary.
      if(time.gt.tlistq) goto 10
*
* Check new KS & output, advance KS/chain and predict at end of block.
      if(tprev.ne.time)then
* Save current block time (used for termination of regularizations).
         tblock = time
         if(ireg.gt.0)then
            time = tprev
            ICOMP = ISAVE
            JCOMP = JSAVE
            iphase = 1
            goto 200
         endif
* Check next adjust time at the end of each integration cycle.
         if(time.gt.tadj)then
            time = tadj
            iphase = 3
            goto 200
         endif
*
         nblock = nblock + 1
*
* See whether to advance any close encounters at first new time.
         if(time.gt.tprev)then
            CALL subint(iq,ncpert,icpert,i10)
            if(iq.lt.0) goto 5
         endif
*
* Make full N prediction to emulate GRAPE hardware.
      DO 20 J = IFIRST,NTOT
          S = TIME - T0(J)
          S1 = 1.5*S
          S2 = 2.0*S
          X(1,J) = ((FDOT(1,J)*S + F(1,J))*S + X0DOT(1,J))*S + X0(1,J)
          X(2,J) = ((FDOT(2,J)*S + F(2,J))*S + X0DOT(2,J))*S + X0(2,J)
          X(3,J) = ((FDOT(3,J)*S + F(3,J))*S + X0DOT(3,J))*S + X0(3,J)
          XDOT(1,J) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
          XDOT(2,J) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
          XDOT(3,J) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
   20 CONTINUE
*
* Copy the first block of particles to be integrated. 
         next = MIN(nxtlen,npipe)
         do ii = 1,next
            j = nxtlst(ii)
            s = time - t0(j)
            s1 = 1.5d0*s
            s2 = 2.d0*s
            s3 = s1*s2*s
            s4 = 0.75d0*s**4
* Improve prediction to second force derivative (for GRAPE emulator only).
            do kk = 1,3
*              x(kk,j) = ((fdot(kk,j)*s + f(kk,j))*s + x0dot(kk,j))*s
*    &                                               + x0(kk,j)
*              xdot(kk,j) = (fdot(kk,j)*s1 + f(kk,j))*s2 + x0dot(kk,j)
               xgp(kk,ii) = d2(kk,j)*s4 + x(kk,j)
               vgp(kk,ii) = d2(kk,j)*s3 + xdot(kk,j)
            enddo
         enddo
      endif
*
* Define indices for updating next and correcting previous block.
 40   kp1 = 1
      kp2 = next
      kcorr = 0
      nspert = 0
      jnext = 0
      np = 0
*     CALL g6_set_ti(gpid,time)
      tmin = 1.0d+10
*
* Set pointers and copy all block values of F, FDOT & PHI to GRAPE variables.
*     do k = 1,nxtlen
*        j = nxtlst(k)
*        gpaddr(k) = j - ifirst
*        gpindx(k) = j
*        do kk = 1,3
*           fgp(kk,k) = 2.d0*f(kk,j)
*           f1gp(kk,k) = 6.d0*fdot(kk,j)
*        enddo
*        phgp(k) = phi(j)
*     enddo
*
* Loop over new block (calculate first, predict next & correct previous).
      do ii = 1,nxtlen,npipe
*
*        nn = n - npairs
         np0 = np
         np = nxtlen - ii + 1
         np = MIN(np,npipe)
* Copy address of all block members for the emulator.
         do j = 1,np
            jj = ii + j - 1
            gpindx(j) = nxtlst(jj)
         end do
*
         CALL g6acc(np,gpindx,xgp(1,1),vgp(1,1),
     &                        gpacc(1,1),gpjerk(1,1),gppot(1))
*        CALL g6calc_firsthalf(gpid,nn,np,gpindx(ii),
*    &                         xgp(1,ii),vgp(1,ii),
*    &                         fgp(1,ii),f1gp(1,ii),phgp(ii),
*    &                         eps2,gph2)
*
* Predict next block (if any) while GRAPE is busy (order FDOT or D2).
         if(kp2.ge.nxtlen) goto 70
         kp1 = kp2 + 1
         kp2 = kp2 + MIN(nxtlen-kp2,npipe)
         l = 0
         do k = kp1,kp2
            l = l + 1
            j = nxtlst(k)
            s = time - t0(j)
            s1 = 1.5d0*s
            s2 = 2.d0*s
            s3 = s1*s2*s
            s4 = 0.75d0*s**4
            do kk = 1,3
*              x(kk,j) = ((fdot(kk,j)*s + f(kk,j))*s + x0dot(kk,j))*s
*    &                                               + x0(kk,j)
*              xdot(kk,j) = (fdot(kk,j)*s1 + f(kk,j))*s2 + x0dot(kk,j)
               xgp(kk,l) = d2(kk,j)*s4 + x(kk,j)
               vgp(kk,l) = d2(kk,j)*s3 + xdot(kk,j)
            enddo
         enddo
*
* Predict X & XDOT of active c.m. & KS components during first loop.
 70      if(kcorr.eq.0)then
            do l = 1,ncpert
               j = icpert(l)
               CALL xvpred(j,-2)
               zz = 1.d0
* Distinguish between low and high-order prediction of U & UDOT.
               if(gamma(j-n).gt.1.0d-04) zz = 0.d0
               CALL ksres2(j-n,j1,j2,zz)
            enddo
*
* See whether short time-step list needs updating (skip block > 32).
            if(nxtlen.le.32)then
               CALL short(nxtlen,nxtlst)
            endif
            kcorr = 1
         else
* Correct previous block and set new steps.
            do k = 1,np0
               jnext = jnext + 1
               i = nxtlst(jnext)
               phi(i) = potz(k)
               CALL nbcorr(i,ireg,ncpert,icpert,a(1,k),adot(1,k))
               tmin = MIN(tnext(i),tmin)
            enddo
         endif
*
*        gperr = g6calc_lasthalf(gpid,nn,np,gpindx(ii),
*    &                           xgp(1,ii),vgp(1,ii),
*    &                           eps2,gph2,gpacc,gpjerk,gppot)
*
* Copy results from GRAPE.
         do k = 1,np
            do kk = 1,3
               a(kk,k) = gpacc(kk,k)
               adot(kk,k) = gpjerk(kk,k)
            enddo
            potz(k) = gppot(k)
         enddo
      enddo
*
* Correct last particles in block and set new steps.
      do k = 1,np
         jnext = jnext + 1
         i = nxtlst(jnext)
         phi(i) = potz(k)
         CALL nbcorr(i,ireg,ncpert,icpert,a(1,k),adot(1,k))
         tmin = MIN(tnext(i),tmin)
      enddo
*
* Send corrected X & XDOT and prediction variables & mass to GRAPE.
*     gpt0 = time
      do ii = 1,nxtlen
         i = nxtlst(ii)
*        gpaddr(ii) = i - ifirst
*        gpindx(ii) = i
* Copy all corrected coordinates & velocities (only at the end of cycle).
         do kk = 1,3
            x(kk,i) = x0(kk,i)
            xdot(kk,i) = x0dot(kk,i)
         enddo
*        CALL g6_set_j_particle(gpid,gpaddr(ii),gpindx(ii),gpt0,step(i),
*    &                          body(i),d2(1,i),fdot(1,i),f(1,i),
*    &                          x0dot(1,i),x0(1,i))
      enddo
*
* Update time of current block and copy any new KS members at end of cycle.
      tprev = time
      if (ireg.gt.0) then
          ISAVE = ICOMP
          JSAVE = JCOMP
      end if
*
* Update integration of any tidal tail members.
      IF (NTAIL.GT.0) THEN
* Allow large quantized intervals with internal sub-integration.
          IF (DMOD(TIME,0.25D0).EQ.0.0D0) THEN
              DO 80 I = ITAIL0,NTTOT
                  IF (TNEXT(I).LE.TIME) THEN
                      CALL NTINT(I)
                  END IF
 80           CONTINUE
          END IF
      END IF
*
* Check optional histogram for active pipes (excluding KS).
      if(kz(33).gt.1)then
         CALL pipes(nxtlen)
      endif
*
* Exit on KS/merger termination, new multiple regularization or merger.
      if(iq.ne.0)then
         nbprev = 0
         if(iq.ge.4.and.iq.ne.7)then
            CALL delay(iq,-1)
         else
* Ensure correct KS index (KSPAIR may denote second termination).
            kspair = kvec(i10)
            iphase = iq
         endif
         goto 200
      endif
*
* Perform optional high-velocity check on time-steps at major times.
      if(kz(37).gt.0.and.nhi.gt.0)then
         if(DMOD(time,stepm).eq.0.d0)then
            CALL shrink
            if(nhi.gt.0)then
               CALL hivel(-1)
            endif
         endif
      endif
*
* Include optional movie at commensurate times of STEPY.
*     if(kz(40).gt.1)then
*        if(DMOD(time,stepy).eq.0.d0)then
*           CALL movie
*        endif
*     endif
*
* Check integration of interstellar clouds.
      if(kz(12).lt.0)then
         if(DMOD(time,stepcl).eq.0.d0)then
            CALL clint
         endif
      endif
*
* Include optional integration of cluster guiding centre.
      if(kz(14).gt.2)then
         if(DMOD(time,stepx).eq.0.d0)then
            CALL gcint
         endif
      endif
*
* Check optional disk shocking and mass loss time at end of block.
      if(kz(19).ne.0)then
* Delay until time commensurate with 1000-year step (new polynomials).
         if(DMOD(time,stepx).eq.0.d0)then
*
* Include optional updating of the astrophysical time scale.
            if(kz(13).ne.0)then
               CALL tnew
            endif
*
* Check astrophysical time for optional disk shock (next integer TIME).
            if(kz(12).gt.0.and.DMOD(time,1.d0).eq.0.d0)then
               if((time+toff)*tscale.gt.tshock)then
                  CALL shock
                  iq = -1
               endif
            endif
*
* Check next time for mass loss.
            if(time.gt.tmdot)then
               if(kz(19).gt.4)then
                  CALL mdot2
               elseif(kz(19).ge.3)then
                  CALL mdot
               elseif(kz(19).gt.0)then
                  CALL mloss
               endif
            endif
*
* Ensure full sorting after significant changes (shock or new NPAIRS).
            if(iq.lt.0.or.iphase.lt.0)then
               goto 6
            endif
         endif
      else
         if(kz(12).gt.1.and.DMOD(time,1.d0).eq.0.d0)then
            if((time+toff)*tscale.gt.tshock)then
               CALL shock
               goto 6
            endif
         endif
      endif
*
* Advance counters and check timer & optional COMMON save (NSUB = 0).
      ntimer = ntimer + nxtlen
      nsteps = nsteps + nxtlen
*
* Check optional time sharing (release GRAPE if elapsed CPU > RATION).
*     if(kz(38).gt.0.and.ntimer.gt.nmax)then
*        CALL cputim(tcomp)
*        share = tcomp - tcomp0
*        if(share.gt.ration)then
*           tcomp0 = tcomp
*           goto 165
*        else
*           goto 10
*        endif 
*     endif
*
      if(n.gt.5000.and.ntimer.lt.10*nmax) goto 10
      if(n.le.5000.and.ntimer.lt.nmax) goto 10
 165  ntimer = 0
*
* Check optional safety dump on unit #1.
      if(nsteps.ge.250*nmax.and.nsub.eq.0)then
         nsteps = 0
         if(kz(1).gt.1) CALL mydump(1,1)
      endif
*
* Include facility for termination of run (create dummy file STOP).
      OPEN(99,file='STOP',status='old',form='formatted',iostat=io)
      if(io.eq.0)then
         CLOSE(99)
         if(nsub.eq.0) WRITE(6,170)
 170     FORMAT(/,9x,'termination by manual intervention')
         cpu = 0.d0
      endif
*
* See whether GRAPE should be released for another user (then wait).
*     if(kz(38).gt.0)then
*        CALL gpfree
*        CALL mysleep(1)
*        gpstat = 0
*        iphase = -2
*        isend = -1
*        if(cpu.gt.0.d0) goto 1
*     endif
*
* Repeat cycle until elapsed computing time exceeds the limit.
      CALL cputim(tcomp)
      if(tcomp.lt.cpu) goto 10
*
* Do not terminate during triple, quad or chain regularization.
      if(nsub.gt.0) goto 10
*
* Terminate run with optional COMMON save.
      if(kz(1).gt.0)then
         cputot = cputot + tcomp - cpu0
         CALL mydump(1,1)
         WRITE(6,190)time+toff,tcomp,cputot/60.0,errtot,detot
  190    FORMAT(//,9x,'COMMON SAVED AT TIME =',f8.2,'  TCOMP =',f7.1,
     &                '  CPUTOT =',f6.1,'  ERRTOT =',f10.6,
     &                '  DETOT =',f10.6)
      endif
*
* Liberate GRAPE for next run before stopping (unless done above).
*     if(kz(38).eq.0.or.cpu.gt.0.d0)then
*        CALL gpfree
*     endif
      STOP
*
  200 TTOT = TIME + TOFF
      TPHYS = TTOT*TSTAR
*
      RETURN
      END
***
