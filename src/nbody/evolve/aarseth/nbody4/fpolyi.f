***
      SUBROUTINE FPOLYI(IBODY)
*
*
*       Force & first derivative on GRAPE.
*       ---------------------------------
*
      INCLUDE 'common4.h'
      INTEGER i,j,k,l,ibody,i1,iph0,j1,ks,ierr,nn,np,nnb
      REAL*8 f2dot(3),dtmax,firr(3),fd(3)
*
*
* Find the largest block time-step commensurate with current time. 
      dtmax = dtk(1)
      CALL dtchck(time,dtmax,dtk(40))
*
* Define basic variables for body #IBODY (and possible KS component).
      i = ibody
 1    t0(i) = time
      step(i) = dtmax
* Copy X0 and initialize local and global second derivative.
      do k = 1,3
         x0(k,i) = x(k,i)
         f2dot(k) = 0.d0
         d2(k,i) = 0.d0
      enddo
      if (i.eq.ifirst.and.iphase.gt.0.and.iphase.lt.8)then
         i = ifirst + 1
         if (iphase.ne.4)then
            goto 1
         endif
      endif
      i = ibody
*
* Send variables to GRAPE and form neighbour list (iphase=-3 done in MDOT).
      if (iphase.ne.-3)then
         CALL gpsend
         h2 = 4.0*cmsep2*rmin2
         h2 = MIN(h2,0.25d0*rscale**2)
         iph0 = iphase
* Obtain F & FDOT from dominant component for initial guess in NBLIST.
         if (i.eq.ifirst.and.iphase.gt.0.and.iphase.lt.8.and.
     &       iphase.ne.4)then
             ilist(1) = 1
             ilist(2) = ifirst + 1
             phi(i) = -vc**2
             CALL ffdot(i)
         end if
         iphase = 0
    5    CALL nblist(i,h2,nnb)
         if ((nnb.lt.10.AND.N.GT.100).OR.
     &       (nnb.LT.1 + SQRT(FLOAT(N)/100))) then
            h2 = 2.0*h2
            go to 5
         end if
         iphase = iph0
         ilist(1) = nnb
      end if
*
* Consider the different cases (some are treated the same way).
*
      if (iphase.eq.1)then
* Combine mass-weighted PHI and F for new c.m. particle. 
         zm12 = body(icomp) + body(jcomp)
         phi(i) = (body(icomp)*phi(icomp) + body(jcomp)*phi(jcomp))/zm12
         do k = 1,3
            f(k,i) = body(icomp)*f(k,icomp) + body(jcomp)*f(k,jcomp)
            f(k,i) = 2.0d0*f(k,i)/zm12
         enddo
* Estimate FDOT from the neighbour list.
         CALL ffdot3(i,1)
*
      else if(iphase.eq.6)then
* Obtain F & FDOT from dominant terms and improve PHI.
         phi(i) = -vc**2
         CALL ffdot3(i,2)
*
      else if(iphase.eq.2)then
* Add body #I and form dominant F, FDOT & PHI for components.
         nnb = nnb + 1
         ilist(nnb+1) = i
         ilist(1) = nnb
         phi(i) = -vc**2
         phi(i+1) = -vc**2
         CALL ffdot(i)
         CALL ffdot(i+1)
*
      else if (iphase.eq.-1.or.iphase.eq.0.or.iphase.ge.7)then
* Treat remaining cases as for two standard KS components.
         phi(i) = -vc**2
         CALL ffdot(i)
*
      else if (iphase.eq.-3)then
* Scale current values from MDOT, COAL & CMBODY.
         nnb = ilist(1)
         do k = 1,3
            f(k,i) = 2.0*f(k,i)
            fdot(k,i) = 6.0*fdot(k,i)
         enddo
* Send all variables to GRAPE on first call only.
         if (isend.lt.0)then
            CALL gpsend
         endif
*
      else if (iphase.eq.4)then
         phi(i) = -vc**2
         CALL ffdot(i)
*
      else
          CALL FPOLY1(I,I,0)
          WRITE (6,10) iphase, i, name(i)
   10     FORMAT (' FPOLYI    BAD CASE    IPH I NAM ',I4,2I7)
          GO TO 80
*         CALL gpfree
*         STOP
      endif
*
* Specify nominal membership of new KS to indicate perturbed motion.
      if(i.gt.n.and.iphase.gt.0)then
         i1 = 2*(i - n) - 1
         list(1,i1) = 1
         list(2,i1) = n
      endif
*
* Calculate F & FDOT on body #I using GRAPE.
      iter = 0
*     ierr = 0
*     eps2 = 0.d0
*     nn = ntot - ifirst + 1
*
*40   CALL g6_set_ti(gpid,t0(i))
*     CALL g6calc_firsthalf(gpid,nn,1,i,x(1,i),xdot(1,i),
*    &                      f(1,i),fdot(1,i),phi(i),eps2,gph2)
*     gperr = g6calc_lasthalf(gpid,nn,1,i,x(1,i),xdot(1,i),
*    &                        eps2,gph2,gpacc,gpjerk,gppot)
40    CALL g6acc(1,i,x(1,i),xdot(1,i),f(1,i),fdot(1,i),phi(i))
*
* Check for hardware error.
*     if(gperr.ne.0)then
*        WRITE(6,*)'GRAPE HARDWARE ERROR: RESET REQUIRED',time
*        WRITE(6,*)'g6calc_lasthalf returned ',gperr,' in fpolyi'
*        CALL gpwipe(gpid,time)
*        ierr = ierr + 1
*        if(ierr.ge.10)then
*           WRITE(6,*)' TOO MANY RESETS '
*           WRITE(6,*)' STOP AT TIME = ',time
*           CALL gpfree
*           STOP
*        endif
*        CALL gpsend
*        goto 40
*     endif
*
* Copy F, FDOT & PHI into COMMON variables and initialize FIRR & FD.
      do k = 1,3
*        f(k,i) = gpacc(k,1)
*        fdot(k,i) = gpjerk(k,1)
         firr(k) = 0.d0
         fd(k) = 0.d0
      enddo
*     phi(i) = gppot(1)
*
* Distinguish between c.m. body and singles or old KS components.
      if(i.gt.n)then
* Save perturber list and copy neighbour list (ignore redundancy).
         i1 = 2*(i - n) - 1
         np = list(1,i1)
         do l = 2,np+1
            jpert(l) = list(l,i1)
         enddo
*
* Form temporary perturber list for new c.m. polynomial (also COAL).
         if(iphase.ge.0)then
* Copy members from ILIST (body #I not a member with CALL NBLIST).
            do l = 2,nnb+1
               list(l,i1) = ilist(l)
            enddo
            list(1,i1) = nnb
         endif
*
* Add differential force & first derivative corrections to c.m. body.
         CALL dfcm2(i,i1,firr,fd)
*
* Check large derivatives on GRAPE or host.
         if(fd(1).gt.-1.0d+12.and.fd(1).lt.1.0d+12.and.
     &       fdot(1,i).gt.-1.0d+15.and.fdot(1,i).lt.1.0d+15) goto 62
         np1 = list(1,i1) + 1
         WRITE(6,61)i,iphase,ttot,r(i-n),(list(k,i1),k=1,np1)
 61      FORMAT(' WARNING!   FPOLYI   I IPH T R LST ',
     &                                i6,i4,f9.2,1p,e9.1,0p,4(2x,8i6,/))
* Obtain F, FDOT & PHI on the host instead (tested OK 21/5/99).
         CALL fpoly1(i,i,0)
         goto 80
*
* Restore perturber list (redundant for call from KSINIT & KSIN2).
 62      list(1,i1) = np
         do l = 2,np+1
             list(l,i1) = jpert(l)
         enddo
      else
* Correct GRAPE force & derivative on single body due to active KS.
         do l = 2,nnb+1
            j = ilist(l)
            if(j.gt.n)then
               j1 = 2*(j - n) - 1
               if(list(1,j1).gt.0)then
                  CALL dfsp2(i,j1,firr,fd)
               endif
            endif
         enddo
      endif
*
* Add differential c.m. corrections from neighbours.
      do k = 1,3
         f(k,i) = f(k,i) + firr(k)
         fdot(k,i) = fdot(k,i) + fd(k)
      enddo
*
* Check option for external force.
      if(kz(14).gt.0)then
         CALL xtrnld(i,i,1)
      endif
*
* Obtain new time-step, scale F & FDOT by factorials and update TNEXT.
      CALL steps(i,i,1)
*
* Send new F, FDOT & T0 of body #I to GRAPE.
80    CONTINUE
*80   gpaddr(1) = i - ifirst
*     gpindx(1) = i
*
*     CALL g6_set_j_particle(gpid,gpaddr,gpindx,t0(i),step(i),
*    &                       body(i),f2dot,fdot(1,i),f(1,i),
*    &                       x0dot(1,i),x0(1,i))
*
* See whether the second component should be treated (IPHASE > 0 but < 8).
      if(i.eq.ifirst.and.iphase.gt.0.and.iphase.lt.8)then
         i = ifirst + 1
         if (iphase.ne.4)then
            goto 40
         endif
      endif
*
* Set ISEND = 0 for skipping GPSEND in INTGRT.
      isend = 0
*
      RETURN
      END
***
