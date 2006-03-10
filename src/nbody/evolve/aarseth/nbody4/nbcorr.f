***
      SUBROUTINE NBCORR(I,IKS,NCPERT,ICPERT,FIRR,FD)
*
*
*       N-body corrector.
*       -----------------
*
      INCLUDE 'common4.h'
*     INCLUDE 'grape6.h'
      INTEGER icpert(kmax)
      REAL*8 firr(3),fd(3),f2dot(3),xi(3),xidot(3)
      REAL*8 vi2,ri2,resc2,one18,frx(3),fdx(3)
      PARAMETER (one18=1.d0/18.d0)
      COMMON /EXTRA/ rpert2,i1,i2
      COMMON /CHAINC/ xc(3,ncmax),uc(3,ncmax),bodyc(ncmax),ich,
     &                listc(lmax),listcm(lmax)
*
* Check regularization criterion for single particles.
      if(step(i).lt.dtmin.and.i.le.n)then
* See whether dominant body can be regularized.
         if(iks.eq.0)then
            CALL search(i,iks)
         endif
         if(iks.gt.0.and.npairs.eq.kmax)then
            WRITE(6,1)i,jcomp,step(i)
            iks = 0
         endif
      endif
 1    FORMAT(5x,'WARNING!    KS LIMIT    I JCOMP DT ',2i6,1p,e9.1)
*
* Include close encounter search for low-eccentric massive binaries.
      if(iks.eq.0.and.step(i).lt.4.d0*dtmin)then
* Consider massive single bodies in absence of subsystems.
         if(i.le.n.and.body(i).gt.2.d0*bodym.and.nsub.eq.0)then
* Obtain two-body elements and relative perturbation.
            if(lshort(1).gt.2)then
               jmin = 0
               CALL orbit(i,jmin,semi,ecc,gi)
               if(jmin.ge.ifirst.and.jmin.le.n)then
                  eb = -0.5d0*body(i)*body(jmin)/semi
                  if(eb.lt.ebh.and.gi.lt.0.25d0)then
                     apo = semi*(1.d0 + ecc)
* Check eccentricity (cf. max perturbation) and neighbour radius.
                     if(ecc.lt.0.5d0)then
                        iks = iks + 1
                        icomp = i
                        jcomp = jmin
                     endif
                  endif
               endif
            endif
         endif
      endif
*
* Correct F & FDOT due to perturbers on c.m. or active KS on I <= N.
      if(npairs.gt.0)then
* Distinguish between c.m. and single particles.
         if(i.gt.n)then
            i1 = 2*(i-n) - 1
* Include force & first derivative corrections on perturbed pair.
            if(list(1,i1).gt.0)then
               CALL dfcm2(i,i1,firr,fd)
            else
* Search the perturber list of each active KS pair for body #I.
               do l = 1,ncpert
                  j1 = 2*(icpert(l)-n) - 1
                  np1 = list(1,j1) + 1
                  do k = np1,2,-1
                     j = list(k,j1)
                     if(j.lt.i)then
                        goto 10
                     elseif(j.eq.i)then
* Add force & first derivative corrections due to KS components.
                        CALL dfsp2(i,j1,firr,fd)
                        goto 10
                     endif
                  enddo
 10               continue
               enddo
            endif
         else
* Search the perturber list of each active KS pair for single body #I.
            do l = 1,ncpert
               j1 = 2*(icpert(l)-n) - 1
               np1 = list(1,j1) + 1
               do k = 2,np1
                  j = list(k,j1)
                  if(j.gt.i)then
                     goto 30
                  elseif(j.eq.i)then
* Add force & first derivative corrections due to KS components.
                     CALL dfsp2(i,j1,firr,fd)
                     goto 30
                  endif
               enddo
 30            continue
            enddo
         endif
      endif
*
* Include force correction due to regularized chain.
      if(nch.gt.0)then
* Distinguish between chain c.m. and any other particle.
         if(name(i).eq.0)then
            do k = 1,3
               xi(k) = x(k,i)
               xidot(k) = xdot(k,i)
            enddo
            CALL chf(i,xi,xidot,firr,fd)
         else
            nnb1 = listc(1) + 1
* Search the chain perturber list for #I.
            do l = 2,nnb1
               j = listc(l)
               if(j.gt.i) goto 50
               if(j.eq.i)then
                  do k = 1,3
                     xi(k) = x(k,i)
                     xidot(k) = xdot(k,i)
                  enddo
* Treat perturbed KS more carefully via routine kspert.f.
                  if (i.le.n)then
                     CALL fchain(i,xi,xidot,firr,fd)
                  else if(list(1,i1).gt.0)then
                     CALL kcpert(i,i1,firr,fd)
                  else
                     CALL fchain(i,xi,xidot,firr,fd)
                  endif
               endif
            enddo
         endif
      endif 
*
* Check option for external tidal field.
 50   if(kz(14).eq.1.or.kz(14).eq.2)then
         CALL xtrnlf(i,firr,fd)
      elseif(kz(14).eq.3)then
         do 55 k = 1,3
            frx(k) = firr(k)
            fdx(k) = fd(k)
 55      continue
         CALL xtrnlf(i,firr,fd)
         WDOT = 0.0
         W2DOT = 0.0
         W3DOT = 0.0
         DO 56 K = 1,3
             PX = FIRR(K) - FRX(K)
             DPX = FD(K) - FDX(K)
             WDOT = WDOT + XDOT(K,I)*PX
             W2DOT = W2DOT + FIRR(K)*PX + XDOT(K,I)*DPX
             W3DOT = W3DOT + 2.0*FIRR(K)*DPX + FD(K)*PX
 56      CONTINUE
         DT = TIME - T0(I)
*       Note: Taylor series at end of interval with negative argument.
         ETIDE = ETIDE - BODY(I)*((ONE6*W3DOT*DT - 0.5*W2DOT)*DT +
     &                                                 WDOT)*DT
      endif
*
* Check optional contribution from interstellar clouds.
      if(kz(12).lt.0)then
         CALL fcloud(i,firr,fd)
      endif
*
* Display formal procedure for time-step convergence test.
*     DO 60 K = 1,3
*         DVP(K) = XDOT(K,I) - X0DOT(K,I)
*         DVC(K) = X0DOT(K,I)
*  60 CONTINUE
*
* Include the corrector and set new T0, F, FDOT, D2 & D3.
      dt = time - t0(i)
      dtsq = dt**2
      dt6 = 6.d0/(dt*dtsq)
      dt2 = 2.d0/dtsq
      dtsq12 = one12*dtsq
      dt13 = one3*dt
      t0(i) = time
*
      do k = 1,3
	 df = 2.d0*f(k,i) - firr(k)
	 fd6 = 6.d0*fdot(k,i)
	 sum = fd6 + fd(k)
	 at3 = 2.d0*df + dt*sum
	 bt2 = -3.d0*df - dt*(sum + fd6)
	 x0(k,i) = x(k,i) + (0.6d0*at3 + bt2)*dtsq12
	 x0dot(k,i) = xdot(k,i) + (0.75d0*at3 + bt2)*dt13
	 f(k,i) = 0.5d0*firr(k)
	 fdot(k,i) = one6*fd(k)
         f2dot(k) = (3.d0*at3 + bt2)*dt2
	 d2(k,i) = one18*f2dot(k)
	 d3(k,i) = at3*dt6
* NOTE: D3 is a real derivative but F, FDOT & D2 are scaled for GRAPE!
      enddo
*
* Specify new time-step (use original expression instead of routine STEPI).
      ttmp = TSTEP(firr,fd,f2dot,d3(1,i),eta)
*     ttmp = STEPI(firr,fd,f2dot,d3(1,i),eta)
*
* Check for step reduction of hierarchical configurations.
      if(i.gt.n)then
         if(h(i-n).lt.-eclose.and.kz(36).gt.0)then
            if(gamma(i-n).gt.1.0d-04)then
               CALL kepler(i,ttmp)
            endif
         endif
      endif
*
*       Check convergence for large steps (cf. Makino, Ap.J. 369, 200).
      IF (TTMP.GT.STEPJ) THEN
         DV = 0.0
         FI = 0.0
         DO 70 K = 1,3
*           DVC(K) = X0DOT(K,I) - DVC(K)
*           DV = DV + (DVP(K) - DVC(K))**2
            DV = DV + (XDOT(K,I) - X0DOT(K,I))**2
            FI = FI + FIRR(K)**2
   70    CONTINUE
*       Employ Jun's criterion to avoid over-shooting (cf. Book, 2.16).
         DTJ = STEP(I)*(1.0D-06*STEP(I)**2*FI/DV)**0.1
         TTMP = MIN(TTMP,DTJ)
      END IF
      DT0 = TTMP
*
* Select discrete value (increased by 2, decreased by 2 or unchanged).
      if(ttmp.gt.2.d0*step(i))then
         if(DMOD(time,2.d0*step(i)).eq.0.d0)then 
            ttmp = MIN(2.d0*step(i),1.d0)
* Check for another factor 2 increase (large DT0 at major TIME values).
            if(dt0.gt.2.d0*ttmp)then
               if(DMOD(time,2.d0*ttmp).eq.0.d0)then
                  ttmp = MIN(2.d0*ttmp,1.d0)
               endif
            endif
         else
            ttmp = step(i) 
         endif
      elseif(ttmp.lt.step(i))then
         ttmp = 0.5d0*step(i)
* Check further decrease by 2 (new closest body switch to approach).
         if(ttmp.gt.dt0)then
            ttmp = 0.5d0*ttmp
         endif
         ttmp = MAX(ttmp,dtk(40))
      else
         ttmp = step(i)
      endif
*
* Try to increase step for high velocity escaper. 
*     if(vi2.gt.100.d0)then
*        ri2 = (x0(1,i) - rdens(1))**2 + (x0(2,i) - rdens(2))**2 +
*    &                                   (x0(3,i) - rdens(3))**2
*        resc2 = 4.d0*rtide**2
*        if(ri2.gt.resc2.and.ri2.lt.1.0d+10)then
*           ttmp = dtk(1)
*           CALL dtchck(time,ttmp,dtk(40))
*           write(6,*)' escaper step increased in nbcorr ',i,ttmp
*        endif
*     endif
*
* Set new block step and update next time.
      step(i) = ttmp
      tnext(i) = step(i) + t0(i)
*
* See whether any KS candidates are in the same block as body #I.
      if(iks.gt.0.and.i.eq.icomp)then
* Accept same time, otherwise reduce STEP(ICOMP) and/or delay.
         if(t0(jcomp).eq.t0(icomp))then
            icomp = MIN(icomp,jcomp)
            jcomp = MAX(i,jcomp)
         elseif(t0(jcomp)+step(jcomp).lt.
     &             t0(icomp)+step(icomp))then
            step(icomp) = 0.5d0*step(icomp)
            tnext(icomp) = t0(icomp) + step(icomp)
            iks = 0
         else
            iks = 0
         endif
      endif
*
* Check for boundary reflection.
*     IF (KZ(29).GT.0) THEN
*         RI2 = X(1,I)**2 + X(2,I)**2 + X(3,I)**2
*         IF (RI2.GT.RSPH2) THEN
*             CALL REFLCT(I,RI2)
*         END IF
*     END IF
*
      nstepi = nstepi + 1
*
      RETURN
      END
***
