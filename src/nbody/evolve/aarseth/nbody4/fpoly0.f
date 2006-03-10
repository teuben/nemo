***
      SUBROUTINE FPOLY0
*
*
*       Force & first derivative on GRAPE.
*       -----------------------------------
*
      INCLUDE 'common4.h'
      INTEGER i,ii,k,nn,np,ip,jj,iloop,gpindx(48)
      REAL*8 f2dot(3),a(9),h2,gpacc(3,48),gpjerk(3,48),gppot(48),
     &       xgp(3,48),vgp(3,48)
*
*
* Initialize GRAPE and obtain current number of pipes.
*     if(gpstat.eq.0)then
*        CALL gpinit(gpid)
*        gpstat = 1
*        do k = 1,48
*           gph2(k) = 0.d0
*        enddo
*     endif
*     npipe = g6_npipes()
      npipe = 48
      nn = n - npairs
*
* Send all single particles to GRAPE with no prediction (t0 = time).
*     do k = 1,3
*        f2dot(k) = 0.d0
*     enddo
*     rsoft = 0.1
*     do i = 1,48
*        gpindx(i) = ifirst
*     enddo
      do i = ifirst,ntot
*        gpaddr(i) = i - ifirst
*        gpindx(i) = i
         t0(i) = time
         step(i) = dtk(2)
         phi(i) = -1.0
*        ri2 = x(1,i)**2 + x(2,i)**2 + x(3,i)**2
*        fdum = 1.0/(ri2 + rsoft**2)
*        fddum = fdum/(SQRT(ri2) + rsoft)
*        fdum = 1000.0*fdum
*        fddum = 1000.0*fddum
         do k = 1,3
            f(k,i) = 0.0
            fdot(k,i) = 0.0
         enddo
*        CALL g6_set_j_particle(gpid,gpaddr(i),gpindx(i),t0(i),step(i),
*    &                          body(i),f2dot,fdot(1,i),f(1,i),
*    &                          xdot(1,i),x(1,i)) 
      enddo
*
* Improve initial guess for primordial binary components.
      do ipair = 1,nbin0
         i = 2*ipair - 1
         j = i + 1
         do k = 1,3
            a(k) = x(k,j) - x(k,i)
            a(k+3) = xdot(k,j) - xdot(k,i)
         enddo
*
         a(7) = 1.0/(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))
         a(8) = a(7)*sqrt(a(7))
         a(9) = 3.0*(a(1)*a(4) + a(2)*a(5) + a(3)*a(6))*a(7)
         do k = 1,3
            f(k,i) = f(k,i) + a(k)*a(8)*body(j)
            f(k,j) = f(k,j) - a(k)*a(8)*body(i)
            fdot(k,i) = fdot(k,i) + (a(k+3) - a(k)*a(9))*a(8)*body(j)
            fdot(k,j) = fdot(k,j) - (a(k+3) - a(k)*a(9))*a(8)*body(i)
         enddo
         phi(i) = phi(i) - body(j)*SQRT(a(7))
         phi(j) = phi(j) - body(i)*SQRT(a(7))
      enddo
*
* Include primordial hierarchies (saved sequentially after binaries).
      do l = 1,nhi0
         i = 2*nbin0 + l
         do kcomp = 1,2
            j = 2*l + kcomp - 2
            do k = 1,3
               a(k) = x(k,j) - x(k,i)
               a(k+3) = xdot(k,j) - xdot(k,i)
            enddo
*
            a(7) = 1.0/(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))
            a(8) = body(j)*a(7)*sqrt(a(7))
            a(9) = 3.0*(a(1)*a(4) + a(2)*a(5) + a(3)*a(6))*a(7)
            do k = 1,3
               f(k,i) = f(k,i) + a(k)*a(8)
               fdot(k,i) = fdot(k,i) + (a(k+3) - a(k)*a(9))*a(8)
            enddo
            phi(i) = phi(i) - body(j)*SQRT(a(7))
         enddo
      enddo
*
* Copy new values to GRAPE.
*     do i = 1,2*nbin0+nhi0
*        gpaddr(i) = i - ifirst
*        gpindx(i) = i
*        CALL g6_set_j_particle(gpid,gpaddr(i),gpindx(i),t0(i),step(i),
*    &                          body(i),f2dot,fdot(1,i),f(1,i),
*    &                          xdot(1,i),x(1,i)) 
*     enddo
*
* Set radius of neighbour sphere for GRAPE. 
*     h2 = cmsep2*rmin**2
*     do i = 1,npipe
*        gph2(i) = h2 + 1.0d-10
*     enddo
*
* Calculate F, FDOT & PHI on all single particles twice using GRAPE.
*     iloop = 0
*     ierr = 0
*     CALL g6_set_ti(gpid,time)
*20   i = 1
      i = 1
      do ii = 1,n,npipe
*
         np = n - ii + 1
         if(np.gt.npipe) np = npipe
         do j = 1,np
            jj = ii + j - 1
            gpindx(j) = jj
            do k = 1,3
               xgp(k,j) = x(k,jj)
               vgp(k,j) = xdot(k,jj)
            enddo
         enddo
*
         CALL g6acc(np,gpindx,xgp(1,1),vgp(1,1),
     &                 gpacc(1,1),gpjerk(1,1),gppot(1))
*        CALL g6calc_firsthalf(gpid,nn,np,gpindx(ii),x(1,ii),xdot(1,ii),
*    &                         f(1,ii),fdot(1,ii),phi(ii),eps2,gph2)
*        gperr = g6calc_lasthalf(gpid,nn,np,gpindx(ii),
*    &                           x(1,ii),xdot(1,ii),
*    &                           eps2,gph2,gpacc,gpjerk,gppot)
*
* Check for hardware error. 
*        if(gperr.ne.0)then
*           WRITE(6,*)'GRAPE HARDWARE ERROR: RESET REQUIRED'
*           WRITE(6,*)'g6calc_lasthalf returned ',gperr,' in fpoly0'
*           CALL gpwipe(gpid,time)
*           ierr = ierr + 1
*           if(ierr.ge.50)then
*              WRITE(6,*)' TOO MANY RESETS '
*              WRITE(6,*)' STOP AT TIME = ',time
*              CALL gpfree
*              STOP
*           endif
*
* Send all particles to GRAPE again.
*
*           do j = ifirst,ntot
*              CALL g6_set_j_particle(gpid,
*    &                                gpaddr(j),gpindx(j),t0(j),step(j),
*    &                                body(j),f2dot,fdot(1,j),f(1,j),
*    &                                xdot(1,j),x(1,j)) 
*           enddo
*           goto 20
*        endif
*
* Copy F, FDOT and PHI into COMMON variables.
*
         do ip = 1,np
            do k = 1,3
               f(k,i) = gpacc(k,ip)
               fdot(k,i) = gpjerk(k,ip)
            enddo
            phi(i) = gppot(ip)
            i = i + 1
         enddo
*
      enddo
*     iloop = iloop + 1
*     if(iloop.eq.1) goto 20
*
* Check option for external force.
      if(kz(14).gt.0)then
         do i = 1,n
            CALL xtrnld(i,i,1)
         enddo
      endif
*
* Obtain time-step for all single particles.
      CALL steps(1,n,1)
*
      RETURN
      END
***
