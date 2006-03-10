***
      SUBROUTINE POTN2(SUMP)
*
*
*       Total potential energy on GRAPE (and host for KS).
*       --------------------------------------------------
*
      INCLUDE 'common4.h'
      INTEGER i,k,ii,jj,jmin,ip,nn,np,ierr
      REAL*8 sump,potj,dtmax,a(9)
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
*
* Find the largest block time-step commensurate with current time. 
*
*     dtmax = dtk(1)
*     CALL dtchck(time,dtmax,dtk(40))
*
* Send all single and c.m. particles to GRAPE (no prediction).
*     gpt0 = time
*     do i = ifirst,ntot
*        gpaddr(i) = i - ifirst
*        gpindx(i) = i
*        gpdtj = MIN(dtmax,step(i))
*        CALL g6_set_j_particle(gpid,gpaddr(i),gpindx(i),gpt0,gpdtj,
*    &                          body(i),d2(1,i),fdot(1,i),f(1,i),
*    &                          xdot(1,i),x(1,i))
*     enddo
*
* Sum the potential energy of single particles using GRAPE.
      sump = 0.d0
*     nn = n - ifirst + 1
*
* Calculate F & FDOT on all single particles using GRAPE.
*     CALL g6_set_ti(gpid,time)
*     ierr = 0
*40   i = ifirst
*     do ii = ifirst,n,npipe
      do i = ifirst,n
*
*        np = n - ii + 1
*        if(np.gt.npipe) np = npipe
*
*        CALL g6calc_firsthalf(gpid,nn,np,gpindx(ii),x(1,ii),xdot(1,ii),
*    &                         f(1,ii),fdot(1,ii),phi(ii),eps2,gph2)
*        gperr = g6calc_lasthalf(gpid,nn,np,gpindx(ii),
*    &                           x(1,ii),xdot(1,ii),
*    &                           eps2,gph2,gpacc,gpjerk,gppot)
*
* Check for hardware error.
*        if(gperr.ne.0)then
*           WRITE(6,*)'GRAPE HARDWARE ERROR: RESET REQUIRED'
*           WRITE(6,*)'g6calc_lasthalf returned ',gperr,' in potn2'
*           CALL gpwipe(gpid,time)
*           ierr = ierr + 1
*           if(ierr.ge.10)then
*              WRITE(6,*)' TOO MANY RESETS '
*              WRITE(6,*)' STOP AT TIME = ',time
*              CALL gpfree
*              STOP
*           endif
*           do i = ifirst,ntot
*              gpaddr(i) = i - ifirst
*              gpindx(i) = i
*              CALL g6_set_j_particle(gpid,
*    &                                gpaddr(i),gpindx(i),gpt0,step(i),
*    &                                body(i),d2(1,i),fdot(1,i),f(1,i),
*    &                                xdot(1,i),x(1,i))
*           enddo
*           goto 40
*        endif
*
         phi(i) = 0.d0
         do 50 j = ifirst,n
            if(j.eq.i) goto 50
            do k = 1,3
               a(k) = x(k,j) - x(k,i)
            enddo
            a(7) = 1.d0/(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))
            phi(i) = phi(i) - body(j)*SQRT(a(7))
 50      continue
* Accumulate the sum and copy potential energy to COMMON variable.
         sump = sump - body(i)*phi(i)
      enddo
*
* Set total potential energy from double summation on GRAPE.
      sump = 0.5d0*sump
*
* Calculate the potential of all KS components.
      do i = 1,ifirst-1
         ip = kvec(i)
* Binding energy of regularized pairs included explicitly in EBIN.
         jmin = 2*ip + 1
         potj = 0.d0
         do j = jmin,n
            if(body(j).gt.0.d0)then
               a1 = x(1,i) - x(1,j)
               a2 = x(2,i) - x(2,j)
               a3 = x(3,i) - x(3,j)
               potj = potj + body(j)/SQRT(a1*a1 + a2*a2 + a3*a3)
            endif
         enddo
         phi(i) = -1.d0*potj
         sump = sump + body(i)*potj
      enddo
*
* Add missing terms above IPAIR to the potential of KS components.
      do ip = 1,npairs
         potj = 0.d0
         jmin = 2*ip - 2
         i1 = 2*ip - 1
* Evaluate the potential using first component (sufficient accuracy).
         do j = 1,jmin
            if(body(j).gt.0.d0)then
               a1 = x(1,i1) - x(1,j)
               a2 = x(2,i1) - x(2,j)
               a3 = x(3,i1) - x(3,j)
               potj = potj + body(j)/SQRT(a1*a1 + a2*a2 + a3*a3)
            endif
         enddo
         phi(i1) = phi(i1) - potj
         phi(i1+1) = phi(i1+1) - potj
* Form the c.m. potential by averaging over components.
         phi(n+ip) = 0.5d0*(phi(i1) + phi(i1+1))
      enddo
*
* Add missing KS terms for the single particles (bug 01/2001).
      if(npairs.gt.0)then
         do i = ifirst,n
            if(body(i).eq.0.d0) goto 70
            potj = 0.d0
            do ip = 1,npairs
               j = n + ip
               if(body(j).gt.0.d0)then
                  a1 = x(1,i) - x(1,j)
                  a2 = x(2,i) - x(2,j)
                  a3 = x(3,i) - x(3,j)
                  potj = potj + body(j)/SQRT(a1*a1 + a2*a2 + a3*a3)
               endif
            enddo
            phi(i) = phi(i) - potj
 70         continue
         enddo
      endif
*
* Include optional external potential (note XTRNLV returns energy).
      if(kz(14).gt.0.and.kz(14).lt.3)then
         do i = 1,ntot
            CALL xtrnlv(i,i)
            if(body(i).gt.0.d0)then
               phi(i) = phi(i) + ht/body(i)
            endif
         enddo
      endif
*
* Restore all data on GRAPE.
      isend = 1
      CALL gpsend
*
      RETURN
      END
***
