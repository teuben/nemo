***
      SUBROUTINE POTN0(SUMP)
*
*
*       Total potential energy on host.
*       ---------------------------------
*
      INCLUDE 'common4.h'
      INTEGER i,k,jmin,ip
      REAL*8 sump,potj,a(10)
*
* Sum the potential energy of single particles.
*
      sump = 0.d0
*
      do i = ifirst,n
*
* Accumulate the sum and copy potential energy to COMMON variable.
*
         phi(i) = 0.d0
         do j = ifirst,n
            if(j.eq.i) goto 10
            do k = 1,3
               a(k) = x(k,j) - x(k,i)
               a(k+3) = xdot(k,j) - xdot(k,i)
            enddo
            a(7) = 1.d0/(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))
            a(8) = body(j)*a(7)*SQRT(a(7))
            a(9) = 3.d0*(a(1)*a(4) + a(2)*a(5) + a(3)*a(6))*a(7)
            phi(i) = phi(i) - body(j)*SQRT(a(7))
 10         continue
         enddo
         sump = sump - body(i)*phi(i)
*
      enddo
*
* Set total potential energy from double summation.
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
      if(kz(14).gt.0)then
         do i = 1,ntot
            CALL xtrnlv(i,i)
            if(body(i).gt.0.d0)then
               phi(i) = phi(i) + ht/body(i)
            endif
         enddo
      endif
*
      RETURN
      END
***
