

      function pot(s,z)
      parameter (pi=3.1415926535)
      common /potconstants/ apot(20,0:20000), frad(20,0:20000), 
     +     dr, nr, lmax, potcor
      common /legendre/ plcon(0:40)
      common /flags/ idiskflag, ibulgeflag, ihaloflag
c     Returns the potential, being the sum of the spherical harmonics
c     with coefficients in apot, and the potential appdiskpot which
c     approximates the high-frequency components of the disk.
      r=sqrt(s*s+z*z)
      ihi=int(r/dr)+1
      if (ihi.lt.1) ihi=1
      if (ihi.gt.nr) ihi=nr
      r1=dr*(ihi-1)
      r2=dr*ihi
      t=(r-r1)/(r2-r1)
      tm1 = 1.0 - t
      if (r.eq.0.) then
         lmaxx=0
         costheta=0
      else
         costheta=z/r
         lmaxx=lmax
      endif
      p=0
      do l=lmaxx,0,-2
         p=p+plgndr1(l,costheta)*plcon(l)*
     +        (t*apot(l/2+1,ihi)+ tm1*apot(l/2+1,ihi-1))
      enddo
      if( idiskflag .eq. 1 ) then
         p = p + appdiskpot(s,z)
      endif
      
      pot=p
      return
      end
