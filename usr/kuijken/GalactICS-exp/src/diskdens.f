
      function diskdens(r,z,psi) 
      parameter (pi=3.1415926535)

      common /potconstants/ apot(20,0:20000), fr(20,0:20000), 
     +     dr, nr, lmax, potcor
      common /gparameters/  a, b, c, v0, q, psi00, 
     +                      psiout, rho1, sigbulge2, 
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /moreconstants/ v02, v03, rdisk2, diskconst, bulgea
      common /flags/ idiskflag, ibulgeflag, ihaloflag

c the truncation factor is derived assuming that the azimuthal velocity dist
c is gaussian with dispersion sigma_R * 2 omega/kappa, as per epicycle theory.
      if( idiskflag .eq. 0 ) then
         diskdens = 0
         return
      endif
      if( abs(z/zdisk) .gt. 30.) goto 99
      erfarg=(r-outdisk)/1.4142136/drtrunc
      if (erfarg.gt.4.) goto 99
      if (erfarg.lt.-4) then
         trunfac=1
      else
         trunfac=0.5*erfc(erfarg)
      endif
      if (z.eq.0.) then
        con=1
      elseif ( ihaloflag .eq. 0 .and. ibulgeflag .eq. 0 ) then
c only doing a disk potential
         con = exp(-abs(z/zdisk))
         con = (2.0*con/(1.0 + con*con))**2
      else
c make disk vertically isothermal; normalize density so that
c at height 3*zdisk the vertical density drop is that of a sech^2.
        psizh=pot(r,3*zdisk)
        psi0=pot(r,0.)
        dpsizh = psizh - psi0
        dpsi   = psi - psi0
        coeff = 0.0
        if( abs(dpsizh) .gt. 0.0) then
            coeff = dpsi/dpsizh
        endif
		if( coeff .gt. 16. .or. coeff .lt. 0.0 ) then
            con = 0.0
        else
            con = 0.009866**coeff
        endif
      endif
      diskdens = diskconst*exp(-r/rdisk)*con*trunfac
      return
 99   diskdens=0.0
      return
      end
