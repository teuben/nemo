

      function bulgedenspsi(psi)
      parameter (pi=3.1415926535, sq8=2.828427125)
      common /gparameters/  a, b, c, v0, q, psi00, 
     +                      psiout, rho1, sigbulge2, 
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /moreconstants/ v02, v03, rdisk2, diskconst, bulgea

      if (psi.gt.psiout) then
         bulgedenspsi=0
      else
         t1=-(psi-psiout)/sigbulge2
         if( t1 .gt. 80) then
            t1 = 80
         endif
c  bulgedenspsi=bulgea*(exp(t1)*erf(t1)-sqrt(4*t1/pi)*(1-0.6666667*t1))
         bulgedenspsi=bulgea*(exp(t1)*erf(sqrt(t1))-
     +   sqrt(4*t1/pi)*(1+0.6666667*t1))
      endif
      return
      end
