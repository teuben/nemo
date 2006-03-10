      SUBROUTINE INTEV(ZMSTAR,TEVDUM,DLOGM)
*
*
*       Chernoff-Weinberg evolution time-scale.
*       ---------------------------------------
*
      include 'common4.h'
cTEVDUM is time in Myr
      REAL*8 xevolg(13),yevolg(13),zmstar,dlogm,tevdum
      data xevolg/-0.08,-0.01,0.07,0.16,0.27,0.40,0.54,0.72,0.91,1.11,
     &     1.33,1.55,1.79/
      data yevolg/10.18,9.93,9.63,9.28,8.90,8.50,8.11,7.68,7.33,7.02,
     &     6.76,6.57,6.50/
*
*       Obtain total evolution time for ZMSTAR in solar masses.
*      TMS = (2.55d+03 + 6.69d+02*ZMSTAR**2.5d0 + ZMSTAR**4.5d0)/
*     &      (3.27d-02*ZMSTAR**1.5d0 + 3.46d-01*ZMSTAR**4.5d0)
*      TG = 0.15d0*TMS
*      THE = TMS*1.37d0*ZMSTAR**(-0.881d0)
*      TEVDUM = TMS + TG + THE
*
c      if (zmstar.le.0.or.tevdum*tstar.ge.1.d5) then 
      if (zmstar.ge.10**xevolg(13)) then
         WRITE (6,*) 'Main sequence mass too great for intev.f;',
     &        ' stopping'
         call gpfree
         stop
      endif
      if (zmstar.le.10**(xevolg(1))) then 
         if (zmstar.le.10**(xevolg(1)-0.04)) then
            dlogm = 0.d0
            tevdum = 1.0d+10
         else
            k = 2
            tevdum=(yevolg(k)-yevolg(k-1))/(xevolg(k)-xevolg(k-1))
            DLOGM = 1.d0/tevdum
            zmslog=log10(zmstar)
            tevdum=tevdum*(zmslog-xevolg(k-1))+yevolg(k-1)
            tevdum=10.d0**tevdum/1.0d+06
         endif
      else
         zmslog=log10(zmstar)
         do 10 k=2,13
            if(zmslog.lt.xevolg(k)) then
               tevdum=(yevolg(k)-yevolg(k-1))/(xevolg(k)-xevolg(k-1))
               DLOGM = 1.d0/tevdum
               tevdum=tevdum*(zmslog-xevolg(k-1))+yevolg(k-1)
               go to 20
            endif
   10    continue
*
   20    tevdum=10.d0**tevdum/1.0d+06
      endif
*
      RETURN
*
      END
