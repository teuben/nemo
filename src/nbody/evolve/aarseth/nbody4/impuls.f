       SUBROUTINE IMPULS(I,DV)
c     
c     Change in velocity components due to disk shocking
c     --------------------------------------------------
c     
c     2-component disk is assumed (see Chernoff et al. 1986 ApJ, 309,183)
c     ---------------------------
c
c     ***************************************
c     VC is the circular velocity of the cluster in N-body units
c     TH is the angle between the pole of the orbit and the plane of the disk
c     TH=0 for orbits perpendicular to the plane of the disk
c     XOM is the 'oscillation frequency' of the star
c     XK1,XK2,Z1 and Z2 are the parameters of the model of the disk
c     ***************************************
c
c
       INCLUDE 'common4.h'
       DATA PC2GM /7.17706D+10/
       REAL*8 DV(3)
       
c
c      GM = 6.67D-08*1.989D+33
c      PC = 3.0857D+18
c      PC2GM = PC**2/GM
       TH = 0.0

       TEMP=(RBAR**2/ZMBAR)*PC2GM
       XK1=3.47E-9*TEMP*exp(-RG0/3.5)/exp(-8.0/3.5)
       XK2=3.96E-9*TEMP*exp(-RG0/3.5)/exp(-8.0/3.5)
       Z1=175/RBAR
       Z2=550/RBAR
c      VC=220*SQRT(RBAR/ZMBAR)*1E5*SQRT(PC/GM)
       VC=220/VSTAR
 
       PAI=3.1415926
       RR=0.0
       VV=0.0
       VR=0.0
       DO K=1,3
          RR=RR+(X(K,I)-RDENS(K))**2
          VV=VV+XDOT(K,I)**2
          VR=VR+(X(K,I)-RDENS(K))*XDOT(K,I)
       END DO
 
       VR=VR/SQRT(RR)
       VVT=VV-VR**2
       XOM=SQRT(VVT/RR)
c
       XNU1=VC*COS(TH)/(2.*Z1)
       XNU2=VC*COS(TH)/(2.*Z2)
c
       FR1=PAI*XOM/(2.*XNU1)
       FR2=PAI*XOM/(2*XNU2)
       if (fr1.gt.80.) then
          XDS1=0.0
       else if(fr1.lt.1e-32) then
          XDS1=XK1/(XNU1*Z1)
       else   
         XDS1=XK1/(2.*XNU1**2*SINH(FR1)*Z1)*PAI*XOM
       end if
       if (fr2.gt.80.) then
          XDS2=0.0
       else if(fr2.lt.1e-32) then
          XDS2=XK2/(XNU2*Z2)
       else   
         XDS2=XK2/(2.*XNU2**2*SINH(FR2)*Z2)*PAI*XOM
       end if
c
       ZT=(X(3,I)-RDENS(3))*SIN(TH)-(X(2,I)-RDENS(2))*COS(TH)
       DVZT=-ZT*(XDS1+XDS2)
       DV(1)=0.0
       DV(2)=-DVZT*COS(TH)
       DV(3)=DVZT*SIN(TH)
       RETURN
 
       END
