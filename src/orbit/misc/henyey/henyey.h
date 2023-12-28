c***********************************************************************
c===HENYEY.H
      include 'maxorb.h'
c
      include 'model.h'
c
c=/inval/   Initial conditions
      INTEGER type, norbits
      REAL*8 x0,y0,u0,v0,per,step
      COMMON /inval/x0,y0,u0,v0,per,type,norbits,step
c=/druk/    Parameters to control output
      LOGICAL uit
      INTEGER lun1,lun2
      COMMON /druk/uit,lun1,lun2
c=/odata/
      REAL*8   tau(MAXORB),x(MAXORB),y(MAXORB),u(MAXORB),v(MAXORB)
      COMMON /odata/tau,x,y,u,v


