C=======================================================================
C
C                        INCLUDE FILE magalie.h
C
C modified 18.08.2k cmb ... added axihalo [logical], chalo [real*8]
C radius [real*4] list [integer*4] 
C modified 20.08.2k cmb ... added halomass_in [real*8] 
C modified 29.05.2k+2 cmb ... put axihalo in common block
C=======================================================================
C Parameter declarations, allocation of array storage, common block defs
C=======================================================================

        include 'params.h' 
c
        CHARACTER*2 halotype

        INTEGER nbodies,ndstars,ndgas,ndisk,nbulge,nhalo,ntabhalo,nsat,
     &          nsimpson
        integer*4 list(nbodsmax) 

        LOGICAL outpteps,usegas,selfggas,usebulge,selfgbul,usehalo,
     &       selfghal,usesat,selfgsat,addmods,axibulge,axihalo,bulgerot


        REAL*8 h,diskmass,G,z0,rsolar,qsolar,epsdisk,zmax,rmax,gasmass,        
     &       gastemp,z0gas,zmaxgas,rmaxgas,rmingas,x,y,z,vx,vy,vz,
     &       pmass,xgasmass,epsdisk2,radcyl,radsph,bulgmass,abulge,
     &       halomass,halomass_in,ahalo,gamhalo,rthalo,rhalo,xmhalo,
     &       surfd0, surfd,
     &       ax,ay,az,aradcyl,kappa,sigr,sigphi,sigz,sigt,rotcirc,
     &       rotmean,sigr0,acorr,acorrgas,satmass,asat,xsat,ysat,zsat,
     &       vxsat,vysat,vzsat,axsat,aysat,azsat,radsat,rmaxsat,
     &       rmaxhalo,epshalo,epsbulge,epssat,rmaxbulg,xmod1,ymod1,
     &       zmod1,vxmod1,vymod1,vzmod1,xmod2,ymod2,zmod2,vxmod2,
     &       vymod2,vzmod2,uhalo,potsat,pot,zmaxbulg,brotfrac,cbulge,
     &       thetmod1,phimod1,thetmod2,phimod2, dadrtab, chalo

        REAL*4 radius(nbodsmax), rsquareIS, stretch_fact

        COMMON/nparams/nbodies,ndstars,ndgas,ndisk,nbulge,nhalo,nsat,
     *                 iseed                 
        COMMON/unitspar/h,diskmass,G
        COMMON/lparams/outpteps,usegas,selfggas,usebulge,selfgbul,
     &                 usehalo,selfghal,usesat,selfgsat,addmods,
     &                 axibulge,axihalo,bulgerot

        COMMON/diskscom/z0,rsolar,qsolar,epsdisk,zmax,rmax,epsdisk2
        COMMON/diskgcom/gasmass,gastemp,z0gas,zmaxgas,rmaxgas,rmingas,
     &                  xgasmass,surfd0,surfd(nbodsmax),sigr0,acorr,
     &                  acorrgas

        COMMON/partcom/x(nbodsmax),y(nbodsmax),z(nbodsmax),vx(nbodsmax),
     &                 vy(nbodsmax),vz(nbodsmax),pmass(nbodsmax),
     &                 radcyl(nbodsmax),radsph(nbodsmax),ax(nbodsmax),
     &                 ay(nbodsmax),az(nbodsmax),aradcyl(nbodsmax),
     &                 kappa(nbodsmax),sigr(nbodsmax),sigt(nbodsmax),
     &                 sigz(nbodsmax),sigphi(nbodsmax),pot(nbodsmax),
     &                 rotcirc(nbodsmax),rotmean(nbodsmax)


        COMMON/bulgecom/bulgmass,abulge,epsbulge,rmaxbulg,zmaxbulg,
     &                  brotfrac,cbulge,nsimpson
        
        COMMON/halocom/halomass,ahalo,chalo,gamhalo,rthalo,
     &                 rhalo(maxtabh),xmhalo(maxtabh),rmaxhalo,epshalo,
     &                 uhalo(maxtabh), rsquareIS(maxtabh) 
        COMMON/halocom2/halotype
        COMMON/halocom3/ntabhalo

        COMMON/dvelcom/dadrtab(maxtab)

        COMMON/satcom/satmass,asat,xsat,ysat,zsat,vxsat,vysat,vzsat,
     &                axsat,aysat,azsat,radsat,rmaxsat,epssat,potsat

        COMMON/addmodc/xmod1,ymod1,zmod1,vxmod1,vymod1,vzmod1,
     &                 xmod2,ymod2,zmod2,vxmod2,vymod2,vzmod2,
     &                 thetmod1,phimod1,thetmod2,phimod2














