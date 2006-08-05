C=======================================================================
C
C
C                        INCLUDE FILE scf.h
C
C
C=======================================================================
C
C
C     Parameter declarations, allocation of array storage, common
C     block definitions.
C
C
C=======================================================================

        INTEGER nbodsmax,nmax,lmax
 
        PARAMETER(nbodsmax=100000,nmax=6,lmax=4)

        CHARACTER*50 headline
        INTEGER nsteps,noutbod,nbodies,noutlog
        LOGICAL selfgrav,inptcoef,outpcoef,zeroodd,zeroeven,fixacc
        REAL*8 tnow,x,y,z,vx,vy,vz,mass,pot,dtime,G,ax,ay,az,one,pi,
     &         twoopi,onesixth,tpos,tvel,cputime0,cputime1,cputime,
     &         potext,two,zero,tiny

        COMMON/bodscom/x(nbodsmax),y(nbodsmax),z(nbodsmax),vx(nbodsmax),
     &                 vy(nbodsmax),vz(nbodsmax),mass(nbodsmax),
     &                 pot(nbodsmax),ax(nbodsmax),ay(nbodsmax),
     &                 az(nbodsmax),potext(nbodsmax)
        COMMON/parcomi/nbodies,nsteps,noutbod,noutlog
        COMMON/parcomr/dtime,G,one,pi,twoopi,onesixth,two,tiny,zero
        COMMON/parcomc/headline
        COMMON/parcoml/selfgrav,inptcoef,outpcoef,zeroodd,zeroeven,
     &                 fixacc
        COMMON/timecom/tpos,tnow,tvel
        COMMON/cpucom/cputime0,cputime1,cputime

C=======================================================================
C   Definitions specific to input/output.
C=======================================================================
        INTEGER uterm,upars,ulog,ubodsin,ubodsout,utermfil,uoutcoef,
     &          uincoef,ubodsel
        CHARACTER*8 parsfile,logfile,ibodfile,obodfile,termfile,
     &              outcfile,incfile,elfile

        PARAMETER(uterm=6,upars=10,ulog=11,ubodsin=12,ubodsout=13,
     &            utermfil=15,uoutcoef=16,uincoef=17,ubodsel=18)
        PARAMETER(parsfile='SCFPAR',logfile='SCFLOG',
     &            ibodfile='SCFBI',obodfile='SNAPxxxx',
     &            termfile='SCFOUT',outcfile='SCFOCOEF',
     &            incfile='SCFICOEF',elfile='SCFELxxx')
