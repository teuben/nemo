C=======================================================================
C
C
C                        INCLUDE FILE treedefs.h                       
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

        CHARACTER*50 headline
        DOUBLE PRECISION mass,tol,tol2inv,eps,rsize,rmin,phi,pos,vel,
     &                   acc,quad,tnow,tpos,dtime,dtime2,tiny,one,two,
     &                   cputime0,cputime1,cputime,zero,cellsize,
     &                   mtot,etot,ektot,eptot,minustwo,
     &                   phsmooth,acsmooth
        INTEGER root,subp,nbodies,incells,nttot,ntmin,ntmax,ntavg,
     &          nsteps,noutbod,noutlog,ndim,nsubcell,nbodsmax,ncells,
     &          nbodcell,nbods1,ninterp
        LOGICAL usequad
     
        PARAMETER(ndim=3,nsubcell=2**ndim)
        PARAMETER(nbodsmax=8192,ncells=8192)
        PARAMETER(nbodcell=nbodsmax+ncells,nbods1=nbodsmax+1)
        PARAMETER(ninterp=30000)

cdeg: tol, tol2inv, eps
c        COMMON/paramcom/nbodies,tol,tol2inv,eps,usequad
        COMMON/paramcom/usequad,nbodies,tol,tol2inv,eps
        COMMON/msgcom/headline
        COMMON/cellcom/rsize,rmin(ndim),incells
        COMMON/pointers/root,subp(nbods1:nbodcell,1:nsubcell)
        COMMON/bodycell/mass(1:nbodcell),phi(1:nbodsmax),
     &                  pos(1:nbodcell,1:ndim),cellsize(1:nbodcell),
     &                  vel(1:nbodsmax,1:ndim),acc(1:nbodsmax,1:ndim)
        COMMON/quadcom/quad(nbods1:nbodcell,1:2*ndim-1)
        COMMON/forcecom/nttot,ntmin,ntmax,ntavg
cdeg: tnow, tpos, dtime, dtime2
        COMMON/timecom/tnow,tpos,dtime,dtime2,nsteps,noutbod,noutlog
        COMMON/cpucom/cputime0,cputime1,cputime
        COMMON/misccom/tiny,zero,one,two,minustwo
        COMMON/enrgycom/mtot,etot,ektot,eptot
        COMMON/gravcom/phsmooth(0:1+ninterp),acsmooth(0:1+ninterp)

C=======================================================================
C   Definitions specific to input/output.
C=======================================================================
        INTEGER uterm,upars,ulog,ubodsin,ubodsout,ubodsasc
        CHARACTER*8 parsfile,logfile,ibodfile,obodfile,oascfile

        PARAMETER(uterm=6,upars=10,ulog=11,ubodsin=12,ubodsout=13,
     &            ubodsasc=14)
        PARAMETER(parsfile='TREEPAR',logfile='TREELOG',
     &            ibodfile='TREEBI',obodfile='TREEBO',
     &            oascfile='TREEBOA')

C=======================================================================
C   Definitions specific to vectorized tree construction, vectorized
C   and tree walk.
C=======================================================================
        DOUBLE PRECISION workvect
        INTEGER subindex,bodlist,templist,parent,asubp,celllist,
     &          subpvect(nsubcell*ncells),isubset,nworkvec

        PARAMETER(nworkvec=9000)

        EQUIVALENCE (subpvect(1),subp(nbods1,1))

        COMMON/concom/celllist(ncells),parent(nbodsmax),asubp(nbodsmax),
     &                templist(nbodsmax),bodlist(nbodsmax),
     &                isubset(nbodsmax),subindex(nbodsmax)
        COMMON/workcom/workvect(nworkvec)

