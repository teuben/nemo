C=======================================================================
C
C
C                    INCLUDE file treedefs.inc
C
C=======================================================================
C       includes:
C               nmax            max number of bodies
C               ndim            dimension of problem (2 or 3)
        include 'nmax.inc'
C
C=======================================================================
C
C
C     Parameter declarations, allocation of array storage, common
C     block definitions.
C
C
C=======================================================================

        CHARACTER*20 headline
        INTEGER root,pskip,subp
        LOGICAL usequad
        REAL mass
        
        PARAMETER(nsubcell=2**ndim)
        PARAMETER(nbodsmax=NMAX,ncells=nbodsmax,maxnterm=NMAX)
        PARAMETER(nbits=30,null=0)
        PARAMETER(nbodcell=nbodsmax+ncells,nbods1=nbodsmax+1,
     &            nboddim=nbodsmax*ndim)

        COMMON/paramcom/nbodies,tol,eps,eps2,epsinv,usequad
        COMMON/msgcom/headline
        COMMON/cellcom/rsize,rmin(ndim),incells,imax,imax2,nindex(ndim)
        COMMON/pointers/root,pskip,subp(nbods1:nbodcell,1:nsubcell)
        COMMON/bodycell/mass(1:nbodcell),phi(1:nbodsmax),
     &                  pos(1:nbodcell,1:ndim),sizetol2(1:nbodcell),
     &                  vel(1:nbodsmax,1:ndim),acc(1:nbodsmax,1:ndim)
        COMMON/tempstor/dr(1:maxnterm,1:ndim),drdotdr(1:maxnterm),
     &                  pmass(1:maxnterm),pos0(1:ndim)
        COMMON/forcecom/nterms,nttot,ntmin,ntmax,ntavg
        COMMON/timecom/nsteps,nout,tnow,tpos,tvel,dtime,dtime2
        COMMON/cpucom/cputime0,cputime1,cputime

C=======================================================================
C   Definitions specific to input/output.
C=======================================================================
        INTEGER uterm,upars,ubods,ulog
        CHARACTER*8 parsfile,bodsfile,logfile

        PARAMETER(uterm=6,upars=10,ubods=11,ulog=12)
        PARAMETER(parsfile='treepars',bodsfile='treebods',
     &            logfile='treelog')

C=======================================================================
C   Definitions specific to the quadrupole calculation.
C=======================================================================
        PARAMETER(nquad=2*ndim-1)

        COMMON/quadcom/quad(nbods1:nbodcell,1:nquad)
        COMMON/tempquad/nqterms,qquad(1:maxnterm,1:nquad),
     &                  qdrdotdr(1:maxnterm),qdr(1:maxnterm,1:ndim)
C-----------------------------------------------------------------------


