C***********************************************************************
C
C
                             PROGRAM treecode
C
C
C                      Version 2: November 18, 1988
C             NEMO version:  apr-1990 (ascii/binary output switched)
C
C             Lars Hernquist, Institute for Advanced Study
C
C
C***********************************************************************
C
C
C     A code to evolve self-gravitating systems using the hierarchical
C     tree method developed by Barnes and Hut (Nature 324, 446 [1986])
C     and implemented in FORTRAN by Hernquist (Ap. J. Suppl. 64, 715
C     [1987]; Comp. Phys. Comm. 48, 107 [1988]).  This version has 
C     been optimized for supercomputers and is fully vectorized.  The 
C     code is written in standard FORTRAN, although CRAY-specific 
C     vector intrinsic functions have been used (see below).
C
C     In this version, vectorization of the tree walks is achieved by 
C     simultaneously processing all cells at the same level in the 
C     tree, as discussed by Hernquist (J. Comp. Phys, submitted 
C     [1988]).  The gravitational force calculation proceeds for a 
C     single particle at a time, in serial order.
C
C     The gravitational field is smoothed according to a cubic spline
C     kernel.  The kernel is computed by linear interpolation from a 
C     look-up table.
C
C     A self-starting leap-frog integrator is used, as described by
C     Hernquist and Katz (Ap. J. Suppl., in press [1988]).
C     
C     The computational system of units is determined by the input
C     data, with the assumption that G=1.  Particles are not
C     required to have identical masses.
C
C     This file contains the complete source code for CRAYs.  In
C     order to run on other machines the file treeutil.f must also
C     be linked along with either treeutilsun.f (UNIX machines) or
C     treeutilvax.f (VMS machines).  Common-block variables are
C     defined in the include file treedefs.h.
C
C     Two input data files are required to run this code: a parameter
C     file, which is read in through the subroutine inparams, and a
C     body data file, read by subroutine inbods.  Both are ASCII and
C     their structure is defined in the subroutines which read them.
C     Three output files are created: an ASCII log file, an ASCII
C     body data file containing the final state of the system, and
C     a binary body data file.  The log and binary body files are
C     updated every noutlog and noutbod steps, respectively.
CNEMO-BEGIN
C     For reasons of easy update this file isexactly the same except
C     the names of subroutines 'outbods' and 'outascb' have been
C     swapped, as to make frequent ascii dumps (atos/stoa) and only
C     one binary dump - which is useless for NEMO anyhow.
C                   -- Peter Teuben
CNEMO-END
C
C     WARNINGS -- To avoid excessive overhead, noutlog should be
C                 larger than 1, typically ~ 10, depending on the
C                 number of steps.
C
C                 When compiling on VAX's, avoid using optimization.
C
C     Please report all problems or suggestions for improvements to
C     this code to lars@helios.ucsc.edu
C
C
C=======================================================================
C
C
C     This is the top-level evolution program treecode.  Its tasks are:
C
C          1) to initialize file structures and global variables;
C          2) to input parameters and the initial system state;
C          3) to advance the state of the system for a given number
C             of timesteps;
C          4) to perform a diagnostic analysis of the system at
C             each time step (energy, angular momentum, etc.);
C          5) to periodically record the state of the system;
C          6) and to terminate the simulation and close data files.
C
C
C=======================================================================
C
C
C     Basic global variables/parameters:
C
C          acc         : acceleration components of a body.
C          acsmooth    : table of smoothed gravitational acceleration.
C          cellsize    : linear sizes of cells.
C          cputime     : cpu time (secs) used during the simulation.
C          cputime0    : cumulative cpu time at start of run.
C          cputime1    : cumulative cpu time at end of run.
C          dtime       : the timestep.
C          dtime2      : timestep/2.
C          eps         : gravitational smoothing parameter.
C          ektot       : total system kinetic energy.
C          eptot       : total system gravitational potential energy.
C          etot        : total energy of the system.
C          headline    : identification string for the run.
C          incells     : number of cells currently in use.
C          mass        : masses of bodies and cells.
C          minustwo    : the constant -2.
C          mtot        : total mass of the system.
C          nbodies     : total number of bodies.
C          nbodsmax    : maximum number of bodies.
C          ncells      : maximum number of cells.
C          ndim        : number of spatial dimensions.
C          ninterp     : number of values in look-up tables.
C          noutbod     : frequency of system record outputs.
C          noutlog     : frequency of outputs to log file.
C          nsteps      : total number of timesteps in the simulation.
C          nsubcell    : number of subcells per cell.
C          ntavg       : average length of interaction lists.
C          ntmax       : largest interaction list in current time step.
C          ntmin       : shortest interaction list in current time step.
C          nttot       : sum of interaction lists in current time step.
C          one         : the constant 1.
C          phi         : gravitational potential.
C          phsmooth    : table of look-up values for grav. potential.
C          pos         : coordinates of bodies, c.m. coords. of cells.
C          quad        : quadrupole moments of cells.
C          rmin        : coords. of lower-left corner of system box.
C          root        : pointer to the top of the tree.
C          rsize       : length of the system box.
C          subp        : pointers to descendents of a cell.
C          tiny        : a small number used to prevent divergences.
C          tnow        : current system time.
C          tol         : accuracy parameter.
C          tol2inv     : 1. / (tol * tol).
C          tpos        : current position time.
C          two         : the constant 2.
C          usequad     : option to use (.TRUE.) quadrupole terms.
C          vel         : velocity components of a body.
C          zero        : the constant 0.
C
C-----------------------------------------------------------------------
C
C   Definitions specific to input/output.
C
C          uterm, upars, ulog, ubodsin,   : logical i/o unit numbers.
C             ubodsout,ubodsasc
C          parsfile, logfile, ibodfile,   : character names of files.
C             obodfile,oascfile
C
C
C-----------------------------------------------------------------------
C
C   Definitions specific to vectorized tree construction, vectorized 
C   tree walk, and vectorized tree search for nearest neighbors.  Note
C   that bottom is equivalenced to pos, which is defined above.
C
C          asubp       : subpointers for active bodies or cells.
C          bodlist     : list of active bodies (i.e. not yet leaves).
C          celllist    : list of cells.
C          isubset     : indices of subsets of active bodies or cells.
C          nworkvec    : length of temporary work array workvect.  It
C                        should be set to 9*max length of grav
C                        interaction list.
C          parent      : parents of active bodies or cells.
C          subindex    : subindices for active bodies or cells.
C          subpvect    : vector equivalenced to subp.
C          templist    : temporary vector to swap arrays.
C          workvect    : temporary work array.
C
C   
C=======================================================================
C
C
C     Data structure used to compute gravitational field:
C
C          The body/cell tree structure is assumed to be of the
C          form discussed by Barnes and Hut.  Schematically, for
C          three dimensions (i.e. eight subcells per cell):
C
C         +-------------------------------------------------+
C  root-->| CELL:  mass, pos, quad, /, o, /, /, /, /, o, /  |
C         +----------------------------|--------------|-----+
C                                      |              |
C     +--------------------------------+              |
C     |                                               |
C     |   +----------------------------------+        |
C     +-->| BODY:  mass, pos, vel, acc, phi  |        |
C         +----------------------------------+        |
C                                                     |
C     +-----------------------------------------------+
C     |
C     |   +--------------------------------------------------+
C     +-->| CELL:  mass, pos, quad,  o, /, /, o, /, /, o, /  |
C         +--------------------------|--------|--------|-----+
C                                    |        |        |
C                                   etc.     etc.     etc.
C
C
C          The body/cell information is stored in arrays which
C          incorporate both bodies and cells.  For physical
C          quantities relevant to both bodies and cells, such as
C          mass and position, the array indices range from
C          1 --> nbodsmax + ncells.  For those quantities defined
C          only for bodies, such as velocity, the array indices
C          range from 1 --> nbodsmax.  For information relevant
C          to cells alone, such as pointers to descendants, the
C          array indices range from nbodsmax + 1 --> nbodsmax +
C          ncells.  With this convention, pointers can refer to
C          either bodies or cells without conflict.
C
C          The identification of a given unit as a body or a cell
C          is determined by the pointer to the body/cell.  For a
C          body, p is less than or equal to nbodsmax, while for a
C          cell, p > nbodsmax.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER n

C=======================================================================

C   Initialize state of the system.
C   -------------------------------
        CALL initsys
C            -------

C   Advance system state for a given number of steps.
C   -------------------------------------------------

        DO 100 n=1,nsteps

           CALL stepsys(n)
C               -------

 100    CONTINUE

C   Terminate the simulation.
C   -------------------------
        CALL endrun
C            ------

        STOP
        END
C***********************************************************************
C
C
                        SUBROUTINE accgrav(option)
C
C
C***********************************************************************
C
C
C     Subroutine to compute the gravitational acceleration for all of
C     the bodies.  Vectorization is achieved by processing all of the
C     cells at a given level in the tree simultaneously.  The local
C     variable option indicates whether the code is to compute the
C     potential and/or acceleration.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*4 option
        INTEGER p,nterms

C=======================================================================

C   Initialize the interaction list diagnostics.
C   --------------------------------------------
        nttot=0
        ntmin=nbodsmax
        ntmax=0

C   Main loop over all bodies.
C   --------------------------

        DO 100 p=1,nbodies

C   Establish interaction lists.
C   ----------------------------
           CALL treewalk(p,nterms)
C               --------

C   Compute potential and/or acceleration.
C   --------------------------------------
           CALL gravsum(p,nterms,option)
C               -------

C   Update diagnostics, subtracting self-interaction term.
C   ------------------------------------------------------
           nterms=nterms-1
           nttot=nttot+nterms
           ntmin=MIN(ntmin,nterms)
           ntmax=MAX(ntmax,nterms)

 100    CONTINUE

C   Compute average number of force terms per body.
C   -----------------------------------------------
        ntavg=nttot/nbodies

        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE corrpos(rc)
C
C
C***********************************************************************
C
C
C     Subroutine to apply a correction factor to the positions to
C     maintain second order accuracy when outputting particle data
C     to body data file or when computing energy diagnostics.  The
C     argument rc indicates whether the correction factor is to be
C     applied (correct) or removed (reset).
C  
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*7 rc
        DOUBLE PRECISION dt2,rcsign
        INTEGER p,k

C=======================================================================

        IF(rc.EQ.'correct') THEN
           rcsign=-1.
        ELSE
           rcsign=1.
        ENDIF

C   Loop over all spatial coordinates for all bodies.
C   -------------------------------------------------
  
        dt2=dtime**2

        DO 200 k=1,ndim
           DO 100 p=1,nbodies
              pos(p,k)=pos(p,k)+rcsign*acc(p,k)*dt2/8.
 100       CONTINUE
 200    CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                           SUBROUTINE endrun
C
C
C***********************************************************************
C
C
C     Subroutine to end the simulation.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C=======================================================================

C   Correct positions, output phase space data to ASCII file.
C   ---------------------------------------------------------
        CALL corrpos('correct')
C            -------
        CALL zeropot
C            -------
        CALL gravity('pot ')
C            -------
        CALL outascib
C            --------

C   Finish timing, close data files.
C   --------------------------------
        CALL SECOND(cputime1)

        CALL outcpu
C            ------
        CALL stopout
C            -------

        RETURN
        END
C***********************************************************************
C
C
                            SUBROUTINE energy
C
C
C***********************************************************************
C
C
C     Subroutine to compute diagnostics for the system: total energy,
C     total kinetic energy, total potential energy, angular momentum,
C     center of mass coordinates, and center of mass velocity.  The 
C     local variables cmpos, cmvel, and amvec are center of mass 
C     position and velocity, and total angular momentum of the system, 
C     respectively.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        DOUBLE PRECISION cmpos(ndim),cmvel(ndim),amvec(3)
        INTEGER p,k

C=======================================================================

C   Zero the accumulators for system diagnostics.
C   ---------------------------------------------
        mtot=0.
        ektot=0.
        eptot=0.

        DO 100 k=1,ndim
           cmpos(k)=0.
           cmvel(k)=0.
 100    CONTINUE

        DO 120 k=1,3
           amvec(k)=0.
 120    CONTINUE

C-----------------------------------------------------------------------
C   Loop over bodies to compute system mass and potential energy.
C-----------------------------------------------------------------------

        DO 150 p=1,nbodies
           mtot=mtot+mass(p)
           eptot=eptot+.5*mass(p)*phi(p)
 150    CONTINUE

C-----------------------------------------------------------------------
C   Compute system kinetic energy, components of center of mass
C   position and velocity.
C-----------------------------------------------------------------------

        DO 250 k=1,ndim
           DO 200 p=1,nbodies
              ektot=ektot+.5*mass(p)*vel(p,k)*vel(p,k)
              cmpos(k)=cmpos(k)+mass(p)*pos(p,k)
              cmvel(k)=cmvel(k)+mass(p)*vel(p,k)
 200       CONTINUE
           cmvel(k)=cmvel(k)/mtot
           cmpos(k)=cmpos(k)/mtot
 250    CONTINUE

C   Compute total system energy.
C   ----------------------------
        etot=ektot+eptot

C-----------------------------------------------------------------------
C   Compute angular momentum of the system.
C-----------------------------------------------------------------------

        IF(ndim.EQ.2) THEN

           DO 300 p=1,nbodies
              amvec(3)=amvec(3)+mass(p)*(pos(p,1)*vel(p,2)-
     &                 pos(p,2)*vel(p,1))
 300       CONTINUE

        ELSE IF(ndim.EQ.3) THEN

           DO 400 p=1,nbodies
              amvec(1)=amvec(1)+mass(p)*(pos(p,2)*vel(p,3)-
     &                 pos(p,3)*vel(p,2))
              amvec(2)=amvec(2)+mass(p)*(pos(p,3)*vel(p,1)-
     &                 pos(p,1)*vel(p,3))
              amvec(3)=amvec(3)+mass(p)*(pos(p,1)*vel(p,2)-
     &                 pos(p,2)*vel(p,1))
 400       CONTINUE

        ENDIF

C   Write diagnostics to the log file.
C   ----------------------------------
        CALL outenrgy(amvec,cmpos,cmvel)
C            --------

        RETURN
        END
C***********************************************************************
C
C
                       SUBROUTINE gravity(option)
C
C
C***********************************************************************
C
C
C     Subroutine to compute gravitational potential and acceleration.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*4 option

C=======================================================================

        CALL maketree
C            --------
        CALL accgrav(option)
C            -------

        RETURN
        END
C***********************************************************************
C
C
                  SUBROUTINE gravsum(p,nterms,option)
C
C
C***********************************************************************
C
C
C     Subroutine to compute the monopole and quadrupole contributions
C     to the potential and acceleration components for body p.  The
C     interaction list is contained in the vector iterms, which is
C     equivalenced to the common array bodlist.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER maxnterm

        PARAMETER(maxnterm=nworkvec/9)

        CHARACTER*4 option
        DOUBLE PRECISION r3inveff(maxnterm),rinveff(maxnterm),drdeldrg,
     &                   pmass,drdotdr(maxnterm),phsm,drsm,accsm,
     &                   dx(maxnterm),dy(maxnterm),dz(maxnterm),
     &                   qr5inv(maxnterm),phiquad(maxnterm),sdrdotdr,
     &                   r2inveff(maxnterm),acci,epsp,epsi
        INTEGER p,i,qindex(nbodsmax),qterms(nbodsmax),smindex(nbodsmax),
     &          nterms,iterms(nbodsmax),nqterms

        EQUIVALENCE (iterms(1),bodlist(1)),(qindex(1),templist(1)),
     &              (qterms(1),isubset(1)),(smindex(1),parent(1)),
     &              (dx(1),workvect(1)),(dy(1),workvect(maxnterm+1)),
     &              (dz(1),workvect(2*maxnterm+1)),(r3inveff(1),
     &              workvect(3*maxnterm+1)),(rinveff(1),
     &              workvect(4*maxnterm+1)),(drdotdr(1),
     &              workvect(5*maxnterm+1)),(qr5inv(1),
     &              workvect(6*maxnterm+1)),(phiquad(1),
     &              workvect(7*maxnterm+1)),(r2inveff(1),
     &              workvect(8*maxnterm+1))

C=======================================================================

        IF(nterms.GT.maxnterm)
     &     CALL terror(' array overflow in gravsum ')
C               ------

C-----------------------------------------------------------------------
C   Compute monopole contribution; temporarily set mass of body p to
C   zero to avoid possible self-interaction contribution.
C-----------------------------------------------------------------------

        pmass=mass(p)
        mass(p)=0.
        epsp=eps

C   Loop over interaction list.
C   ---------------------------

CDIR$ IVDEP
        DO 30 i=1,nterms
           dx(i)=pos(p,1)-pos(iterms(i),1)
           dy(i)=pos(p,2)-pos(iterms(i),2)
           dz(i)=pos(p,3)-pos(iterms(i),3)
           drdotdr(i)=dx(i)**2+dy(i)**2+dz(i)**2
           sdrdotdr=SQRT(drdotdr(i))
           rinveff(i)=1./(sdrdotdr+tiny)
           r3inveff(i)=rinveff(i)/(drdotdr(i)+tiny)
           epsi=eps
           drdeldrg=sdrdotdr*ninterp/(epsp+epsi)
           smindex(i)=drdeldrg
           smindex(i)=MIN(ninterp,smindex(i))
           drsm=MIN(one,drdeldrg-smindex(i))
           phsm=(1.-drsm)*phsmooth(smindex(i))+
     &          drsm*phsmooth(1+smindex(i))
           accsm=(1.-drsm)*acsmooth(smindex(i))+
     &           drsm*acsmooth(1+smindex(i))
           rinveff(i)=phsm*rinveff(i)
           r3inveff(i)=accsm*r3inveff(i)
 30     CONTINUE

        IF(option.NE.'acc ') THEN

CDIR$ IVDEP
           DO 40 i=1,nterms
              phi(p)=phi(p)-mass(iterms(i))*rinveff(i)
 40        CONTINUE

        ENDIF

        IF(option.NE.'pot ') THEN

CDIR$ IVDEP
           DO 50 i=1,nterms
              acci=mass(iterms(i))*r3inveff(i)
              acc(p,1)=acc(p,1)-dx(i)*acci
              acc(p,2)=acc(p,2)-dy(i)*acci
              acc(p,3)=acc(p,3)-dz(i)*acci
 50        CONTINUE

        ENDIF

C   Reset mass of body p.
C   ---------------------
        mass(p)=pmass
 
C   If required, compute quadrupole contribution.
C   ---------------------------------------------
        IF(usequad) THEN

C   Filter out bodies.
C   ------------------

           CALL WHENIGT(nterms,iterms,1,nbodsmax,qindex,nqterms)
 
C   Compute quadrupole interaction from cells.
C   ------------------------------------------
CDIR$ IVDEP
           DO 60 i=1,nqterms
              qterms(i)=iterms(qindex(i))
              r2inveff(i)=rinveff(qindex(i))*rinveff(qindex(i))
              qr5inv(i)=r3inveff(qindex(i))*r2inveff(i)
              phiquad(i)=(-.5*((dx(qindex(i))**2-dz(qindex(i))**2)*
     &              quad(qterms(i),1)+(dy(qindex(i))**2-
     &              dz(qindex(i))**2)*quad(qterms(i),4))-
     &              (dx(qindex(i))*dy(qindex(i))*quad(qterms(i),2)+
     &              dx(qindex(i))*dz(qindex(i))*quad(qterms(i),3)+
     &              dy(qindex(i))*dz(qindex(i))*quad(qterms(i),5)))*
     &              qr5inv(i)
 60        CONTINUE

           IF(option.NE.'acc ') THEN

CDIR$ IVDEP
              DO 70 i=1,nqterms
                 phi(p)=phi(p)+phiquad(i)
 70           CONTINUE

           ENDIF

           IF(option.NE.'pot ') THEN

CDIR$ IVDEP
              DO 80 i=1,nqterms
                 phiquad(i)=5.*phiquad(i)*r2inveff(i)
                 acc(p,1)=acc(p,1)+dx(qindex(i))*phiquad(i)+
     &                  (dx(qindex(i))*quad(qterms(i),1)+
     &                  dy(qindex(i))*quad(qterms(i),2)+
     &                  dz(qindex(i))*quad(qterms(i),3))*qr5inv(i)
                 acc(p,2)=acc(p,2)+dy(qindex(i))*phiquad(i)+
     &                  (dy(qindex(i))*quad(qterms(i),4)+
     &                  dx(qindex(i))*quad(qterms(i),2)+
     &                  dz(qindex(i))*quad(qterms(i),5))*qr5inv(i)
                 acc(p,3)=acc(p,3)+dz(qindex(i))*phiquad(i)+
     &                  (dz(qindex(i))*(-quad(qterms(i),1)-
     &                  quad(qterms(i),4))+dx(qindex(i))*
     &                  quad(qterms(i),3)+dy(qindex(i))*
     &                  quad(qterms(i),5))*qr5inv(i)
 80           CONTINUE

           ENDIF

        ENDIF
 
        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE hackcell
C
C
C***********************************************************************
C
C
C     Subroutine to compute masses, center of mass coordinates,
C     and optional quadrupole moments of cells, processing cells
C     in order of increasing size.  The permutation vector is
C     returned in the common variable celllist.  Vectorization is
C     achieved by simultaneously processing all cells at the
C     same level in the hierarchy (see Hernquist, J. Comput. Phys.,
C     submitted [1988]).
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,fcell,lcell,i,j,k,l,m,n,nsubc,nnodes

C=======================================================================
        
C   Generate permutation of cells, according to cellsize.
C   -----------------------------------------------------

        DO 5 i=1,incells
           celllist(i)=nbodsmax+incells-(i-1)
 5      CONTINUE

C   Initialize properties of cells.
C   -------------------------------

        DO 10 p=nbodsmax+1,nbodsmax+incells
           mass(p)=0.
           pos(p,1)=0.
           pos(p,2)=0.
           pos(p,3)=0.
 10     CONTINUE

        IF(usequad) THEN
           DO 30 k=1,2*ndim-1
              DO 20 p=nbodsmax+1,nbodsmax+incells
                 quad(p,k)=0.
 20           CONTINUE
 30        CONTINUE
        ENDIF

C   Process cells in order of increasing size.
C   ------------------------------------------

        fcell=1

 40     CONTINUE

        IF(fcell.LE.incells) THEN

C   Determine which cells to process.
C   ---------------------------------

           DO 50 i=fcell,incells
              IF(ABS(cellsize(celllist(i))-cellsize(celllist(fcell)))
     &           .LT.0.01*cellsize(celllist(fcell))) THEN

                 lcell=i
              ELSE
                 GO TO 60
              ENDIF
 50        CONTINUE                    

 60        CONTINUE

C   Compute properties of the selected cells, looping over subcells.
C   ----------------------------------------------------------------

           DO 110 j=1,nsubcell

              DO 70 i=fcell,lcell
                 asubp(i-fcell+1)=subp(celllist(i),j)
 70           CONTINUE

              CALL WHENIGT(lcell-fcell+1,asubp,1,0,isubset,nnodes)

CDIR$ IVDEP
              DO 80 i=1,nnodes
                 parent(i)=celllist(isubset(i)+fcell-1)
                 asubp(i)=subp(parent(i),j)
                 mass(parent(i))=mass(parent(i))+mass(asubp(i))
                 pos(parent(i),1)=pos(parent(i),1)+mass(asubp(i))*
     &                            pos(asubp(i),1)
                 pos(parent(i),2)=pos(parent(i),2)+mass(asubp(i))*
     &                            pos(asubp(i),2)
                 pos(parent(i),3)=pos(parent(i),3)+mass(asubp(i))*
     &                            pos(asubp(i),3)
 80           CONTINUE

 110       CONTINUE

CDIR$ IVDEP
           DO 120 i=fcell,lcell
              pos(celllist(i),1)=pos(celllist(i),1)/mass(celllist(i))
              pos(celllist(i),2)=pos(celllist(i),2)/mass(celllist(i))
              pos(celllist(i),3)=pos(celllist(i),3)/mass(celllist(i))
 120       CONTINUE

C   Compute optional quadrupole moments.
C   ------------------------------------

           IF(usequad) THEN

              DO 210 j=1,nsubcell

                 DO 130 i=fcell,lcell
                    asubp(i-fcell+1)=subp(celllist(i),j)
 130             CONTINUE

                 CALL WHENIGT(lcell-fcell+1,asubp,1,0,isubset,nnodes)

CDIR$ IVDEP
                 DO 140 i=1,nnodes
                    parent(i)=celllist(isubset(i)+fcell-1)
                    asubp(i)=subp(parent(i),j)
 140             CONTINUE

                 CALL WHENIGT(nnodes,asubp,1,nbodsmax,isubset,nsubc)

                 DO 200 m=1,MIN(2,ndim)
                    DO 190 n=m,ndim

                       l=(m-1)*(ndim-1)+n

CDIR$ IVDEP
                       DO 150 i=1,nnodes
                          quad(parent(i),l)=quad(parent(i),l)+
     &                       mass(asubp(i))*(3.*(pos(asubp(i),m)-
     &                       pos(parent(i),m))*(pos(asubp(i),n)-
     &                       pos(parent(i),n)))
 150                   CONTINUE

                       IF(m.EQ.n) THEN
                          DO 170 k=1,ndim
CDIR$ IVDEP
                             DO 160 i=1,nnodes
                                quad(parent(i),l)=quad(parent(i),l)-
     &                             mass(asubp(i))*(pos(asubp(i),k)-
     &                             pos(parent(i),k))**2
 160                         CONTINUE
 170                      CONTINUE
                       ENDIF

CDIR$ IVDEP
                       DO 180 i=1,nsubc
                          templist(i)=parent(isubset(i))
                          quad(templist(i),l)=quad(templist(i),l)+
     &                       quad(asubp(isubset(i)),l)
 180                   CONTINUE

 190                CONTINUE
 200             CONTINUE

 210          CONTINUE

           ENDIF

           fcell=lcell+1
             
           GO TO 40

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                           SUBROUTINE inbods
C
C
C***********************************************************************
C
C
C     Subroutine to read in the data associated with the bodies from
C     the ASCII input file.  The records are assumed to be in the 
C     following order:
C
C                    nbodies
C                    ndim
C                    tnow
C                    mass(1)
C                      .
C                      .
C                      .
C                    mass(nbodies)
C                       x(1)       y(1)       z(1)
C                        .          .          .
C                        .          .          .
C                        .          .          .
C                    x(nbodies) y(nbodies) z(nbodies)
C                       vx(1)       vy(1)       vz(1)
C                         .           .           .
C                         .           .           .
C                         .           .           .
C                    vx(nbodies) vy(nbodies) vz(nbodies)
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,ndimi,k

C=======================================================================
 
C   Read in body data.
C   ------------------

        READ(ubodsin,*) nbodies,ndimi,tnow
 
        IF(nbodies.GT.nbodsmax.OR.ndimi.NE.ndim)
     &     CALL terror(' error in inbods--inconsistent inputs ')
C               ------

        READ(ubodsin,*) (mass(p),p=1,nbodies),((pos(p,k),k=1,
     &                  ndim),p=1,nbodies),((vel(p,k),k=1,ndim),
     &                  p=1,nbodies)

        CLOSE(ubodsin)
 
        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initpars
C
C
C***********************************************************************
C
C
C     Subroutine to initialize system parameters that depend on
C     either the input data or defined PARAMETERS.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        DOUBLE PRECISION xw,xw2,deldrg,xw3,xw4,cvmgt
        INTEGER p,i

C=======================================================================

C   Initialize misc. useful numbers.
C   --------------------------------
        minustwo = -2.
        tiny=1.e-20
        zero=0.
        one=1.
        two=2.

C   Initialize position and velocity times, 1/2 timestep.
C   -----------------------------------------------------
        tpos=tnow
        dtime2=.5*dtime
        tol2inv=1./(tol*tol)

C   Initialize size parameter for bodies.
C   -------------------------------------
        DO 5 p=1,nbodies
           cellsize(p)=0.
 5      CONTINUE

C-----------------------------------------------------------------------
C   Initialize variables and arrays for gravitational field smoothing 
C   interpolation.  Interpolation performed in distance.
C-----------------------------------------------------------------------
        deldrg=2./ninterp

        DO 30 i=0,1+ninterp
           xw=i*deldrg
           xw2=xw*xw
           xw3=xw2*xw
           xw4=xw2*xw2
           phsmooth(i)=CVMGT(-2.*xw3*(one/3.-3.*xw2/20.+xw3/20.)+
     &                       7.*xw/5.,-one/15.+8.*xw/5.-xw3*(4./3.-xw+
     &                       3.*xw2/10.-xw3/30.),xw.LE.one)
           phsmooth(i)=CVMGT(one,phsmooth(i),xw.GE.two)
           acsmooth(i)=CVMGT(xw3*(4./3.-6.*xw2/5.+0.5*xw3),-one/15.+
     &                       8.*xw3/3.-3.*xw4+6.*xw3*xw2/5.-xw4*xw2/6.,
     &                       xw.LE.one)
           acsmooth(i)=CVMGT(one,acsmooth(i),xw.GE.two)
 30     CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initpos
C
C
C***********************************************************************
C
C
C     Subroutine to apply correction factor to make the leap-frog
C     algorithm self-starting (Hernquist and Katz, Ap. J. Suppl.,
C     in press [1988]).
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        DOUBLE PRECISION dt2
        INTEGER p,k

C=======================================================================

C   Loop over all spatial coordinates for all bodies.
C   -------------------------------------------------
  
        dt2=dtime**2

        DO 200 k=1,ndim
           DO 100 p=1,nbodies
              pos(p,k)=pos(p,k)+acc(p,k)*dt2/8.
 100       CONTINUE
 200    CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initsys
C
C
C***********************************************************************
C
C
C     Subroutine to initialize the state of the system.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C=======================================================================

C   Begin timing.
C   -------------
        CALL SECOND(cputime0)

C-----------------------------------------------------------------------
C   Open data files, read input parameters and initial system state,
C   and initialize system parameters.
C-----------------------------------------------------------------------
        CALL startout
C            --------
        CALL inparams
C            --------
        CALL inbods
C            ------
        CALL initpars
C            --------

C   Compute gravitational potential and acceleration.
C   -------------------------------------------------
        CALL zeroacc
C            -------
        CALL zeropot
C            -------
        CALL gravity('both')
C            -------

C   Output system state.
C   --------------------
        CALL outstate(0)
C            --------

C   Correct positions.
C   ------------------
        CALL initpos
C            -------

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE inparams
C
C
C***********************************************************************
C
C
C     Subroutine to read in parameters.
C
C     Input parameters:
C
C        headline  : identification string for the run.
C        nsteps    : number of timesteps.
C        noutbod   : output system state once every nsteps/noutbod 
C                    steps.
C        noutlog   : output logfile data once every nsteps/noutlog
C                    steps.
C        dtime     : the timestep.
C        tol       : error tolerance; 0.0 => exact (PP) calculation.
C        eps       : potential softening parameter.
C        usequad   : option to include (.TRUE.) quadrupole terms.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C=======================================================================
 
C   Read parameters, close the file.
C   --------------------------------
        READ(upars,'(a)') headline
        READ(upars,*) nsteps
        READ(upars,*) noutbod
        READ(upars,*) noutlog
        READ(upars,*) dtime
        READ(upars,*) tol
        READ(upars,*) eps
        READ(upars,*) usequad

        CLOSE(UNIT=upars)
 
        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE loadtree
C
C
C***********************************************************************
C
C
C     Subroutine to insert the bodies into the tree.  The process is
C     vectorized over active bodies (Makino, J. Comput. Phys., 
C     submitted [1988]).  Active bodies are those which are not yet in 
C     place in the tree, as leaves.  The local variables pm1 and nindex
C     are used to convert back and forth between physical coordinates 
C     and subcell coordinates.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        DOUBLE PRECISION pm1(nsubcell,ndim)
        INTEGER k,p,nindex(ndim),j,i,nbodlist,nclist,nclist2,
     &          nsubset,indcell,nsubbod1,nbodtemp

        SAVE nindex,pm1

        DATA pm1/4*-1.,4*1.,2*-1.,2*1.,2*-1.,2*1.,-1.,1.,-1.,1.,
     &            -1.,1.,-1.,1./,nindex/4,2,1/

C=======================================================================

C   Deallocate old tree, compute coordinates of center of root cell.
C   ----------------------------------------------------------------
        incells=1
        root=nbodsmax+1

        DO 5 j=1,nsubcell
           subp(root,j)=0
 5      CONTINUE

        cellsize(root)=rsize

        DO 10 k=1,ndim
           pos(root,k)=rmin(k)+0.5*rsize
 10     CONTINUE

C-----------------------------------------------------------------------
C   Place all bodies on active body list, having root as parent; place
C   root on active cell list.
C-----------------------------------------------------------------------

        DO 20 i=1,nbodies
           parent(i)=root
           bodlist(i)=i
 20     CONTINUE

        nbodlist=nbodies
        celllist(1)=root
        nclist=1

C   Loop until no bodies are left active.
C   -------------------------------------

 200    CONTINUE

        IF(nclist.GT.0) THEN

C   Compute subindices for all active bodies.
C   -----------------------------------------
           DO 30 i=1,nbodlist
              subindex(i)=1
 30        CONTINUE

           DO 50 k=1,ndim
              DO 40 i=1,nbodlist
                 IF(pos(bodlist(i),k).GE.pos(parent(i),k)) 
     &                  subindex(i)=subindex(i)+nindex(k)
 40           CONTINUE
 50        CONTINUE

C   Compute number of bodies in each subcell.
C   -----------------------------------------
           DO 60 i=1,nbodlist
              subp(parent(i),subindex(i))=subp(parent(i),subindex(i))+1
 60        CONTINUE

C-----------------------------------------------------------------------
C   Open all subcells with more than one body, placing them on active 
C   cell list.
C-----------------------------------------------------------------------
           nclist2=0

           DO 110 j=1,nsubcell

              DO 70 i=1,nclist
                 asubp(i)=subp(celllist(i),j)
 70           CONTINUE

              CALL WHENIGT(nclist,asubp,1,1,isubset,nsubset)

              incells=incells+nsubset
              IF(incells.GT.ncells) CALL terror(' overflow in loadtree')
C                                        ------
              indcell=incells-nsubset+nbodsmax

              DO 90 k=1,nsubcell
                 DO 80 i=1,nsubset
                    subp(indcell+i,k)=0
 80              CONTINUE
 90           CONTINUE

CDIR$ IVDEP
              DO 100 i=1,nsubset
                 p=indcell+i
                 asubp(i)=celllist(isubset(i))
                 subp(asubp(i),j)=p
                 cellsize(p)=cellsize(asubp(i))*0.5
                 templist(nclist2+i)=p
                 pos(p,1)=pos(asubp(i),1)+pm1(j,1)*0.5*cellsize(p)
                 pos(p,2)=pos(asubp(i),2)+pm1(j,2)*0.5*cellsize(p)
                 pos(p,3)=pos(asubp(i),3)+pm1(j,3)*0.5*cellsize(p)
 100          CONTINUE

              nclist2=nclist2+nsubset

 110       CONTINUE

           nclist=nclist2
        
           DO 120 i=1,nclist
              celllist(i)=templist(i)
 120       CONTINUE

C   Find all subcells with one body; add bodies to tree.
C   ----------------------------------------------------
           DO 130 i=1,nbodlist
              templist(i)=ncells*(subindex(i)-1)+(parent(i)-nbodsmax)
              asubp(i)=subpvect(templist(i))
 130       CONTINUE

           CALL WHENEQ(nbodlist,asubp,1,1,isubset,nsubbod1)

CDIR$ IVDEP
           DO 140 i=1,nsubbod1
              subpvect(templist(isubset(i)))=bodlist(isubset(i))
 140       CONTINUE

C   Place bodies in cells with more than one body on active list.
C   ------------------------------------------------------------

           CALL WHENIGT(nbodlist,asubp,1,1,isubset,nbodtemp)

           nbodlist=nbodtemp

           DO 150 i=1,nbodlist
              parent(i)=asubp(isubset(i))
              templist(i)=bodlist(isubset(i))
 150       CONTINUE

           DO 160 i=1,nbodlist
              bodlist(i)=templist(i)
 160       CONTINUE

           GO TO 200
        
        ENDIF 

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE maketree
C
C
C***********************************************************************
C
C
C     Main routine to control initialization of the tree structure 
C     for computing the gravitational interaction.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C=======================================================================

C   Set box properties.
C   -------------------
        CALL setbox
C            ------
 
C   Load bodies into the tree.
C   --------------------------
        CALL loadtree
C            --------

C   Compute properties of cells.
C   ----------------------------
        CALL hackcell
C            --------
 
        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE outbods
C                              (old: outascib)
C
C
C***********************************************************************
C
C
C     Subroutine to output the body data to an ascii data file.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,ndimo

C=======================================================================
 
C   Output system state.
C   --------------------

        ndimo=ndim

        WRITE(ubodsasc,*) nbodies
        WRITE(ubodsasc,*) ndimo
        WRITE(ubodsasc,*) tnow

        DO 10 p=1,nbodies
           WRITE(ubodsasc,*) mass(p)
 10     CONTINUE

        DO 20 p=1,nbodies
           WRITE(ubodsasc,*) pos(p,1),pos(p,2),pos(p,3)
 20     CONTINUE

        DO 30 p=1,nbodies
           WRITE(ubodsasc,*) vel(p,1),vel(p,2),vel(p,3)
 30     CONTINUE

        DO 40 p=1,nbodies
           WRITE(ubodsasc,*) phi(p)
 40     CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE outascib
C                              (old: outbods)
C
C
C***********************************************************************
C
C
C     Subroutine to output the body data to binary output file.
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,ndimo,k

C=======================================================================
 
C   Output system state.
C   --------------------

        ndimo=ndim

        WRITE(ubodsout) nbodies,ndimo,tnow
 
        WRITE(ubodsout) (mass(p),p=1,nbodies),((pos(p,k),p=1,
     &                  nbodies),k=1,ndim),((vel(p,k),p=1,nbodies),
     &                  k=1,ndim),(phi(p),p=1,nbodies)

        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE outcpu
C
C
C***********************************************************************
C
C
C     Subroutine to output cpu timing data to the log file.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C=======================================================================
 
C   Output timing data to the log file.
C   -----------------------------------
        cputime=cputime1-cputime0
        WRITE(ulog,10) cputime
 
 10     FORMAT(//,' Total cpu time used (seconds) : ',1pe12.4)
 
        RETURN
        END
C***********************************************************************
C
C
                  SUBROUTINE outenrgy(am,cmpos,cmvel)
C
C
C***********************************************************************
C
C
C     Subroutine to output diagnostic data to the log file.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        DOUBLE PRECISION am(3),cmpos(ndim),cmvel(ndim),cpunew,
     &                   cpuold,cpustep
        INTEGER k

        SAVE cpuold

        DATA cpuold/0.0/

C=======================================================================

C-----------------------------------------------------------------------
C   Write mass, energy, angular momentum, center of mass quantities.
C-----------------------------------------------------------------------
 
        WRITE(ulog,10) mtot
        WRITE(ulog,20) etot,ektot,eptot
        WRITE(ulog,40) am(1),am(2),am(3)
        WRITE(ulog,50) (cmpos(k),k=1,ndim)
        WRITE(ulog,60) (cmvel(k),k=1,ndim)

        CALL SECOND(cpunew)
        cpustep=cpunew-cpuold
        cpuold=cpunew
        WRITE(ulog,70) cpustep
 
 10     FORMAT(/,7x,'mtot = ',4(1pe12.4))
 20     FORMAT(7x,'e, ek, ep = ',4(1pe12.4))
 40     FORMAT(7x,'amx, amy, amz = ',3(1pe12.4))
 50     FORMAT(7x,'cmpos = ',3(1pe12.4))
 60     FORMAT(7x,'cmvel = ',3(1pe12.4))
 70     FORMAT(/,10x,'cpu time per step = ',1pe12.4,/)
 
        RETURN
        END
C***********************************************************************
C
C
                       SUBROUTINE outerror(message)
C
C
C***********************************************************************
C
C
C     Subroutine to output error messages to the log file.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*(*) message

C=======================================================================

C   Write the message.
C   ------------------
        WRITE(ulog,40)
 40     FORMAT(/,1x,72('*'))
        WRITE(ulog,50) message
 50     FORMAT(/,a)
        WRITE(ulog,40)
 
        RETURN
        END
C***********************************************************************
C
C
                      SUBROUTINE outhead(outunit)
C
C
C***********************************************************************
C
C
C     Subroutine to output a standard header to a logical device
C     unit specified by outunit.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER outunit

C=======================================================================
 
        WRITE(outunit,10)
        WRITE(outunit,20)
        WRITE(outunit,20)
        WRITE(outunit,25) headline
        WRITE(outunit,20)
        WRITE(outunit,10)
        WRITE(outunit,20)
        WRITE(outunit,28)
        WRITE(outunit,29)
        WRITE(outunit,30) nbodies,nsteps,noutbod,noutlog
        WRITE(outunit,20)
        WRITE(outunit,40) dtime,eps,usequad,tol
        WRITE(outunit,20)
        WRITE(outunit,10)
 
 10     FORMAT(1x,72('*'))
 20     FORMAT(1x,'*',70(' '),'*')
 25     FORMAT(1x,'*',10x,1a50,10x,'*')
 28     FORMAT(1x,'*',4x,'Input parameters:',49x,'*')
 29     FORMAT(1x,'*',4x,'----------------',50x,'*')
 30     FORMAT(1x,'*',8x,'nbodies=',1i7,2x,'nsteps=',1i4,4x,
     &         'noutbods=',1i4,2x,'noutlog=',1i4,3x,'*')
 40     FORMAT(1x,'*',8x,'dtime=',1pe10.3,1x,'eps=',1pe10.3,1x,
     &         'usequad=',1l1,6x,'tol=',0pf5.2,6x,'*')

        RETURN
        END
C***********************************************************************
C
C
                        SUBROUTINE outlog(istep)
C
C
C***********************************************************************
C
C
C     Subroutine to monitor status of the program by writing to
C     the log file.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER istep

C=======================================================================
 
C   If first call, write header.
C   ----------------------------
        IF(istep.EQ.0) CALL outhead(ulog)
C                           -------
 
C-----------------------------------------------------------------------
C   Output system time and force evaluation diagnostics.
C-----------------------------------------------------------------------
 
        WRITE(ulog,75) tnow,incells
        WRITE(ulog,80) nttot,ntmin,ntmax,ntavg
 
 75     FORMAT(//,2x,'time: ',1pe12.4,5x,'ncells: ',1i5,/)
 80     FORMAT(7x,'nttot, min, max, avg = ',1i8,5x,1i5,5x,1i5,5x,1i5)

        RETURN
        END
C***********************************************************************
C
C
                        SUBROUTINE outstate(n)
C
C
C***********************************************************************
C
C
C     Subroutine to output information about the system state to
C     the log and body data files.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER n

C=======================================================================

        CALL outterm(' step completed: ',n)
C            -------

        IF(n.EQ.0) THEN

           CALL outlog(0)
C               ------
           CALL energy
C               ------
           CALL outbods
C               -------

        ELSE

           IF(MOD(n,noutbod).EQ.0.OR.MOD(n,noutlog).EQ.0) THEN

              CALL corrpos('correct')
C                  -------
              CALL zeropot
C                  -------
              CALL gravity('pot ')
C                  -------

              IF(MOD(n,noutlog).EQ.0) THEN

                 CALL outlog(n)
C                     ------
                 CALL energy
C                     ------
              ENDIF

              IF(MOD(n,noutbod).EQ.0) CALL outbods
C                                          -------
              CALL corrpos('reset  ')
C                  -------
           ENDIF

        ENDIF
 
        RETURN
        END
C***********************************************************************
C
C
                      SUBROUTINE outterm(message,n)
C
C
C***********************************************************************
C
C
C     Subroutine to output a message to the terminal.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*(*) message
        INTEGER n

C=======================================================================
 
C   Write the message.
C   ------------------
        WRITE(uterm,*) message,n
 
        RETURN
        END
C***********************************************************************
C
C
                           SUBROUTINE setbox 
C
C
C***********************************************************************
C
C
C     Subroutine to adjust system box so that it contains all bodies.
C     The local variable rebox indicates whether a resizing of the
C     system box is to take place.  The variables posmin and posmax
C     are the minimum and maximum coordinates of bodies in each
C     dimension.  
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        DOUBLE PRECISION posmin(ndim),posmax(ndim),posx(nbodsmax),
     &                   posy(nbodsmax),posz(nbodsmax)
        INTEGER k,ismin,ismax
        LOGICAL rebox

        EQUIVALENCE (posx(1),pos(1,1)),(posy(1),pos(1,2)),
     &              (posz(1),pos(1,3))

        SAVE rebox

        DATA rebox/.TRUE./,rsize/0./,rmin/0.,0.,0./

C=======================================================================

C   Determine minimum and maximum coordinates of bodies.
C   ----------------------------------------------------
        posmin(1)=posx(ISMIN(nbodies,posx,1))
        posmin(2)=posy(ISMIN(nbodies,posy,1))
        posmin(3)=posz(ISMIN(nbodies,posz,1))
        posmax(1)=posx(ISMAX(nbodies,posx,1))
        posmax(2)=posy(ISMAX(nbodies,posy,1))
        posmax(3)=posz(ISMAX(nbodies,posz,1))

C   Determine if a resizing is required.
C   ------------------------------------
        DO 50 k=1,ndim
           IF(rmin(k).GT.posmin(k).OR.rmin(k)+rsize.LT.posmax(k)) 
     &        rebox=.TRUE.
 50     CONTINUE

C   If a resizing is necessary, recompute rsize and rmin.
C   -----------------------------------------------------

        IF(rebox) THEN

           DO 70 k=1,ndim
              rsize=MAX(rsize,posmax(k)-posmin(k))
 70        CONTINUE

           DO 80 k=1,ndim
              rmin(k)=0.5*(posmin(k)+posmax(k))-0.5*rsize
 80        CONTINUE

        ENDIF

        rebox=.FALSE.

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE startout
C
C
C***********************************************************************
C
C
C     Subroutine to initialize disk files for subsequent input/output.
C     All files, other than the binary body data file, are assumed to 
C     be of ASCII (text) form.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C=======================================================================
 
C   Open parameter file.
C   --------------------
        OPEN(UNIT=upars,FILE=parsfile,STATUS='OLD')
 
C   Open log file.
C   --------------
        OPEN(UNIT=ulog,FILE=logfile,STATUS='NEW')
 
C   Open input body data file.
C   --------------------------
        OPEN(UNIT=ubodsin,FILE=ibodfile,STATUS='UNKNOWN')

C   Open output binary body data file.
C   ----------------------------------
        OPEN(UNIT=ubodsout,FILE=obodfile,STATUS='NEW',
     &       FORM='UNFORMATTED')

C   Open output ascii body data file.
C   ---------------------------------
        OPEN(UNIT=ubodsasc,FILE=oascfile,STATUS='NEW')
 
        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE steppos
C
C
C***********************************************************************
C
C
C     Subroutine to advance the positions of the bodies for a
C     timestep dtime/2.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,k

C=======================================================================

C   Loop over all spatial coordinates for all bodies.
C   -------------------------------------------------

        DO 200 k=1,ndim
           DO 100 p=1,nbodies
              pos(p,k)=pos(p,k)+vel(p,k)*dtime2
 100       CONTINUE
 200    CONTINUE

C   Update position time, system time.
C   ----------------------------------
        tpos=tpos+dtime2
        tnow=tpos

        RETURN
        END

C***********************************************************************
C
C
                         SUBROUTINE stepsys(n)
C
C
C***********************************************************************
C
C
C     Subroutine to advance the state of the system by one large
C     timestep.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER n

C=======================================================================

C   Update positions by 1/2 step.
C   -----------------------------
           CALL steppos
C               -------

C   Zero out acceleration, compute acceleration, advance velocities.
C   ----------------------------------------------------------------

           CALL zeroacc
C               -------
           CALL gravity('acc ')
C               -------
           CALL stepvel
C               -------

C   Update positions by 1/2 step.
C   -----------------------------
           CALL steppos
C               -------

C   Output system state.
C   --------------------
        CALL outstate(n)
C            --------

        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE stepvel
C
C
C***********************************************************************
C
C
C     Subroutine to advance the velocities of the bodies for a
C     timestep dtime.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,k

C=======================================================================

C   Loop over all velocity components for all bodies.
C   -------------------------------------------------

        DO 200 k=1,ndim

CDIR$ IVDEP
           DO 100 p=1,nbodies
              vel(p,k)=vel(p,k)+acc(p,k)*dtime
 100       CONTINUE

 200    CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                           SUBROUTINE stopout
C
C
C***********************************************************************
C
C
C     Subroutine to close the output files.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C=======================================================================
 
C   Close the open files.
C   ---------------------
        CLOSE(UNIT=ubodsout)
        CLOSE(UNIT=ulog)
        CLOSE(UNIT=ubodsasc)
 
        RETURN
        END
C***********************************************************************
C
C
                       SUBROUTINE terror(message)
C
C
C***********************************************************************
C
C
C     Subroutine to terminate the program as the result of a fatal
C     error, close the output files, and dump timing information.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*(*) message

C=======================================================================

C   Write error message to the log file.
C   ------------------------------------
        CALL outerror(message)
C            --------

C   Output system state.
C   --------------------
        CALL corrpos('correct')
C            -------
        CALL zeropot
C            -------
        CALL gravity('pot ')
C            -------
        CALL outascib
C            --------
        CALL outbods
C            -------

C-----------------------------------------------------------------------
C   Stop timing, output timing data, close files, terminate the
C   simulation.
C-----------------------------------------------------------------------

        CALL SECOND(cputime1)

        CALL outcpu
C            ------
        CALL stopout
C            -------

        STOP
        END

C***********************************************************************
C
C
                     SUBROUTINE treewalk(p,nterms)
C
C
C***********************************************************************
C
C
C     Subroutine to walk through the tree and accumulate the list of
C     interactions for body p.  The interaction list is passed back
C     to the calling subroutine accgrav in the vector iterms, which
C     is equivalenced to the common array bodlist.  Vectorization is
C     achieved by processing all cells at the same level in the tree
C     simultaneously (Hernquist, J. Comput. Phys., submitted [1988]).
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        DOUBLE PRECISION cvmgt,testkeep
        INTEGER p,i,nnodes,nkeep,nsubdiv,nterms,iterms(nbodsmax),
     &          nodelist(nbodsmax),keepterm(nbodsmax)
        LOGICAL tolcrit

        EQUIVALENCE (iterms(1),bodlist(1)),(nodelist(1),celllist(1)),
     &              (keepterm(1),parent(1))

C=======================================================================

C   Initialize list of cells to examine.
C   ------------------------------------
        nterms=0
        nnodes=1
        nodelist(1)=root
     
 10     CONTINUE

C   Loop until no cells are left to examine.
C   ----------------------------------------
        IF(nnodes.GT.0) THEN

C   Apply tolerance criterion to list of cells.
C   -------------------------------------------

CDIR$ IVDEP
           DO 20 i=1,nnodes
              tolcrit=((pos(p,1)-pos(nodelist(i),1))**2+(pos(p,2)-
     &                pos(nodelist(i),2))**2+(pos(p,3)-
     &                pos(nodelist(i),3))**2).GE.
     &                (cellsize(nodelist(i))**2*tol2inv)
              testkeep=CVMGT(two,minustwo,tolcrit)
              keepterm(i)=INT(testkeep)
 20        CONTINUE

C-----------------------------------------------------------------------
C   Add cells which satisfy criterion to interaction list.  Note that,
C   depending on theta, self-interaction term will be included.
C-----------------------------------------------------------------------

           CALL WHENIGT(nnodes,keepterm,1,0,isubset,nkeep)

           IF(nterms+nkeep.GT.nbodsmax) 
     &        CALL terror(' array overflow in treewalk ')
C                  ------

CDIR$ IVDEP
           DO 30 i=1,nkeep
              iterms(nterms+i)=nodelist(isubset(i))
 30        CONTINUE

           nterms=nterms+nkeep

C-----------------------------------------------------------------------
C   Add subcells of cells which fail tolerance criterion to list of
C   cells to examine.
C-----------------------------------------------------------------------

           CALL WHENILT(nnodes,keepterm,1,0,isubset,nsubdiv)

           IF(8*nsubdiv.GT.nbodsmax.OR.8*nsubdiv.GT.ncells)
     &        CALL terror(' asubp overflow in treewalk ')
C                  ------

CDIR$ IVDEP
           DO 40 i=1,nsubdiv
              asubp(i)=subp(nodelist(isubset(i)),1)
              asubp(i+nsubdiv)=subp(nodelist(isubset(i)),2)
              asubp(i+2*nsubdiv)=subp(nodelist(isubset(i)),3)
              asubp(i+3*nsubdiv)=subp(nodelist(isubset(i)),4)
              asubp(i+4*nsubdiv)=subp(nodelist(isubset(i)),5)
              asubp(i+5*nsubdiv)=subp(nodelist(isubset(i)),6)
              asubp(i+6*nsubdiv)=subp(nodelist(isubset(i)),7)
              asubp(i+7*nsubdiv)=subp(nodelist(isubset(i)),8)
 40        CONTINUE

           CALL WHENNE(8*nsubdiv,asubp,1,0,isubset,nnodes)

CDIR$ IVDEP
           DO 60 i=1,nnodes
              nodelist(i)=asubp(isubset(i))
 60        CONTINUE

           GO TO 10

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE zeroacc
C
C
C***********************************************************************
C
C
C     Subroutine to zero out acceleration.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p

C=======================================================================

CDIR$ IVDEP
        DO 10 p=1,nbodies
           acc(p,1)=0.
           acc(p,2)=0.
           acc(p,3)=0.
 10     CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE zeropot
C
C
C***********************************************************************
C
C
C     Subroutine to zero out potential.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p

C=======================================================================

CDIR$ IVDEP
        DO 10 p=1,nbodies
           phi(p)=0.
 10     CONTINUE

        RETURN
        END

