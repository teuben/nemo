C***********************************************************************
C
C
C                            PROGRAM treecode
	SUBROUTINE treecode
C
C
C***********************************************************************
C
C
C     Fortran version of the Barnes/Hut hierarchical N-body code,
C     adapted from the original C version, for use on the MFE CRAY
C     system.  It is assumed that the dialect of Fortran in use is
C     recursive.  Otherwise, the code is written in standard
C     Fortran, with the proviso that force computations have been
C     optimized for a vector machine.
C
C     The computational system of units is determined by the input
C     data, with the assumption that G=1 .  Particles are not
C     required to have identical masses.
C
C     Version 2 allows quadrupole terms to be included in the
C     potential and acceleration calculations.
C
C     Version 2.1 contains copies of the force evaluation routines
C     that are optimized for either two or three dimensional
C     systems.
C
C     Version 2.2 contains copies of the force evaluation routines
C     that allow for arbitrarily large values of the tolerance
C     parameter, tol, by forcing subdivision of cells containing
C     the body for which to compute the acceleration components.
C
C        Version 2.2:  September 1, 1986   Lars Hernquist, U.C.B.
C
C
C=======================================================================
C
C
C     This is the top-level evolution program treecode.  Its tasks
C     are to:
C
C          1) initialize file structures and global variables;
C          2) input parameters and the initial system state;
C          3) advance the state of the system for a given number
C             of timesteps;
C          4) perform a diagnostic analysis of the system at
C             each time step (energy, angular momentum, etc.);
C          5) periodically record the state of the system;
C          6) terminate the simulation and close data files.
C
C
C=======================================================================
C
C
C     Global variables/parameters:
C
C          ndim        : the number of spatial dimensions.
C          nsubcell    : the number of subcells per cell.
C          nbodsmax    : the maximum number of allowed bodies.
C          nbodies     : the number of bodies in the system.
C          ncells      : the maximum number of allowed cells.
C          maxnterm    : maximum number of terms allowed in
C                        force evaluations.
C          nbits       : number of bits in an integer word.
C          null        : alternate name for zero (0).
C          incells     : number of cells currently in use.
C          imax        : highest significant bit in an integer
C                        word, used to create body/cell tree.
C          imax2       : imax, right bit-shifted by one bit.
C          nindex      : vector used by subindex to create tree.
C          root        : pointer to the top of the tree.
C          rmin        : cartesian coordinates of lower-left
C                        corner of the system box.
C          rsize       : length of the system box.
C          tol         : accuracy parameter.
C          eps         : potential softening parameter.
C          eps2        : square of eps.
C          epsinv      : 1. / eps .
C          headline    : identification string for the run.
C          usequad     : option to include (.TRUE.) quadrupole
C                        terms.
C
C          mass        : masses of bodies, total masses of cells.
C          pos         : coordinates of bodies, center of mass
C                        coordinates of a cell.
C          vel         : velocity components of a body.
C          acc         : acceleration components of a body.
C          phi         : potential field at a body.
C          subp        : pointers to descendents of a cell.
C          quad        : quadrupole moments of cells.
C          sizetol2    : square of ((size of a body/cell) / tol).
C
C          pskip       : pointer to body on which to compute force.
C          pos0        : coordinates of body pskip.
C          nterms      : number of terms in force evaluation of
C                        body pskip.
C          nttot       : running sum of nterms over all bodies in
C                        a force evaluation.
C          ntmin       : minimum nterms for a body in a given step.
C          ntmax       : maximum nterms for a body in a given step.
C          ntavg       : average number of nterms/body in a given step.
C          pmass       : masses of bodies/cells in force evaluation of
C                        body pskip.
C          dr          : displacement vector from bodies/cells to
C                        body pskip.
C          drdotdr     : the dot product dr * dr.
C
C          nqterms     : number of quadrupole terms in force
C                        evaluation of body pskip.
C          qdrdotdr,   : temporary storage for
C          qquad,qdr     quadrupole calculation.
C
C          nsteps      : total number of timesteps in the simulation.
C          nout        : frequency of system record outputs.
C          tnow        : current system time.
C          tpos        : current position time.
C          tvel        : current velocity time.
C          dtime       : the timestep.
C          dtime2      : timestep/2.
C
C          uterm, upars, ubods, ulog   : logical i/o unit numbers.
C          parsfile, bodsfile, logfile : character names of files.
C
C          cputime     : cpu time (secs) used during the simulation.
C          cputime0    : cumulative cpu time at start of run.
C          cputime1    : cumulative cpu time at end of run.
C
C
C=======================================================================
C
C
C     Data structure:
C
C          The body/cell tree structure is assumed to be of the
C          form discussed by Barnes and Hut.  Schematically, for
C          three dimensions (i.e. eight subcells per cell):
C
C         +-----------------------------------------------------------+
C  root-->| CELL:  mass, pos, quad, sizetol2, /, o, /, /, /, /, o, /  |
C         +--------------------------------------|--------------|-----+
C                                                |              |
C     +------------------------------------------+              |
C     |                                                         |
C     |   +--------------------------------------------+        |
C     +-->| BODY:  mass, pos, vel, acc, phi, sizetol2  |        |
C         +--------------------------------------------+        |
C                                                               |
C     +---------------------------------------------------------+
C     |
C     |   +-----------------------------------------------------------+
C     +-->| CELL:  mass, pos, quad, sizetol2, o, /, /, o, /, /, o, /  |
C         +-----------------------------------|--------|--------|-----+
C                                             |        |        |
C                                            etc.     etc.     etc.
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
C=======================================================================

C-----------------------------------------------------------------------
C   Definition of macro for including the Cliche file--specific to
C   the CTSS operating system.
C-----------------------------------------------------------------------

        include 'treedefs.inc'

C-----------------------------------------------------------------------
C   Create a dropfile, associate logical i/o unit 6 with terminal.
C-----------------------------------------------------------------------

CPJT        CALL LINK("unit6=tty//")

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
C-----------------------------------------------------------------------
C   Generate tree structure and compute acceleration components.
C-----------------------------------------------------------------------

        CALL maketree
C            --------
        CALL accel
C            -----
C-----------------------------------------------------------------------
C   Output header for log file, initial force evaluation and system
C   diagnostics.
C-----------------------------------------------------------------------
        CALL outlog(0)
C            ------
        CALL energy
C            ------
C-----------------------------------------------------------------------
C   Advance body velocities by 1/2 timestep, positions by 1 timestep.
C-----------------------------------------------------------------------

        CALL stepvel(dtime2)
C            -------
        CALL steppos(dtime)
C            -------
C=======================================================================
C   Primary loop to advance system state for a given number of steps.
C=======================================================================

        DO 100 n=1,nsteps

C   Write message to terminal.
C   --------------------------
           CALL outterm(' istep= ',n)
C               -------

C   Create tree, compute forces.
C   ----------------------------
           CALL maketree
C               --------
           CALL accel
C               -----

C   Output force evaluation diagnostics.
C   ------------------------------------
           CALL outlog(n)
C               ------

C   Advance velocity 1/2 step to synchronize tvel, tpos.
C   ----------------------------------------------------
           CALL stepvel(dtime2)
C               -------

C   Output system diagnostics and, if appropriate, system state.
C   ------------------------------------------------------------
           CALL energy
C               ------
           IF(MOD(n,nout).EQ.0) CALL outbods
C                                    -------

C   Advance velocity additional 1/2 step.
C   -------------------------------------
           CALL stepvel(dtime2)
C               -------

C   Advance position 1 full step.
C   -----------------------------
           CALL steppos(dtime)
C               -------

 100    CONTINUE

C-----------------------------------------------------------------------
C   Stop timing, write timing data, close files, terminate simulation.
C-----------------------------------------------------------------------

        CALL SECOND(cputime1)

        CALL outcpu
C            ------
        CALL stopout
C            -------

        CALL EXIT(0)
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
C     either the input data or defined PARAMETERS.  The local
C     variable p is a pointer to the bodies/cells.
C
C
C=======================================================================

C-----------------------------------------------------------------------
C   Standard CTSS macro for including the Cliche file.
C-----------------------------------------------------------------------

        include 'treedefs.inc'

        INTEGER p

C   Initialize position and velocity times, 1/2 timestep.
C   -----------------------------------------------------
        tpos=tnow
        tvel=tnow
        dtime2=.5*dtime

C   Compute square, inverse of the softening parameter.
C   ---------------------------------------------------
        eps2=eps*eps
        epsinv=1./eps

C   Initialize imax, imax2, nindex for creating the tree.
C   -----------------------------------------------------

        imax=1

        DO 10 i=1,nbits-2
           imax=2*imax
 10     CONTINUE

        imax2=imax/2


        nindex(1)=2**(ndim-1)

        DO 20 k=2,ndim
           nindex(k)=nindex(k-1)/2
 20     CONTINUE

C=======================================================================
C   Compute initial box properties from input body data.
C=======================================================================

        rtemp=0.

        DO 40 k=1,ndim

           DO 30 p=1,nbodies
              if(abs(pos(p,k)).GT.rtemp) rtemp=abs(pos(p,k))
 30        CONTINUE

 40     CONTINUE

        DO 50 k=1,ndim
           rmin(k)=-2.*rtemp
 50     CONTINUE

        rsize=-2.*rmin(1)

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
C     Initialize the tree structure for the force/potential
C     calculation.  The local variable p is a pointer to the
C     bodies.
C
C
C=======================================================================

C-----------------------------------------------------------------------
C   Standard CTSS macro for including the Cliche file.
C-----------------------------------------------------------------------

        include 'treedefs.inc'

        INTEGER p

C   Deallocate current tree.
C   ------------------------
        root=null
        incells=0

C=======================================================================
C   Load bodies into the new tree, one at a time.
C=======================================================================

        DO 100 p=1,nbodies
           CALL hackload(p)
C               --------

 100    CONTINUE

C-----------------------------------------------------------------------
C   Compute masses, center of mass coordinates, and
C   (sizes of bodies/cells / tol) **2 for the tree.
C-----------------------------------------------------------------------
        CALL hackcofm(root,rsize*rsize/(tol*tol))
C            --------

C   Compute quadrupole moments of cells, if required.
C   -------------------------------------------------
        IF(usequad) CALL hackquad(root)
C                        --------

        RETURN
        END

C***********************************************************************
C
C
                           SUBROUTINE accel
C
C
C***********************************************************************
C
C
C     Subroutine to compute the acceleration components and potential
C     for all of the bodies.  The local variable p is a pointer to
C     the bodies.
C
C
C=======================================================================

C-----------------------------------------------------------------------
C   Standard CTSS macro for including the Cliche file.
C-----------------------------------------------------------------------

        include 'treedefs.inc'

        INTEGER p

C   Initialize the force evaluation diagnostics.
C   --------------------------------------------
        nttot=0
        ntmin=maxnterm
        ntmax=0

C=======================================================================
C   Loop over bodies--compute acceleration for body labeled pskip.
C=======================================================================

        DO 200 p=1,nbodies
           pskip=p

C   Where to evaluate force, potential.
C   -----------------------------------
           DO 100 k=1,ndim
              pos0(k)=pos(p,k)
 100       CONTINUE

C   Perform force computation for body pskip.
C   -----------------------------------------
           CALL hackgrav
C               --------

C   Update diagnostics, avoiding self-force term.
C   ---------------------------------------------
           nttot=nttot+nterms-1
           IF(nterms-1.LT.ntmin) ntmin=nterms-1
           IF(nterms-1.GT.ntmax) ntmax=nterms-1

 200    CONTINUE

C   Compute average number of force terms per body.
C   -----------------------------------------------
        ntavg=nttot/nbodies

        RETURN
        END

C***********************************************************************
C
C
                         SUBROUTINE stepvel(dt)
C
C
C***********************************************************************
C
C
C     Subroutine to advance the velocities of the bodies for a
C     timestep dt.  The argument dt allows the velocities to be
C     stepped for either 1/2 timestep (initially, or to synchronize
C     tvel and tpos) or a full timestep.  The local variable p is
C     a pointer to the bodies.
C
C
C=======================================================================

C-----------------------------------------------------------------------
C   Standard CTSS macro for including the Cliche file.
C-----------------------------------------------------------------------

        include 'treedefs.inc'

        INTEGER p

C=======================================================================
C   Loop over all velocity components for all bodies.
C=======================================================================

        DO 200 k=1,ndim
           DO 100 p=1,nbodies
              vel(p,k)=vel(p,k)+acc(p,k)*dt
 100       CONTINUE
 200    CONTINUE

C   Update velocity time, system time.
C   ----------------------------------
        tvel=tvel+dt
        tnow=tvel

        RETURN
        END

C***********************************************************************
C
C
                         SUBROUTINE steppos(dt)
C
C
C***********************************************************************
C
C
C     Subroutine to advance the positions of the bodies for a
C     timestep dt.  The local variable p is a pointer to the bodies.
C
C
C=======================================================================

C-----------------------------------------------------------------------
C   Standard CTSS macro for including the Cliche file.
C-----------------------------------------------------------------------

        include 'treedefs.inc'

        INTEGER p

C=======================================================================
C   Loop over all spatial coordinates for all bodies.
C=======================================================================

        DO 200 k=1,ndim
           DO 100 p=1,nbodies
              pos(p,k)=pos(p,k)+vel(p,k)*dt
 100       CONTINUE
 200    CONTINUE

C   Update position time, system time.
C   ----------------------------------
        tpos=tpos+dt
        tnow=tpos

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
C     local variable p is a pointer to the bodies.  The local
C     variables mtot, etot, ektot, eptot, cmpos, cmvel, and amvec are
C     the total mass, total energy, total kinetic energy, total
C     potential energy, center of mass position and velocity, and
C     total angular momentum of the system, respectively.
C
C
C=======================================================================

C-----------------------------------------------------------------------
C   Standard CTSS macro for including the Cliche file.
C-----------------------------------------------------------------------

        include 'treedefs.inc'

        INTEGER p
        REAL mtot,cmpos(ndim),cmvel(ndim),amvec(3)

C-----------------------------------------------------------------------
C   Terminate the simulation if tvel, tpos not synchronized.
C-----------------------------------------------------------------------

        IF(ABS(tpos-tvel).GT.1.e-3*dtime2)
     &     CALL terror(' error in energy--times not synchronized')
C               ------

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
        CALL outenrgy(mtot,etot,ektot,eptot,amvec,cmpos,cmvel)
C            --------

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

        CHARACTER*(*) message

C   Write error message to the log file.
C   ------------------------------------
        CALL outerror(message)
C            --------

C-----------------------------------------------------------------------
C   Stop timing, output timing data, close files, terminate the
C   simulation.
C-----------------------------------------------------------------------

        CALL SECOND(cputime1)

        CALL outcpu
C            ------
        CALL stopout
C            -------

CPJT        CALL EXIT(0)
        END


