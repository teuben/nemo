      SUBROUTINE DEFINE
*
*
*       Definition of input parameters, options & counters.
*       ---------------------------------------------------
*
*
*       Input parameters
*       ****************
*
*       ---------------------------------------------------------------------
*       KSTART  Control index (1: new run; >1: restart; 3, 4, 5: new params).
*       TCOMP   Maximum computing time in minutes (saved in CPU).
*
*       N       Total particle number (<= NMAX).
*       NFIX    Output frequency of data save or binaries (option 3 & 6).
*       NRAND   Random number sequence skip.
*       NNBMAX  Maximum number of neighbours (< LMAX).
*       NRUN    Run identification index.
*
*       ETAI    Time-step parameter for irregular force polynomial.
*       ETAR    Time-step parameter for regular force polynomial.
*       RS0     Initial radius of neighbour sphere.
*       DELTAT  Output time interval in units of the crossing time.
*       TCRIT   Termination time in units of the crossing time.
*       QE      Energy tolerance (restart if DE/E > 5*QE & KZ(2) > 1).
*       EPS     Softening parameter (square saved in EPS2).
*
*       KZ(J)   Non-zero options for alternative paths (see table).
*
*       XTPAR1  Mass of external Plummer model (KZ(15) = 1; scaled units).
*       XTPAR2  Length scale for Plummer model (KZ(15) = 1).
*       ZMGAS   Mass scale for logarithmic potential (KZ(15) = 2).
*       RGAS    Length scale for logarithmic potential (KZ(15) = 2).
*
*       ALPHAS  Power-law index for initial mass function (not KZ(4) = 2)
*       BODY1   Maximum particle mass before scaling.
*       BODYN   Minimum particle mass before scaling.
*
*       Q       Virial ratio (routine SCALE; Q = 0.5 for equilibrium).
*       VXROT   XY-velocity scaling factor (> 0 for solid-body rotation).
*       VZROT   Z-velocity scaling factor (not used if VXROT = 0).
*       RBAR    Virial radius in pc (for scaling to physical units).
*       ZMBAR   Mean mass in solar units.
*
*       XCM     Displacement for subsystem (routine SCALE; option 17).
*       ECC     Eccentricity of relative motion for subsystem (ECC =< 1).
*       ---------------------------------------------------------------------
*
*
*       Options KZ(J)
*       *************
*
*       ---------------------------------------------------------------------
*       1  COMMON save on unit 1 at end of run (=2: every 100*NMAX steps).
*       2  COMMON save on unit 2 at output (=1); restart if DE/E > 5*QE (=2).
*       3  Basic data written to unit 3 at output time (frequency NFIX).
*       4  Initial conditions on unit 4 (=1: output; =2: input).
*       5  Initial conditions (=0: uniform & isotropic; =1: Plummer model).
*       6  Output of significant binaries (=2: frequency NFIX).
*       7  Lagrangian radii (=1: unit 6; =2: unit 7; =3: both types).
*       8  Core radius & density centre (N > 20 only).
*       9  Individual bodies printed at output time (MIN(5**KZ9,N)).
*      10  No unique density centre (skips velocity modification of RS(I)).
*      11  Modification of ETAI & ETAR by tolerance QE.
*      12  Inelastic mergers (>1: diagnostic output).
*      13  Escaper removal (R > 2*RTIDE; RTIDE = 10*RSCALE).
*      14  Skip full predictor loop if NNB > KZ(14) = <NNB> & KZ(14) > 0.
*      15  External potential (=1: Plummer model; =2: logarithmic potential).
*      16  No scaling of initial conditions.
*      17  Generation of two subsystems (merger experiment).
*      18  Adjustment of coordinates & velocities to c.m. condition.
*      19  Not used at present (same for # 20).
*       ---------------------------------------------------------------------
*
*
*       Output counters
*       ***************
*
*       -------------------------------------------------------------------
*       NSTEPI  Irregular integration steps.
*       NSTEPR  Regular integration steps.
*       NNPRED  Coordinate predictions of all particles (cf. option 14).
*       NBCORR  Force polynomial corrections.
*       NBFULL  Too many neighbours with standard criterion.
*       NBVOID  No neighbours inside 1.26 times the basic sphere radius.
*       NMTRY   Merger attempts (option 12).
*       NMERG   Mergers (option 12).
*       NESC    Escaped particles (option 13).
*       -------------------------------------------------------------------
*
*
      RETURN
*
      END
