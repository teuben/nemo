      SUBROUTINE DEFINE
*
*
*       Definition of input parameters & options.
*       -----------------------------------------
*
*
*       Input parameters
*       ****************
*
*       ---------------------------------------------------------------------
*       KSTART  Control index (1: new run; >1: restart; 3: new params).
*       TCOMP   Maximum computing time in minutes (saved in CPU).
*
*       N       Total particle number.
*       NFIX    Output frequency of data save or binaries (option 3 & 6).
*       NRAND   Random number sequence skip.
*       NRUN    Run identification index.
*
*       ETA     Time-step parameter for total force polynomial.
*       DELTAT  Output time interval in units of the crossing time.
*       TCRIT   Termination time in units of the crossing time.
*       QE      Energy tolerance (stop if DE/E > 5*QE & KZ(2) <= 1).
*       EPS     Softening parameter (square saved in EPS2).
*
*       KZ(J)   Non-zero options for alternative paths (see table).
*
*       ALPHAS  Power-law index for initial mass function (routine DATA).
*       BODY1   Maximum particle mass before scaling.
*       BODYN   Minimum particle mass before scaling.
*
*       Q       Virial ratio (routine SCALE; Q = 0.5 for equilibrium).
*       VXROT   XY-velocity scaling factor (> 0 for solid-body rotation).
*       VZROT   Z-velocity scaling factor (not used if VXROT = 0).
*       RBAR    Virial radius in pc (for scaling to physical units).
*       ZMBAR   Mean mass in solar units.
*
*       NFRAME  Maximum number of movie frames (routine SCALE; option 7).
*       DELTAF  Time interval for movie output (in units of crossing time).
*
*       XCM     Displacement for subsystem (routine SCALE; option 8).
*       ECC     Eccentricity of relative motion for subsystem (ECC =< 1).
*
*       SEMI    Semi-major axis of binary orbit (options 10 & 12).
*       ECC     Eccentricity of binary.
*       ---------------------------------------------------------------------
*
*
*       Options KZ(J)
*       *************
*
*       ---------------------------------------------------------------------
*       1  COMMON save on unit 1 if TCOMP > CPU or if TIME > TCRIT.
*       2  COMMON save on unit 2 at output (=1); restart if DE/E > 5*QE (=2).
*       3  Basic data written to unit 3 at output time (frequency NFIX).
*       4  Initial conditions on unit 4 (=1: output; =2: input).
*       5  Initial conditions (=0: uniform & isotropic; =1: Plummer).
*       6  Output of significant binaries.
*       7  Output of movie frames on unit 7.
*       8  Generation of two subsystems (merger experiment).
*       9  Individual bodies printed at output time (MIN(5**KZ9,N)).
*      10  No scaling of initial conditions.
*      11  Modification of ETA by tolerance QE.
*      12  Initial parameters for binary orbit.
*      13  Escaper removal (R > 2*RTIDE; RTIDE = 10*RSCALE).
*      14  Adjustment of coordinates & velocities to c.m. condition.
*       ---------------------------------------------------------------------
*
*
      RETURN
*
      END
