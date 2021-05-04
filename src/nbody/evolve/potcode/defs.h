/*
 * DEFS.H: parameter and structure definitions for potential code.
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <potential.h>

/*
 * GLOBAL: pseudo-keyword for storage class.
 */
 
#if !defined(global)
#  define global extern
#endif

global string infile;		/* input file with initial conds */
global string outfile;		/* output file for simulation results */
global string savefile;		/* output file for system state */

global real freq;		/* fundamental integration frequency */

global int mode;		/* integrator: RK, PC or PC1 */

global real eta;		/* fractional dissipation [0..1] */
global real sigma;              /* diffusion angle (in radians) */

global vector dr;		/* cell size for dissipation */
global real rmax;               /* max. gridsize for dissipation */

global real freqout, minor_freqout;	/* major, minor output frequencies */

global real tstop;		/* time to stop integration */

global string options;		/* misc. options */

global string headline;		/* identification message */

global real tnow;		/* time state is defined */

global real tout, minor_tout;	/* time of next major, minor output */

global real ome, ome2, half_ome2, two_ome;	/* pattern speed + handy numbers */



/*
 * BODY: structure storing fundamental per-particle variables.
 */

typedef struct {
    real mass;			/* particle mass */
    vector phase[2];		/* phase-space coordinates */
    real phi;			/* potential at position */
    vector acc;			/* gravitational acceleration */
#if 1
    real A,B,kappa,nu,xiv0,etav0,zetav0;    /* for now: the epi constants */
#endif
    int key;                    /* some index */
} body, *bodyptr;

#define Body	 body
#define Mass(p)  ((p)->mass)
#define Phase(p) ((p)->phase)
#define Pos(p)   (Phase(p)[0])
#define Vel(p)   (Phase(p)[1])
#define Acc(p)   ((p)->acc)
#define Phi(p)   ((p)->phi)
#define Key(p)   ((p)->key)

#ifndef MBODY
#  define MBODY 250000
#endif

global int nbody;		/* number of bodies simulated */

global body bodytab[MBODY];	/* array representing state */


extern int diffuse(body *btab, int nb, int ndim, real sigma);
extern int dissipate(body *btab, int nb, int ndim, real *dr, real eta, real grid);
/* orbstep.c */
void initstep(bodyptr btab, int nb, real *tptr, proc force);
void orbstep(bodyptr btab, int nb, real *tptr, proc force, real dt, int mode);
void rkstep(bodyptr btab, int nb, real *tptr, proc force, real dt, real atmp1[]);
void pcstep(bodyptr btab, int nb, real *tptr, proc force, real dt);
void eulerstep(bodyptr btab, int nb, real *tptr, proc force, real dt);
void modeulerstep(bodyptr btab, int nb, real *tptr, proc force, real dt);
void leapfrogstep(bodyptr btab, int nb, real *tptr, proc force, real dt);
void rk4step(bodyptr btab, int nb, real *tptr, proc force, real dt);
void epistep(bodyptr btab, int nb, real *tptr, proc force, real dt, int mode);
/* code_io.c */
void inputdata(void);
void initoutput(void);
void stopoutput(void);
void output(void);
void savestate(string file);
void restorestate(string file);
