/*
 * DEFS.H: parameter and structure definitions for potential code.
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <potential.h>

string infile;			/* input file with initial conds */
string outfile;			/* output file for simulation results */

real freq;			/* fundamental integration frequency */

int mode;			/* integrator: RK, PC or PC1 */

real eta;			/* fractional dissipation [0..1] */
real sigma;                     /* diffusion angle (in radians) */
real freqdiff;                  /* diffusion frequency */
real fheat;                     /* diffusion/dissipation */
vector dr;			/* cell size for dissipation */
real rmax;                      /* max. gridsize for dissipation */

real freqout, minor_freqout;	/* major, minor output frequencies */

real tstop;			/* time to stop integration */

string options;			/* misc. options */

string headline;		/* identification message */

real tnow;			/* time state is defined */

real tout, minor_tout;		/* time of next major, minor output */

real ome, ome2, half_ome2, two_ome;	/* pattern speed + handy numbers */



/*
 * BODY: structure storing fundamental per-particle variables.
 */

typedef struct {
    real mass;			/* particle mass */
    vector phase[2];		/* phase-space coordinates */
    real phi;			/* potential at position */
    vector acc;			/* gravitational acceleration */
    int key;                    /* some index */
    real aux;                   /* auxiliary data */
} body, *bodyptr;

#define Body	 body
#define Mass(p)  ((p)->mass)
#define Phase(p) ((p)->phase)
#define Pos(p)   (Phase(p)[0])
#define Vel(p)   (Phase(p)[1])
#define Acc(p)   ((p)->acc)
#define Phi(p)   ((p)->phi)
#define Aux(p)   ((p)->aux)
#define Key(p)   ((p)->key)

#ifndef MBODY
#  define MBODY 4096		/* max number of bodies */
#endif

int nbody;			/* number of bodies simulated */

body bodytab[MBODY];		/* array representing state */

typedef void (*fproc)(bodyptr p, int np, real time, bool Qnew);

/* flowcode.c */

extern void setparams(void);
extern void force(bodyptr btab, int nb, real time, bool Qnew);

/* diffuse.c */
extern void rotate_aux(bodyptr btab, int nb);
extern void diffuse(body *btab, int nb, int ndim, real sigma, bool Qrotate);

/* code_io.c */

extern void inputdata(void);
extern void initoutput(void);
extern void stopoutput(void);
extern void output(void);
extern void savestate(string file);
extern void restorestate(string file);

/* orbstep.c */
extern void initstep(bodyptr btab, int nb, real *tptr, fproc force);
extern void orbstep(bodyptr btab, int nb, real *tptr, fproc force, real dt, int mode);
