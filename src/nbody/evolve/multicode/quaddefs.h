/*
 * QUADDEFS.H: definitions and global variables for quadrupole codes.
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <snapshot/body.h>

#include "quadfield.h"

typedef void (*force_proc)(Body *, int , real);

string infile;			/* input file with initial conds            */
string outfile;			/* output file for simulation results       */
string quadfile;		/* output file for field tables             */
string savefile;		/* output file for system state             */

real freq;			/* fundamental integration frequency        */

int mode;			/* integrator: RK, PC or PC1                */

real eps1, eps2;		/* radial, tangential softening             */

real freqout, minor_freqout;	/* major, minor output frequencies          */

real tstop;			/* time to stop integration                 */

string options;			/* misc. options                            */

string headline;		/* identification message                   */

real tnow;			/* time state is defined                    */

real tout, minor_tout;		/* time of next major, minor output         */

quadfield qfld;			/* tables of field moments                  */

int nbody;			/* number of bodies simulated               */

Body *bodytab;			/* array representing state                 */

#if !defined(MBODY)
#  define MBODY 100000		/* max number of bodies, for orbstep        */
#endif

/* quadcode_io.c */
extern void inputdata(void);
extern void initoutput(void);
extern void stopoutput(void);
extern void output(void);
extern int savestate(string file);
extern int restorestate(string file);

/* orbstep.c */
extern int initstep(body *btab, int nb, real *tptr, force_proc force);
extern int orbstep(body *btab, int nb, real *tptr, force_proc force, real dt, int mode);
extern int rkstep(body *btab, int nb, real *tptr, force_proc force, real dt, real atmp1[]);
extern int pcstep(body *btab, int nb, real *tptr, force_proc force, real dt);
extern int moveaccel(body *btab, int nb);

/* quadforce.c */
extern int quadforce(body *btab, int nb, real eps1, real eps2);

