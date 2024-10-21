/*
 * QUADDEFS.H: definitions and global variables for quadrupole codes.
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <snapshot/body.h>

#include "quadfield.h"
/*
 * GLOBAL: pseudo-keyword for storage class.
 */
 
#if !defined(global)
#  define global extern
#endif

typedef void (*force_proc)(Body *, int , real);

global string infile;			/* input file with initial conds            */
global string outfile;			/* output file for simulation results       */
global string quadfile;		/* output file for field tables             */
global string savefile;		/* output file for system state             */

global real freq;			/* fundamental integration frequency        */

global int mode;			/* integrator: RK, PC or PC1                */

global real eps1, eps2;		/* radial, tangential softening             */

global real freqout, minor_freqout;	/* major, minor output frequencies          */

global real tstop;			/* time to stop integration                 */

global string options;			/* misc. options                            */

global string headline;		/* identification message                   */

global real tnow;			/* time state is defined                    */

global real tout, minor_tout;		/* time of next major, minor output         */

global quadfield qfld;			/* tables of field moments                  */

global int nbody;			/* number of bodies simulated               */

global Body *bodytab;			/* array representing state                 */

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
extern void initstep(body *btab, int nb, real *tptr, force_proc force);
extern void orbstep(body *btab, int nb, real *tptr, force_proc force, real dt, int mode);
extern void rkstep(body *btab, int nb, real *tptr, force_proc force, real dt, real atmp1[]);
extern void pcstep(body *btab, int nb, real *tptr, force_proc force, real dt);
extern void moveaccel(body *btab, int nb);

/* quadforce.c */
extern void quadforce(body *btab, int nb, real eps1, real eps2);

/* quadinter.c */
extern void quadinter(Body *btab, int nb, real eps1, real eps2);
