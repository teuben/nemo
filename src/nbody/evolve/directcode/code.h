/*
 * CODE.H: define various global things for CODE.C.
 */

#include "defs.h"
#include <getparam.h>

global string infile;			/* file name for snapshot input */
global string outfile;			/* file name for snapshot output */
global string savefile;		/* file name for state output */

global real freq;			/* fundamental integration frequency */

global real freqout, minor_freqout;	/* major, minor output frequencies */

global real tstop;			/* time to stop calculation */

global string options;                 /* various option flags */

extern string headline;		/* message describing calculation */

global real tnow;			/* current value of time */

global real tout, minor_tout;		/* time of next major, minor output */

global int nstep;			/* number of micro-steps */

global int nbody;			/* number of bodies in system */

global bodyptr bodytab;		/* points to array of bodies */

global real eps;                       /* grav softening length */

global real gravc;                     /* gravitational constant [1] */

/* code.c */
void nemo_main(void);
void startrun(void);
void testdata(bool cencon);

void stepsystem_leapfrog(void);
void stepsystem_euler(void);
void stepsystem_old(void);

/* code_io.c */
void inputdata(string file);
void initoutput(void);
void stopoutput(void);
void output(void);

/* util.c */
void pickvec(vector x, bool cf);

/* grav.c */
void hackgrav(bodyptr p);
