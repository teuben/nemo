/*
 * CODE.H: define various global things for CODE.C.
 */

#include "defs.h"
#include <getparam.h>

string infile;			/* file name for snapshot input */
string outfile;			/* file name for snapshot output */
string savefile;		/* file name for state output */

real freq;			/* fundamental integration frequency */

real freqout, minor_freqout;	/* major, minor output frequencies */

real tstop;			/* time to stop calculation */

string options;                 /* various option flags */

extern string headline;		/* message describing calculation */

real tnow;			/* current value of time */

real tout, minor_tout;		/* time of next major, minor output */

int nstep;			/* number of micro-steps */

int nbody;			/* number of bodies in system */

bodyptr bodytab;		/* points to array of bodies */

real eps;                       /* grav softening length */

/* code.c */
void nemo_main(void);
void startrun(void);
void testdata(bool cencon);
void stepsystem(void);

/* code_io.c */
void inputdata(string file);
void initoutput(void);
void stopoutput(void);
void output(void);

/* util.c */
void pickvec(vector x, bool cf);

/* grav.c */
void hackgrav(bodyptr p);
