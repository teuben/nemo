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

int nfcalc;			/* number of n-on-1 force calculations */
int n2bcalc;			/* number of 2-body force calculations */
int nbccalc;			/* num of body-cell force calculations */

int nbody;			/* number of bodies in system */
int nrigid;                     /* number of rigid particles */
bodyptr bodytab;		/* points to array of bodies */

proc extpot;			/* external potential(5) */
