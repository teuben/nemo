/*
 * CODE.H: define various global things for CODE.C.
 */

#include "defs.h"
#include "proto.h"
#include <getparam.h>


/* from defs.h */

/*
 * ROOT: origin of tree; declared as nodeptr for tree with only 1 body.
 */

global nodeptr troot;

/*
 * Integerized coordinates: used to mantain body-tree.
 */

global vector rmin;			/* lower-left corner of coord. box */
global real rsize;			/* side-length of int. coord. box */

/*
 * Parameters and results for gravitational calculation.
 */

global real fcells;			/* ratio of cells/bodies allocated */

global real tol;                       /* accuracy parameter: 0.0 => exact */
global real eps;                       /* potential softening parameter */

global int n2bterm;                    /* number 2-body of terms evaluated */
global int nbcterm;			/* num of body-cell terms evaluated */

global bool debug;                     /* control debugging messages */


/* old */

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

global int nfcalc;			/* number of n-on-1 force calculations */
global int n2bcalc;			/* number of 2-body force calculations */
global int nbccalc;			/* num of body-cell force calculations */

global int nbody;			/* number of bodies in system */

global bodyptr bodytab;		/* points to array of bodies */
