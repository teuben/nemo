/*
 * QUADDEFS.H: definitions and global variables for quadrupole codes.
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <snapshot/body.h>

#include "quadfield.h"

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
#  define MBODY 4096		/* max number of bodies, for orbstep        */
#endif
