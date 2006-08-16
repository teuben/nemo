/*
 * QUADCODE.C: main routines for global quadrupole N-body code.
 *
 *	 3-mar-89  V1.2  JEB	Some last version
 *	12-nov-91  V1.3  PJT  Nemo V2.
 *      20-may-94  V1.3a pjt  usage
 */


#include "quaddefs.h"

string defv[] = {		/* DEFAULT INPUT PARAMETERS		    */
    "in=???\n			input file name",
    "out=\n			output file name",
    "quad=\n			field file name",
    "save=\n			state file name",
    "eps_r=0.05\n		radial softening parameter",
    "eps_t=0.07\n		tangential softening parameter",
    "freq=64.0\n		fundamental frequency (inv delta-t)",
    "mode=3\n			integrator: 1 => RK, 2 => PC, 3 => PC1",
    "tstop=2.0\n		time to stop integration",
    "freqout=4.0\n		major data-output frequency",
    "minor_freqout=32.0\n	minor data-output frequency",
    "options=\n			misc options",
    "headline=\n		random mumble for humans",
    "VERSION=1.3b\n	        15-aug-06 PJT",
    NULL,
};

string usage = "Global quadrupole N-body integrator";

string cvsid="$Id$";

local void force(Body *, int , real);


void nemo_main(void)
{
    setparams();
    inputdata();
    initoutput();
    initstep(bodytab, nbody, &tnow, force);
    output();
    while (tnow + 0.1/freq < tstop) {
	orbstep(bodytab, nbody, &tnow, force, 1.0/freq, mode);
	output();
    }
    stopoutput();
}

setparams()
{
    infile = getparam("in");
    outfile = getparam("out");
    quadfile = getparam("quad");
    savefile = getparam("save");
    eps1 = getdparam("eps_r");
    eps2 = getdparam("eps_t");
    freq = getdparam("freq");
    mode = getiparam("mode");
    tstop = getdparam("tstop");
    freqout = getdparam("freqout");
    minor_freqout = getdparam("minor_freqout");
    options = getparam("options");
    headline = getparam("headline");
}

/*
 * FORCE: force calcualtion routine.
 */

local void force(
	   Body *btab,			/* array of bodies */
	   int nb,			/* number of bodies */
	   real time)			/* current time (ignored) */
{
    Body *p;

    for (p = btab; p < btab+nb; p++) {		/* loop over bodies         */
	CLRV(Acc(p));				/*   zero acceleration      */
	Phi(p) = 0.0;				/*   and potential          */
    }
    qfld.nqtab = 0;				/* get new field points     */
    quadforce(btab, nb, eps1, eps2);		/* compute quadpole force   */
}
