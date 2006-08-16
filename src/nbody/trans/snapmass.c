/*
 *  SNAPMASS:  set masses in a snapshot
 *
 *	1-Aug-88	V1.0 created      				PJT
 *	5-dec-88	V1.1 allow all masses same through one keyword  PJT
 *				(proposed - to be coded )
 *	5-feb-89	V1.2 allow bodytrans(3NEMO) expression for mass PJT
 *	13-may-91	V1.3 added helpvec etc.				PJT
 *      23-nov-93       V1.3a cosmetic
 *      21-dec-95       V1.4 process all snapshots
 *	28-may-97	2.0 allow bodytrans and generic function	PJT
 *				(code taken from mkplummer)
 *      24-jul-97       2.1 added norm=
 *	 8-sep-01       a   init_xrandom
 *      15-aug-06       b   prototypes
 */
#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

#include <bodytrans.h>

string defv[] = {
    "in=???\n		      input (snapshot) file",
    "out=???\n                output (snapshot) file",
    "mass=\n		      expression for new masses",
    "inmass=\n                extra file, if to use their masses only",
    "massname=\n              If used,  Mass-function name (mostly n(m))",
    "massexpr=pow(m,p)\n      Mass function expression (e.g. pow(m,p))",
    "masspars=p,0.0\n         Mass function parameters (e.g. p,0.0)",
    "massrange=1,1\n          Range for mass-spectrum (e.g. 1,2)",
    "seed=0\n                 Random seed",
    "norm=\n                  Normalization value for the total mass (if used)",
    "VERSION=2.1b\n           15-aug-06 PJT",
    NULL,
};

string usage="(re)assign masses to a snapshot";

string cvsid="$Id$";

#define TIMEFUZZ	0.0001	/* tolerance in time comparisons */

extern rproc getrfunc(string , string , string , int *);

nemo_main()
{
    stream instr, inmassstr, outstr;
    real   tsnap, tsnapmass, mass_star, mtot, mrange[2], norm;
    Body  *btab = NULL, *bp;
    Body  *bmasstab = NULL, *bmassp;
    int i, n, nbody, nbodymass, bits, bitsmass, seed;
    rproc_body  bfunc;
    rproc       mfunc;
    bool  Qnorm, first = TRUE;

    instr = stropen(getparam("in"), "r");
    if (hasvalue("mass")) {
        dprintf(0,"Using a bodytrans mass expression\n");
        bfunc = btrtrans(getparam("mass"));     /* use bodytrans expression */
        mfunc = NULL;
        inmassstr = NULL;
    } else if (hasvalue("massname")) {
        dprintf(0,"Using a general functional mass expression\n");
        mysymbols(getargv0());
        n=1;
        mfunc = getrfunc(getparam("massname"),
                         getparam("massexpr"),getparam("masspars"),&n);
        if (nemoinpr(getparam("massrange"),mrange,2)!=2)
            error("Need two numbers in massrange=");
        dprintf(1,"Massrange from %f : %f\n",mrange[0],mrange[1]);
        bfunc = NULL;
        inmassstr = NULL;
    } else if (hasvalue("inmass")) {
        inmassstr = stropen(getparam("inmass"),"r");
    } else
    	error("One of: mass=, massname=, inmass= must be given");
    outstr = stropen(getparam("out"), "w");
    seed = init_xrandom(getparam("seed"));
    Qnorm = hasvalue("norm");
    if (Qnorm) norm = getdparam("norm");


    for(;;) {
					/* input main data */
    	get_history(instr);                    /* skip history & comments */
        if (!get_tag_ok(instr, SnapShotTag)) {
            break;
            error("Snapmass (in): Need a snapshot");
        }
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if (bits&MassBit) {
            dprintf(0,"Warning: existing masses overwritten\n");
            if (!Qnorm) {
            	Qnorm = TRUE;
        	for (bp=btab, i=0, norm=0.0; i<nbody; bp++,i++)            
            	    norm += Mass(bp);
            }
        }
        if (!(bits&TimeBit)) {
            dprintf(1,"Warning: time=0.0 assumed\n");
            tsnap = 0.0;
        }
        
        if (inmassstr) { 		/* input mass data from a file */
    	    get_history(inmassstr);
            if (!get_tag_ok(inmassstr, SnapShotTag))
       		error("Snapmass (inmass): Need a snapshot");
            get_snap(inmassstr, &bmasstab, &nbodymass, &tsnapmass, &bitsmass);
            if (nbodymass != nbody) {
	        if (nbodymass < nbody)
                    error("too few bodies (%d < %d)",nbodymass,nbody);
                else
                    warning("too many bodies (%d > %d)",nbodymass,nbody);
            }
            if ((bitsmass & MassBit) == 0) {
                unlink(getparam("out"));
                error("essential mass data missing\tbits = %x\n", bits);
            }
        }
        				/* output */
        if (first) {
            put_history(outstr);
            first = FALSE;
        }
        if (inmassstr) {
	    if (nbody==nbodymass)
	        for (bp=btab, bmassp=bmasstab; bp<btab+nbody; bp++, bmassp++)
	            Mass(bp) = Mass(bmassp);
	    else
	        for (bp=btab, bmassp=bmasstab; bp<btab+nbody; bp++)
	            Mass(bp) = Mass(bmassp);
        } else if (bfunc) {
            for (bp=btab, i=0; i<nbody; bp++,i++)
                Mass(bp) = bfunc(bp, tsnap, i);
        } else if (mfunc) {
            for (bp=btab, i=0; i<nbody; bp++,i++)
	        Mass(bp) = frandom(mrange[0], mrange[1], mfunc);
        } else             
            error("bad flow logic");

        for (bp=btab, i=0, mtot=0.0; i<nbody; bp++,i++)            
            mtot += Mass(bp);
        dprintf(0,"Total mass: %g\n",mtot);
        if (Qnorm) {
            dprintf(0,"Rescaling total mass to %g\n",norm);
            for (bp=btab, i=0; i<nbody; bp++,i++)
                Mass(bp) *= norm/mtot;
        }
	bits |= MassBit;    /* turn mass bit on anyways */
        put_snap(outstr, &btab, &nbody, &tsnap, &bits);
    }

    strclose(instr);
    if (inmassstr) strclose(inmassstr);
    strclose(outstr);
}
