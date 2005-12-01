/*
 * SNAPSCALE.C: rescale an N-body snapshot.
 *
 *	<< also used in benchmarking speed of snap i/o and math operations >>
 *	
 *	16-May-88	V1.0	created 		JEB
 *	 1-Jun-88	V2.0	added old PJT's options	PJT
 *	 6-jun-88	V2.0a	using libT for new getparam/history	PJT
 *	23-aug-88	V2.1	added pscale, accel's will be wrong then PJT
 *      16-Nov-88       V2.2    correct for weird datafiles PJT
 *      31-mar-89       V2.3    nemoind parsing PJT
 *      15-aug-89       V2.4    added ascale, pscale    PJT
 *	15-nov-90	V3.0    NEMO 2.x PJT
 *	27-apr-92	V3.1    toyed with benchmarking; added the Qxxx bool`s	PJT
 *				A 20,000 body snapshot on an Sparc IPX:
 *					  old    new
 *				gcc -g:	  1.82   1.48 sec.
 *				gcc -O2:  0.78   0.69 sec.
 *	16-feb-97	V3.1a   fixed for SINGLEPREC
 *       8-oct-01        3.2    add Dens and Eps
 *       8-aug-05           a   fix aux normalization 
 *       1-dec-05           b   fix softening and density scaling
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

string defv[] = {
    "in=???\n	    Input file name",
    "out=???\n	    Output file name",
    "mscale=1.0\n   Mass scale factor [1.0]",
    "rscale=1.0\n   Position scale factor [1.0[,1.0,1.0]]",
    "vscale=1.0\n   Velocity scale factor [1.0[,1.0,1.0]]",
    "pscale=1.0\n   Potential scale factor [1.0]",
    "ascale=1.0\n   Accelerations scale factor [1.0[,1.0,1.0]]",
    "xscale=1.0\n   Aux scale factor [1.0]",
    "kscale=1\n     Key scale factor [1]",
    "dscale=1\n     Dens scale factor",
    "escale=1\n     Eps scale factor",
    "times=all\n    Times to select snapshots from",
    "VERSION=3.2b\n 1-dec-05 PJT",
    NULL,
};

string usage="rescale a snapshot";

string cvsid="$Id$";

#define TIMEFUZZ	0.000001	/* tolerance in time comparisons */

bool uscalar(real x)
{
    return x==1.0;
}

bool uvector(vector v)
{
    register int i;

    for (i=0; i<NDIM; i++)
        if(v[i] != 1.0) return FALSE;
    return TRUE;
}


void nemo_main()
{
    stream instr, outstr;
    real   mscale, pscale, xscale, tsnap, escale, dscale;
    vector rscale, vscale, ascale;
    string times;
    int i, nbody, bits, nrscale, nvscale, nascale, kscale;
    Body *btab = NULL, *bp;
    bool Qmass, Qphase, Qacc, Qpot, Qkey, Qaux, Qeps, Qdens;

    nrscale = nemoinpr(getparam("rscale"),rscale,NDIM);     /* RSCALE */
    if (nrscale==1) 
    	for (i=1; i<NDIM; i++)
    	   rscale[i] = rscale[0];
    else if (nrscale!=NDIM)
    	error("keyword rscale needs either 1 or %d numbers", NDIM);
    			

    nvscale = nemoinpr(getparam("vscale"),vscale,NDIM);     /* VSCALE */
    if (nvscale==1)
    	for (i=1; i<NDIM; i++)
    	   vscale[i] = vscale[0];
    else if (nvscale!=NDIM)
    	error("keyword vscale needs either 1 or %d numbers", NDIM);    

    nascale = nemoinpr(getparam("ascale"),ascale,NDIM);     /* ASCALE */
    if (nascale==1)
    	for (i=1; i<NDIM; i++)
    	   ascale[i] = ascale[0];
    else if (nascale!=NDIM)
    	error("keyword ascale needs either 1 or %d numbers", NDIM);    

    mscale = getdparam("mscale");
    pscale = getdparam("pscale");
    xscale = getdparam("xscale");
    dscale = getdparam("dscale");
    escale = getdparam("escale");
    kscale = getiparam("kscale");
    times = getparam("times");

    instr = stropen(getparam("in"), "r");   /* open files */
    outstr = stropen(getparam("out"), "w");

    get_history(instr);
    put_history(outstr);		
    for (;;) {
    	get_history(instr);		/* skip over stuff we can forget */
        if (!get_tag_ok(instr, SnapShotTag))
		break;			/* done with work in loop */
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if ((bits & MassBit) == 0 && (bits & PhaseSpaceBit) == 0) {
	    continue;       /* just skip it's probably a diagnostics */
        }

        if ((bits & TimeBit) == 0)
	    tsnap = 0.0;
        else if (!streq(times,"all") && !within(tsnap, times, TIMEFUZZ))
            continue;
        dprintf (1,"Scaling snapshot at time= %f bits=0x%x\n",tsnap,bits);

        Qmass  = MassBit & bits        && !uscalar(mscale);
        Qphase = PhaseSpaceBit & bits  &&(!uvector(rscale) || !uvector(vscale));
        Qacc   = AccelerationBit & bits&& !uvector(ascale);
        Qpot   = PotentialBit & bits   && !uscalar(pscale);
        Qaux   = AuxBit & bits         && !uscalar(xscale);
        Qkey   = KeyBit & bits         && (kscale!=1);
        Qdens  = DensBit & bits        && !uscalar(dscale);
        Qeps   = EpsBit & bits         && !uscalar(escale);

        dprintf(1,"Scaling: ");
        if (Qmass)  dprintf(1," mass");
        if (Qphase) dprintf(1," phase");
        if (Qacc)   dprintf(1," acc");
        if (Qpot)   dprintf(1," pot");
        if (Qaux)   dprintf(1," aux");
        if (Qkey)   dprintf(1," key");
        if (Qdens)  dprintf(1," dens");
        if (Qeps)   dprintf(1," eps");
        dprintf(1,"\n");

        if (Qmass || Qphase || Qacc || Qpot || Qaux || Qkey || Qdens || Qeps) {
            for (bp = btab; bp < btab+nbody; bp++) {
                if(Qmass) Mass(bp) *= mscale;
                if(Qphase) {
                    SMULVV(Pos(bp),rscale);
                    SMULVV(Vel(bp),vscale);
	        }
                if(Qpot) Phi(bp) *= pscale;
                if(Qacc) {
                    SMULVV(Acc(bp),ascale);
                }
                if(Qaux) Aux(bp) *= xscale;
                if(Qkey) Key(bp) *= kscale;
		if(Qdens) Dens(bp) *= dscale;
		if(Qeps) Eps(bp) *= escale;
            }
        } else {
            warning("No scaling applied to snapshot");
	}

        put_snap(outstr, &btab, &nbody, &tsnap, &bits);
    }
}

