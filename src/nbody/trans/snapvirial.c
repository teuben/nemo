/*
 * SNAPVIRIAL.C: rescale an N-body snapshot for virialization
 *		 requires potentials and forces to be present
 *	
 *	11-Mar-89	V1.0	created 		PJT
 *	 6-apr-89	V1.1    2t/w added (did this have a bug?)   PJT
 *	11-aug-89	V1.2	massscale implemented		PJT
 *	29-oct-90	V1.3    fixed decision bug			PJT
 *	 8-nov-91	V1.3a   Clausius virial, if possible		PJT
 *      22-nov-93       V1.4    proper Clausius option                  PJT
 *				(was never implemented yet)
 *      14-aug-96       V1.4    source code moved from nbody/{reduc->trans} PJT
 *       5-jun-97       V1.4a   fixed typo, removed nested externs
 *      18-jul-2012     V2.0    allow out= to be optional, so it only reports    PJT
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
    "in=???\n		Input snapshot",
    "out=\n		(optional) Output snapshot",
    "mscale=f\n		Scale masses ?",
    "rscale=t\n		Scale positions ?",
    "vscale=t\n		Scale velocities ?",
    "times=all\n	Times of snapshots to select",
    "virial=\n		New virial ratio (|2T/W|), if to be changed",
    "VERSION=1.4a\n	6-jun-97 PJT",
    NULL,
};

string usage="rescale snapshot for re-virialization";

#define TIMEFUZZ	0.0001	/* tolerance in time comparisons */

nemo_main()
{
    stream instr, outstr;
    real   mscale, rscale, vscale, e_pot, e_kin, tsnap, virial;
    real   phi_scale, acc_scale;
    string times, vstr;
    Body   *btab = NULL, *bp;
    int    i, nbody, bits, nrscale, nvscale;
    bool   Qvirial, Qmsc, Qrsc, Qvsc;         /* boolean checks for scalings */
    bool   Qout;

    times = getparam("times");
    instr = stropen(getparam("in"), "r");
    Qout = hasvalue("out");

    Qmsc = getbparam("mscale");
    Qrsc = getbparam("rscale");
    Qvsc = getbparam("vscale");
    Qvirial = hasvalue("virial");

    if (Qout)  {
      outstr = stropen(getparam("out"), "w");
      if (Qmsc && Qrsc && Qvsc)
        error("Cannot scale m, r and v all at same time");
      else if (!Qmsc && !Qrsc && !Qvsc)
        error("Nothing being scaled is not very productive");
      if (Qvirial) {
        virial = getdparam("virial");	/* get the actual value */
        if (virial<0.0) {               /* catch illegal ratios */
            warning("virial=%g ?; old virial ratio retained",virial);
	    Qvirial = FALSE;
        }
      }
    } else {
      Qmsc = Qrsc = Qvsc = FALSE;
      if (Qvirial) error("Cannot specify virial in query mode");
    }

    get_history(instr);			        /* get history */
    if (Qout) put_history(outstr);

    for (;;) {				/* loop through snapshots */
        if (!get_tag_ok(instr, SnapShotTag))
		break;			         /* until done */
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if ((bits & MassBit) == 0 && (bits & PhaseSpaceBit) == 0) {
            dprintf (2,"Time= %f auto skipping ",tsnap);
	    continue;       /* just skip - it maybe diagnostics */
        }
        if ((bits & TimeBit) == 0)
        	tsnap = 0.0;
        else if (!streq(times,"all") && !within(tsnap, times, TIMEFUZZ))
        	continue;		/* skip this snapshot too */
        dprintf (2,"Time= %f ",tsnap);
        if ((bits & PotentialBit)) {
            e_pot = e_kin = 0.0;
            for (bp = btab; bp < btab+nbody; bp++) {
                e_pot += Mass(bp) * Phi(bp);
                for (i=0; i<NDIM; i++)
                    e_kin += Mass(bp)*sqr(Vel(bp)[i]);
            }
        } else
            error("Potentials required");
        e_kin *= 0.5;   /* proper factor 0.5 for kinetic energy */
        e_pot *= 0.5;   /* and correct this for double count */
        if (e_pot >= 0.0)
            error("Positive total potential energy???");
        if (e_kin <= 0.0)
            error("Negative total kinetic energy???");
        mscale = rscale = vscale = 1.0;
        if (!Qmsc) {                /* no mass scale */
            if (Qvsc && !Qrsc) {
	        if(Qvirial)
                    vscale = sqrt(-0.5*e_pot/e_kin * virial);
                else
                    error("virial ratio not set (1: m=f r=f v=t)");
            } else if (Qvsc && Qrsc) {
                rscale = -2.0 * e_pot;
                vscale = 0.5 / sqrt(e_kin);
            } else if (Qrsc && !Qvsc) {
                if (Qvirial)
                    rscale = -0.5*e_pot/e_kin * virial; 
                else
                    error("virial ratio not set (2: m=f r=t v=f)");
	    }
        } else {
	    error("Mscale not implemented");
        }
        phi_scale = mscale/rscale;
        acc_scale = phi_scale/rscale;
	dprintf(0,"U= %g T= %g rscale= %g vscale= %g virial=%g\n",
		e_pot,e_kin,rscale,vscale,-2*e_kin/e_pot);
        for (bp = btab; bp < btab+nbody; bp++) {
            Mass(bp)   *= mscale;
	    Phi(bp)    *= phi_scale;
            for (i=0; i<NDIM; i++) {
               Vel(bp)[i] *= vscale;
               Pos(bp)[i] *= rscale;
               Acc(bp)[i] *= acc_scale;
            }
        }
        if (bits&AuxBit)
            warning("Aux information unscaled");
        if (bits&KeyBit)
            warning("Key information unscaled");
	if (Qout) put_snap(outstr, &btab, &nbody, &tsnap, &bits);
#if 1
        free(btab);
        btab = NULL;
#endif
    }
}
