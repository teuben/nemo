/*
 * SNAPSTAB.C:  report on disk stability / virial etc.
 *              It ignores stuff that's happening in Z (z,vz)
 *              and assumes normal motion is on circular orbits
 *
 *      
 *      29-Mar-89       V1.0    created                 PJT
 *                          b   layout                  PJT
 *	29-oct-90	V1.1    fixed yet another bug	PJT
 *	16-feb-97	V1.2    NEMO V2.x		pjt
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
    "in=???\n           Input file name (snapshot)",
    "times=all\n        Times to select snapshots from",
    "VERSION=1.2\n      16-feb-97 PJT",
    NULL,
};

string usage = "report disk stability/virial";

#define TIMEFUZZ        0.00001 /* tolerance in time comparisons */

nemo_main()
{
    stream instr;
    real   inv_rscale, e_pot, e_kin, tsnap;
    string times;
    Body *btab = NULL, *bp;
    int i, nbody, bits, nrscale, nvscale;
    real   tmporb, tmpran, e_orb, e_ran;

    times = getparam("times");
    instr = stropen(getparam("in"), "r");

    get_history(instr);
    for (;;) {
        if (!get_tag_ok(instr, SnapShotTag))
                break;                  /* done with work in loop */
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if ((bits & MassBit) == 0 && (bits & PhaseSpaceBit) == 0) {
            continue;       /* just skip */
        }

        if ((bits & TimeBit) == 0)
                tsnap = 0.0;
        else if (!streq(times,"all") && !within(tsnap, times, TIMEFUZZ))
                continue;               /* now skip this snapshot */
        dprintf (1,"Snapshot time=%f\n",tsnap);
        if ((bits & PotentialBit)) {
            e_pot = e_kin = e_orb = 0.0;
            for (bp = btab; bp < btab+nbody; bp++) {
                e_pot += Mass(bp) * Phi(bp);
                split_e (Mass(bp), Pos(bp), Vel(bp), &tmporb, &tmpran);
                e_kin += tmporb + tmpran;
                e_orb += tmporb;
            }
        } else
            error("snapstab: Potentials required\n");
        e_kin *= 0.5;   /* proper factor 0.5  */
        e_pot *= 0.5;   /* double count */
        e_orb *= 0.5;
        if (e_pot > 0.0)
            error("positive total potential energy???");
        if (e_kin < 0.0)
            error("negative total kinetic energy???");

        printf("t= %g E=T+U:  %g %g %g   2T+U= %g t_OP= %g\n",
            tsnap,e_pot+e_kin,e_kin,e_pot,2*e_kin+e_pot,e_orb/ABS(e_pot));
    }
}

split_e (m, r, v, eorb, eran)           /* really a kludge */
real m, *eorb, *eran;
vector r, v;
{
    real rad, rv, vtot, vran;
    int    i;

    rad = sqr(r[0]) + sqr(r[1]);      /* only in the XY plane */
    if (rad==0.0) {                   /* put all in random motion */
        *eorb = 0.0;
	*eran = 0.0;
        for (i=0; i<NDIM; i++)
            *eran += m * sqr(v[i]);
    } else {
        rv = r[0]*v[1] - r[1]*v[0];
        vran = rv*rv / rad;
        vtot = sqr(v[0]) + sqr(v[1]) + sqr(v[2]);
        *eran = m * vran;
        *eorb = m * (vtot-vran);
    }
}
