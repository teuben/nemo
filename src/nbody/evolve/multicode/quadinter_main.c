/*
 * QUADINTER_MAIN.C: quadrupole-order force calculation from tabulated field.
 *
 *	 4-mar-89  V1.0 JEB	Some original version
 *	12-nov-91  V1.1 PJT	New NEMO V2.
 *       6-may-92  V1.1a PJT    usage, hasvalue, extra warning if no output
 */

#include "quaddefs.h"

#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

string defv[] = {
    "quad=???\n			input file with field tables",
    "in=???\n			input file with N-body snapshot",
    "out=\n			output file with forces and potentials",
    "VERSION=1.1\n		12-nov-91 PJT",
    NULL,
};

string usage="quadrupole-order force calculation from tabulated field.";

nemo_main()
{
    stream instr, quadstr, outstr;
    Body *btab = NULL, *bp;
    int nbody, bits;
    real eps_r, eps_t, tsnap = 0.0;

    if (!hasvalue("out"))
        warning("No output supplied (out=)");

    quadstr = stropen(getparam("quad"), "r");
    get_history(quadstr);
    instr = stropen(getparam("in"), "r");
    get_history(instr);
    get_quadfield(quadstr, &eps_r, &eps_t);
    get_snap(instr, &btab, &nbody, &tsnap, &bits);
    if (bits & PhaseSpaceBit == 0)
	error("not enuf info: bits = %o", bits);
    for (bp = btab; bp < btab+nbody; bp++) {
	CLRV(Acc(bp));
	Phi(bp) = 0.0;
    }
    quadinter(btab, nbody, eps_r, eps_t);
    if (hasvalue("out")) {
	outstr = stropen(getparam("out"), "w");
	put_history(outstr);
	bits = bits | PotentialBit | AccelerationBit;
	put_snap(outstr, &btab, &nbody, &tsnap, &bits);
    }
}

get_quadfield(quadstr, eps1, eps2)
stream quadstr;
real *eps1, *eps2;
{
    get_set(quadstr, "QuadField");
    get_data(quadstr, "nqtab", IntType, &qfld.nqtab, 0);
    get_data(quadstr, "radtab", RealType, qfld.radtab, qfld.nqtab, 0);
    get_data(quadstr, "Q00tab", RealType, qfld.Q00tab, qfld.nqtab, 0);
    get_data(quadstr, "P00tab", RealType, qfld.P00tab, qfld.nqtab, 0);
    get_data(quadstr, "Q10tab", RealType, qfld.Q10tab, qfld.nqtab, 0);
    get_data(quadstr, "P10tab", RealType, qfld.P10tab, qfld.nqtab, 0);
    get_data(quadstr, "Q11tab", RealType, qfld.Q11tab, NDIM, qfld.nqtab, 0);
    get_data(quadstr, "P11tab", RealType, qfld.P11tab, NDIM, qfld.nqtab, 0);
    get_data(quadstr, "Q22tab", RealType, qfld.Q22tab,
	     NDIM, NDIM, qfld.nqtab, 0);
    get_data(quadstr, "P22tab", RealType, qfld.P22tab,
	     NDIM, NDIM, qfld.nqtab, 0);
    get_data(quadstr, "eps1", RealType, eps1, 0);
    get_data(quadstr, "eps2", RealType, eps2, 0);
    get_tes(quadstr, "QuadField");
}
