/*
 * QUADFORCE_MAIN.C: quadrupole-order force calculation of an N-body system.
 *
 *	 4-mar-89  V1.0  JEB	Some original version
 *	12-nov-91  V1.1  PJT	New NEMO V2.
 *	 6-may-92      a PJT	extra warning if no output selected
 */

#include "quaddefs.h"

#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

string defv[] = {
    "in=???\n			input file with N-body snapshot",
    "out=\n			output file with forces and potentials",
    "quad=\n			output file with field tables",
    "eps_r=0.05\n		radial softening parameter",
    "eps_t=0.07\n		tangential softening parameter",
    "VERSION=1.1a\n		6-may-92 PJT",
    NULL,
};

string usage="quadrupole-order force calculation of an N-body system.";

nemo_main()
{
    stream instr, outstr, quadstr;
    Body *btab = NULL, *bp;
    int nbody, bits;
    real eps_r, eps_t, tsnap = 0.0;

    if (!hasvalue("out") && !hasvalue("quad"))
        warning("No output created (out=, quad=)");

    instr = stropen(getparam("in"), "r");
    get_history(instr);
    eps_r = getdparam("eps_r");
    eps_t = getdparam("eps_t");
    get_snap(instr, &btab, &nbody, &tsnap, &bits);
    if (bits & MassBit == 0 || bits & PhaseSpaceBit == 0)
	error("not enuf info: bits = 0x%x", bits);
    for (bp = btab; bp < btab+nbody; bp++) {
	CLRV(Acc(bp));
	Phi(bp) = 0.0;
    }
    qfld.nqtab = 0;
    quadforce(btab, nbody, eps_r, eps_t);
    if (hasvalue("out")) {
	outstr = stropen(getparam("out"), "w");
	put_history(outstr);
	bits = bits | PotentialBit | AccelerationBit;
	put_snap(outstr, &btab, &nbody, &tsnap, &bits);
    }
    if (hasvalue("quad")) {
	quadstr = stropen(getparam("quad"), "w");
	put_history(quadstr);
	put_quadfield(quadstr, tsnap, eps_r, eps_t);
    }
}

put_quadfield(quadstr, tsnap, eps1, eps2)
stream quadstr;
real tsnap, eps1, eps2;
{
    put_set(quadstr, "QuadField");
    put_data(quadstr, "Time", RealType, &tsnap, 0);
    put_data(quadstr, "nqtab", IntType, &qfld.nqtab, 0);
    put_data(quadstr, "radtab", RealType, qfld.radtab, qfld.nqtab, 0);
    put_data(quadstr, "Q00tab", RealType, qfld.Q00tab, qfld.nqtab, 0);
    put_data(quadstr, "P00tab", RealType, qfld.P00tab, qfld.nqtab, 0);
    put_data(quadstr, "Q10tab", RealType, qfld.Q10tab, qfld.nqtab, 0);
    put_data(quadstr, "P10tab", RealType, qfld.P10tab, qfld.nqtab, 0);
    put_data(quadstr, "Q11tab", RealType, qfld.Q11tab, NDIM, qfld.nqtab, 0);
    put_data(quadstr, "P11tab", RealType, qfld.P11tab, NDIM, qfld.nqtab, 0);
    put_data(quadstr, "Q22tab", RealType, qfld.Q22tab,
	     NDIM, NDIM, qfld.nqtab, 0);
    put_data(quadstr, "P22tab", RealType, qfld.P22tab,
	     NDIM, NDIM, qfld.nqtab, 0);
    put_data(quadstr, "eps1", RealType, &eps1, 0);
    put_data(quadstr, "eps2", RealType, &eps2, 0);
    put_tes(quadstr, "QuadField");
}
