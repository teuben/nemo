/*
 * PLUMMER.C: tabulate Osipkov-Merrirr generalization of Plummer model.
 *		See also mkommod
 *
 *	22-jan-89  V1.1 				Josh Barnes
 *	15-nov-90  V1.2 helpvec and into NEMO 2.x	PJT
 *	23-mar-97  V1.3 usage, cleanup proto's		pjt
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>

string defv[] = {		/* DEFAULT INPUT PARAMETERS                 */
    "out=???\n			output file, written in binary format",
    "anisorad=-1.0\n		anisotropy radius (if le 0, isotropic)",
    "mcutoff=0.999999\n		fractional mass cutoff for finite radius",
    "VERSION=1.2\n		15-nov-90 PJT",
    NULL,
};

string usage = "tabulate Osipkov-Merrirr generalization of Plummer model";

/*
 * The structure of the model to generate is specified by the following.
 */

real anisorad;			/* anisotropy radius (if le 0, isotropic)   */

real mcutoff;			/* fractional mass cutoff                   */

#define MXTB 512		/* max number of values per table           */

int ntab;			/* count of tabulated values                */

real rtab[MXTB];		/* table of radii                           */

real mtab[MXTB];		/* table of interior masses                 */

real ptab[MXTB];		/* table of potentials                      */

real ftab[MXTB];		/* table of distribution function values    */

nemo_main()
{
    anisorad = getdparam("anisorad");
    mcutoff = getdparam("mcutoff");
    ntab = MXTB;
    plummer();
    writemodel(getparam("out"));
}

/*
 * PLUMMER: tabulate Osipkov-Merrirr generalization of Plummer model.
 */

plummer()
{
    real m0, sig0, phi0, r0, phim, s;
    int i;

    m0 = 1.0;						/* unit mass system */
    sig0 = sqrt(8.0 / (9.0 * PI * m0));			/* binding E = -1/4 */
    phi0 = -6.0 * sqr(sig0);				/* pot. at center   */
    r0 = - m0 / phi0;					/* scaling radius   */
    phim = phi0 * sqrt(1.0 - pow(mcutoff, 2.0/3.0));	/* pot. at rmax     */
    for (i = 0; i < ntab; i++) {			/* loop over tables */
	ptab[i] = phi0 - (phi0 - phim) * i / (ntab - 1.0);
	rtab[i] = r0 * sqrt(sqr(phi0/ptab[i]) - 1.0);
	mtab[i] = m0 * qbe(rtab[i]/r0) * qbe(ptab[i]/phi0);
	if (anisorad > 0.0) {				/*   aniso. model?  */
	    s = 1.0 - (1.0 - 15.75 * sqr(sqr(sig0)/ptab[i])) *
		        sqr(r0/anisorad);
	    if (s < 0.0)
		error("plummer: f < 0  anisorad = %g\tQ = %g",
		      anisorad, ptab[i]);
	} else
	    s = 1.0;
	ftab[i] = s * pow(- ptab[i] / sqr(sig0), 3.5);
    }
}

/*
 * WRITEMODEL: output model tables in binary format.
 */

writemodel(name)
string name;
{
    stream outstr;

    outstr = stropen(name, "w");
    put_history(outstr);
    put_set(outstr, "OsipkovMerrittModel");
    put_data(outstr, "AnisoRadius", RealType, &anisorad, 0);
    put_data(outstr, "Ntab", IntType, &ntab, 0);
    put_data(outstr, "Radius", RealType, rtab, ntab, 0);
    put_data(outstr, "Mass", RealType, mtab, ntab, 0);
    put_data(outstr, "Potential", RealType, ptab, ntab, 0);
    put_data(outstr, "DistFunc", RealType, ftab, ntab, 0);
    put_tes(outstr, "OsipkovMerrittModel");
    strclose(outstr);
}
