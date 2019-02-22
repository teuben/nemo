/*
 * PLUMMER.C: tabulate Osipkov-Merritt generalization of Plummer model.
 *	      See also mkommod
 *
 *	22-jan-89   V1.1 				Josh Barnes
 *	15-nov-90   V1.2 helpvec and into NEMO 2.x	PJT
 *	23-mar-97   V1.3 usage, cleanup proto's		pjt
 *      28-aug-2018 V2.0 ntab=, some code cleanup       pjt
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#include <history.h>

string defv[] = {		/* DEFAULT INPUT PARAMETERS                 */
    "out=???\n			output file, written in snapshot format",
    "anisorad=-1.0\n		anisotropy radius (<=0 means isotropic)",
    "mcutoff=0.999999\n		fractional mass cutoff for finite radius",
    "ntab=512\n                 number of table entries (<=MXTB)",
    "VERSION=2.0\n		28-aug-2018 PJT",
    NULL,
};

string usage = "tabulate Osipkov-Merritt generalization of Plummer model";

/*
 * The structure of the model to generate is specified by the following:
 */

#define MXTB 10000		/* max number of values per table           */

real anisorad;			/* anisotropy radius (if le 0, isotropic)   */
real mcutoff;			/* fractional mass cutoff                   */
int ntab;			/* count of tabulated values                */
real rtab[MXTB];		/* table of radii                           */
real mtab[MXTB];		/* table of interior masses                 */
real ptab[MXTB];		/* table of potentials                      */
real ftab[MXTB];		/* table of distribution function values    */

local void writemodel(string name);
local void plummer(void);

void nemo_main()
{
    anisorad = getdparam("anisorad");
    mcutoff = getdparam("mcutoff");
    ntab = getiparam("ntab");
    if (ntab > MXTB) {
      ntab = MXTB;
      warning("Resetting ntab=%d to MXTB",ntab);
    }
    plummer();
    writemodel(getparam("out"));
}

/*
 * PLUMMER: tabulate Osipkov-Merritt generalization of Plummer model.
 */

void plummer(void)
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
		error("plummer: f < 0  anisorad=%g\tQ=%g",anisorad, ptab[i]);
	} else
	    s = 1.0;
	ftab[i] = s * pow(- ptab[i] / sqr(sig0), 3.5);
    }
}

/*
 * WRITEMODEL: output model tables in binary snapshot format.
 */

void writemodel(string name)
{
    stream outstr;

    outstr = stropen(name, "w");
    put_history(outstr);
    put_set(outstr, "OsipkovMerrittModel");
    put_data(outstr, "AnisoRadius", RealType, &anisorad,  0);
    put_data(outstr, "Ntab",        IntType,  &ntab,      0);
    put_data(outstr, "Radius",      RealType, rtab, ntab, 0);
    put_data(outstr, "Mass",        RealType, mtab, ntab, 0);
    put_data(outstr, "Potential",   RealType, ptab, ntab, 0);
    put_data(outstr, "DistFunc",    RealType, ftab, ntab, 0);
    put_tes(outstr, "OsipkovMerrittModel");
    strclose(outstr);
}
