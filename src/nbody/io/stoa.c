/*
 * STOA.C: convert binary N-body file to ASCII format.
 *         "I did it my way" --  F. Sinatra.
 *
 *  11-may-88    -- JEB
 *  23-oct-91    %21.13 is better default for rformat?
 *  20-feb-92	 2.2 usage, nemo_main; helpvec, local allocate() removed  PJT
 *   7-mar-92    happy gcc2.0
 *   9-nov-92    make it skip SnapShot's which have no Particles 
 *               (e.g. Diagnostics)  - fprintf(stderr, -> dprintf(0,      PJT
 *  22-feb-94    ansi headers
 */

#include <stdinc.h>
#include <assert.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <snapshot/snapshot.h>

string defv[] = {               /* DEFAULT INPUT PARAMETERS                 */
    "in=???\n                     input file name",
    "out=???\n                    output file name",
    "options=mass,phase\n         information to output {mass,phase,phi}",
    "iformat=  %d\n		  output format for integers",
    "rformat= %21.13E\n		  output format for floating-point",
    "VERSION=2.3c\n               22-feb-94 PJT",
    NULL,
};

string usage = "convert binary N-body file to ASCII format";

string iformat, rformat;

void out_int(), out_real(), out_vec();

void
nemo_main()
{
    stream instr, outstr;
    string options;
    int nbody, i, ndiag=0;
    bool scanopt();
    real tsnap, *mbuf, *mp, *cbuf, *cp, *pbuf, *pp;

    instr = stropen(getparam("in"), "r");
    get_history(instr);
    outstr = stropen(getparam("out"), "w");
    options = getparam("options");
    iformat = getparam("iformat");
    rformat = getparam("rformat");
    while (get_tag_ok(instr, SnapShotTag)) {	/* loop reading data frames */
	get_set(instr, SnapShotTag);
	get_set(instr, ParametersTag);
	if (get_tag_ok(instr, TimeTag))
	    get_data_coerced(instr, TimeTag, RealType, &tsnap, 0);
	else {
	    tsnap = 0.0;
	    dprintf(0, "[%s: assuming %s = %f]\n",
		    getargv0(), TimeTag, tsnap);
	}
	get_data_coerced(instr, NobjTag, IntType, &nbody, 0);
	get_tes(instr, ParametersTag);
        if (!get_tag_ok(instr, ParticlesTag)) {
            ndiag++;
	    get_tes(instr, SnapShotTag);
            continue;
        }
        dprintf(0, "[%s: writing %d bodies at time %f]\n",      /*  report  */
		getargv0(), nbody, tsnap);
	out_int(outstr, nbody);			/*   write frame header     */
	out_int(outstr, NDIM);
	out_real(outstr, tsnap);
	get_set(instr, ParticlesTag);
        if (scanopt(options, "mass")) {		/*   mass data to convert?  */
	    mbuf = mp = (real *) allocate(nbody * sizeof(real));
	    get_data_coerced(instr, MassTag, RealType, mbuf, nbody, 0);
	    for (i = 0; i < nbody; i++) {
                out_real(outstr, *mp);
		mp++;
	    }
	    free(mbuf);
        }
        if (scanopt(options, "phase")) {	/*   phase-space data?      */
	    cbuf = (real *) allocate(nbody * 2 * NDIM * sizeof(real));
	    get_data_coerced(instr, PhaseSpaceTag, RealType, cbuf,
			     nbody, 2, NDIM, 0);
	    cp = cbuf;				/*     write out positions  */
	    for (i = 0; i < nbody; i++) {
		out_vec(outstr, cp);
		cp += 2 * NDIM;
	    }
	    cp = cbuf + NDIM;			/*     write out velocities */
	    for (i = 0; i < nbody; i++) {
		out_vec(outstr, cp);
		cp += 2 * NDIM;
	    }
	    free(cbuf);
	}
        if (scanopt(options, "phi")) {		/*   potential data?        */
	    pbuf = pp = (real *) allocate(nbody * sizeof(real));
	    get_data_coerced(instr, PotentialTag, RealType, pbuf, nbody, 0);
	    for (i = 0; i < nbody; i++) {
                out_real(outstr, *pp);
		pp++;
	    }
	    free(pbuf);
        }
	get_tes(instr, ParticlesTag);
	get_tes(instr, SnapShotTag);
    } /* while more snapshots in input file */
    strclose(outstr);
    if (ndiag>0) 
       dprintf(0,"Skipped %d snapshots that did not contain particles\n",ndiag);
}

#define MBUF  64			/* size of temp format buffer       */

void
out_int(str, ival)
stream str;
int ival;
{
    char buf[MBUF];

    sprintf(buf, "%s\n", iformat);
    fprintf(str, buf, ival);
}

void
out_real(str, rval)
stream str;
real rval;
{
    char buf[MBUF];

    sprintf(buf, "%s\n", rformat);
    fprintf(str, buf, rval);
}

void
out_vec(str, vec)
stream str;
vector vec;
{
    char buf[MBUF];

    sprintf(buf, "%s%s%s\n", rformat, rformat, rformat);
    fprintf(str, buf, vec[0], vec[1], vec[2]);
}
