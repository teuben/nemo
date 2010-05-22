/*
 * ATOS.C: convert ASCII N-body file to binary format.
 *         "I did it my way" --  F. Sinatra.
 *
 *	28-mar-89   V2.4    --apparent previous version by JEB? --
 *	14-mar-90   V2.5	made GCC happy, added helpvec	[PJT]
 *	24-oct-91   V2.5a   make it accept interspersed comma's (',') too  PJT
 *                          oeps, didn't work, just a few small changes
 *	 7-mar-92   V2.5b   happy gcc2.0				   pjt
 *	22-feb-94       c   ansi headers
 *	 2-sep-94   V2.6    allow new format
 *      22-may-10   V2.7    remove locally (re)defined Ngas etc., use snapshot.h
 */

#include <stdinc.h>
#include <assert.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>

string defv[] = {               /* DEFAULT INPUT PARAMETERS                 */
    "in=???\n                               Input file name",
    "out=???\n                              Output file name (in append mode)",
#if !defined(SPH)
    "options=mass,phase\n                   Information to output",
#else
    "options=mass,phase,rho,temp,hsph,phi\n Information to output",
#endif
    "hexpack=false\n    If true, read packed hex data",
    "headline=\n        Random mumblage for humans",
#if defined(SPH)
    "mode=1\n		1: header has 4 items      2: header has 5 items",
#endif
    "VERSION=2.7\n      22-may-2010 PJT",
    NULL,
};

string usage = "convert ASCII N-body file to binary format";

stream instr, outstr;

string options;
int mode;

bool hexpack;

int coordsys = CSCode(Cartesian, NDIM, 2);
int linecnt=0;

bool in_header();
void conv_real(), conv_vect(), conv_phase(), 
     in_real(), in_vect(), hexcheck();



void
nemo_main()
{
    int nbody, ngas;
    real tsnap;

    instr = stropen(getparam("in"), "r");
    outstr = stropen(getparam("out"), "w");     /*  "w" or "a" */
#if defined(SPH)    
    mode = getiparam("mode");
#endif    
    if (! streq(getparam("headline"), "")) {
	set_headline(getparam("headline"));
    }
    put_history(outstr);
    options = getparam("options");
    hexpack = getbparam("hexpack");
    while (in_header(&nbody, &ngas, &tsnap)) {	/* loop reading data frames */
        dprintf(0, "[%s: reading %d bodies at time %f]\n",
		getargv0(), nbody, tsnap);
	put_set(outstr, SnapShotTag);		/*   start snapshot output  */
	put_set(outstr, ParametersTag);
	put_data(outstr, TimeTag, RealType, &tsnap, 0);
	put_data(outstr, NobjTag, IntType, &nbody, 0);
#if defined(SPH)
	put_data(outstr, NgasTag, IntType, &ngas, 0);
#endif
	put_tes(outstr, ParametersTag);
	put_set(outstr, ParticlesTag);		/*   start particle output  */
	put_data(outstr, CoordSystemTag, IntType, &coordsys, 0);
        if (scanopt(options, "mass"))		/*   convert mass data      */
	    conv_real(MassTag, nbody);
	if (scanopt(options, "pos")) {		/*   convert position data  */
	    assert(! scanopt(options, "phase"));
	    conv_vect(PosTag, nbody);
	}
	if (scanopt(options, "vel")) {		/*   convert velocity data  */
	    assert(! scanopt(options, "phase"));
	    conv_vect(VelTag, nbody);
	}
        if (scanopt(options, "phase"))		/*   convert phase-space    */
	    conv_phase(PhaseSpaceTag, nbody);
#if defined(SPH)
	if (scanopt(options, "rho"))		/*   convert density data   */
	    conv_real(DensityTag, ngas);
	if (scanopt(options, "temp"))		/*   convert temperature    */
	    conv_real(TemperatureTag, ngas);
	if (scanopt(options, "hsph"))		/*   convert smoothing      */
	    conv_real(SmoothLengthTag, ngas);
#endif
        if (scanopt(options, "phi"))		/*   convert potential data */
	    conv_real(PotentialTag, nbody);
        if (scanopt(options, "acc"))		/*   convert acceleration   */
	    conv_vect(AccelerationTag, nbody);
	put_tes(outstr, ParticlesTag);
	put_tes(outstr, SnapShotTag);
    }
    strclose(outstr);
} /* nemo_main */

#if !defined(SPH)

bool 
in_header(nbptr, ngptr, tsptr)
int *nbptr, *ngptr;
real *tsptr;
{
    int ndim;
    double dval;

    linecnt++;
    if (fscanf(instr, " %d %d %lf\n", nbptr, &ndim, &dval) != 3)
	return (FALSE);
    if (ndim != NDIM)
	error("%s: got ndim = %d, not %d\n", getargv0(), ndim, NDIM);
    *ngptr = 0;
    *tsptr = dval;
    return (TRUE);
}

#else

bool in_header(nbptr, ngptr, tsptr)
int *nbptr, *ngptr;
real *tsptr;
{
    int ndim, ndummy;
    double dval;

    linecnt++;
    if (mode==1) {
        /* old 'standard' format: NBODY, NGAS, NDIM, TIME */
        if (fscanf(instr, " %d %d %d %lf\n", 
            nbptr, ngptr, &ndim, &dval) != 4)
	    return (FALSE);
    } else {
        /* some new format: NBODY, NGAS, NBODY-NGAS, NDIM, TIME */
        if (fscanf(instr, " %d %d %d %d %lf\n", 
    	    nbptr, ngptr, &ndummy, &ndim, &dval) != 5)
	    return (FALSE);
    }
    if (ndim != NDIM)
	error("%s: got ndim = %d, not %d\n", getargv0(), ndim, NDIM);
    *tsptr = dval;
    return (TRUE);
}

#endif

void
conv_real(tag, ndata)
string tag;
int ndata;
{
    real *rbuf, *rp;
    int i;

    rp = rbuf = (real *) allocate(ndata * sizeof(real));
    hexcheck();					/* check for hex data input */
    for (i = 0; i < ndata; i++) {
	in_real(rp);
	rp++;
    }
    put_data(outstr, tag, RealType, rbuf, ndata, 0);
    free(rbuf);
}

void
conv_vect(tag, ndata)
string tag;
int ndata;
{
    real *vbuf, *vp;
    int i;

    vp = vbuf = (real *) allocate(ndata * NDIM * sizeof(real));
    hexcheck();
    for (i = 0; i < ndata; i++) {
	in_vect(vp);
	vp += NDIM;
    }
    put_data(outstr, tag, RealType, vbuf, ndata, NDIM, 0);
    free(vbuf);
}

void
conv_phase(tag, ndata)
string tag;
int ndata;
{
    real *pbuf, *pp;
    int i;

    pp = pbuf = (real *) allocate(ndata * 2 * NDIM * sizeof(real));
    hexcheck();
    for (i = 0; i < ndata; i++) {	/*     read in positions    */
	in_vect(pp);
	pp += 2 * NDIM;
    }
    pp = pbuf + NDIM;			/*     read in velocities   */
    hexcheck();
    for (i = 0; i < ndata; i++) {
	in_vect(pp);
	pp += 2 * NDIM;
    }
    put_data(outstr, tag, RealType, pbuf, ndata, 2, NDIM, 0);
    free(pbuf);
}

local int hexmax;		/* max value of packed hex coordinates      */
local double scale;		/* scale factor for packed hex coords.      */

void 
hexcheck()
{
    if (hexpack) {		        	/* packed hex data to read? */
        linecnt++;
	if (fscanf(instr, " %lf\n", &scale) != 1)
	    error("hexcheck in %s: cannot read scale\n", getargv0());
	linecnt++;
	if (fscanf(instr, " %d\n", &hexmax) != 1)
	    error("hexcheck in %s: cannot read hexmax\n", getargv0());
    }
}

void
in_real(rptr)
real *rptr;
{
    char line[128];
    double dval;
    int ival;

    linecnt++;
    if (fgets(line, 128, instr) == NULL)
	error("in_real in %s: unexpected EOF\n", getargv0());
    if (! hexpack) {
	if (sscanf(line, " %lf\n", &dval) != 1)
	    error("in_real in %s: cannot read line %d: %s", 
	                getargv0(), linecnt, line);
	*rptr = dval;
    } else {
	if (sscanf(line, " %x\n", &ival) != 1)
	    error("in_real in %s: cannot read line %d: %s", 
	            getargv0(), linecnt, line);
	*rptr = (2 * scale / hexmax) * ival - scale;
    }	
}

void
in_vect(vptr)
vector vptr;
{
    char line[128];
    double dx, dy, dz;
    int ix, iy, iz;

    linecnt++;
    if (fgets(line, 128, instr) == NULL)
	error("in_vect in %s: unexpected EOF\n", getargv0());
    if (! hexpack) {
	if (sscanf(line, " %lf%lf%lf\n", &dx, &dy, &dz) != 3)
	    error("in_vect in %s: cannot read line %d: %s", 
	            getargv0(), linecnt, line);
	vptr[0] = dx;
	vptr[1] = dy;
	vptr[2] = dz;
    } else {
	if (sscanf(line, " %x%x%x\n", &ix, &iy, &iz) != 3)
	    error("in_vect in %s: cannot read line %d: %s", 
	            getargv0(), linecnt, line);
	vptr[0] = (2 * scale / hexmax) * ix - scale;
	vptr[1] = (2 * scale / hexmax) * iy - scale;
	vptr[2] = (2 * scale / hexmax) * iz - scale;
    }
}
