/*
 * HACKDENS.C: hackcode local density calculation tool.
 *
 *          1988  V1.0   originally written        Jun Makino
 *	7-jul-89  V1.1   updated with get_history()	PJT
 *     23-oct-90  V1.2   helpvec			PJT
 *     18-jul-92  V1.3   replaced many if(debug)printf(...) by dprintf(1,...)
 *			 Also nemo_main() and usage	PJT
 *	1-apr-01      a  compiler warning
 *     24-may-02  V2.0   fixed for high-N system by using int_hack in load.c
 *     29-may-02  V2.1   add nudge= keyword to nudge overlapping particles
 *     25-apr-06  V2.2b  use global to isolate extern's (for Mac linking)
 *     28-jul-06  V2.2c  default for tag is now Density
 *                V2.2d  clarify D vs. P, working with std snapshot, not archaic
 *     21-dep-23  V2.3   add a slow direct= for benchmark/comparison
 *     12-oct-23  V2.4   scale by mass
 *     11-OCT-23  v3.0   add norm=1 and made it the default
 *
 * NOTE:   for snapshots with unequal masses this program doesn't work
 *
 */

#define global

#include "defs.h"
#include <getparam.h>
#include <history.h>
#include <filestruct.h>
#include <snapshot/snapshot.h>

string defv[] = {	
    "in=???\n			  input snapshot",
    "out=\n			  optional output file with density results ",
    "neib=6\n			  number of neighbours to define local density ",
    "rneib=0.1\n		  initial guess for neighbour sphere radius ",
    "write_at_phi=f\n		  flag to write density with Potential instead of Density tag",
    "rsize=4.0\n		  side-length of initial box",
    "rmin=\n			  lower left corner of initial box",
    "options=phase,mass\n	  misc. control options {phase, mass}",
    "fcells=0.9\n		  cell/body allocation ratio ",
    "nudge=0\n                    nudge overlapping particles with this dispersion",
    "verbose=f\n		  flag to print # of particles finished ",
    "density=t\n                  write density, or distance to Kth particle",
    "ndim=3\n                     3D or 2D computation",
    "direct=f\n                   slower direct density computation",
    "norm=1\n                     normalization mode (0=nothing   1=1/N)",
    "VERSION=3.0\n		  12-oct-2023 PJT",
    NULL,
};

string usage = "hackcode local density calculation tool";


void inputdata(void);
void readsnapshot(bodyptr *btab_ptr, int *nobj_ptr, stream instr);
void dencalc(void);
void outresult(void);
void writesnapshot(void);

// util.c
extern void pickvec(vector x, bool cf);

// load.c
extern void maketree(bodyptr btab, int nbody, double nudge);


void nemo_main()
{
    inputdata();			/* input mass and test data */
    dencalc();				/* find force at test pos */
    outresult();			/* write snap with results */
}

bodyptr massdata;		/* array of mass points */
int nmass;			/* number of mass points */
real totalmass;

bodyptr testdata;		/* array of test points */
int ntest;			/* number of test points */

bool Qdirect;
bool norm;

void inputdata()
{
    string input;
    stream instr;

    input = getparam("in");
    instr = stropen(input, "r");
    Qdensity = getbparam("density");
    Qdirect = getbparam("direct");
    readsnapshot(&massdata, &nmass, instr);	/* read mass coord data */
    strclose(instr);
    testdata = massdata;
    ntest = nmass;
    norm = getiparam("norm");
}

real tsnap;

void readsnapshot(btab_ptr, nobj_ptr, instr)
bodyptr *btab_ptr;				/* gets particle array */
int *nobj_ptr;					/* gets number of bodies */
stream instr;					/* stream to read from */
{
    int nobj, cs, i;
    real *mbuf, *mp, *pbuf, *pp;
    bodyptr bp;

    get_history(instr);

    get_set(instr, SnapShotTag);
    get_set(instr, ParametersTag);
    if (get_tag_ok(instr, TimeTag)) {
	get_data_coerced(instr, TimeTag, RealType, &tsnap, 0);
    }else{
	tsnap=0.0;
    }
    get_data(instr, NobjTag, IntType, &nobj, 0);
    if (nobj < 1)
	error("readsnapshot: %s = %d  is absurd\n", NobjTag, nobj);
    get_tes(instr, ParametersTag);
    get_set(instr, ParticlesTag);
    get_data(instr, CoordSystemTag, IntType, &cs, 0);
    if (cs != CSCode(Cartesian, NDIM, 2))
	error("readsnapshot: cannot handle %s = %d\n", CoordSystemTag, cs);
    mbuf = mp = (real *) malloc(nobj * sizeof(real));
    pbuf = pp = (real *) malloc(nobj * 2 * NDIM * sizeof(real));
    if (mp == NULL || pp == NULL)
	error("readsnapshot: not enuf memory for buffers\n");
    get_data_coerced(instr, MassTag, RealType, mbuf, nobj, 0);
    get_data_coerced(instr, PhaseSpaceTag, RealType, pbuf, nobj, 2, NDIM, 0);
    get_tes(instr, ParticlesTag);
    get_tes(instr, SnapShotTag);
    *btab_ptr = bp = (bodyptr) malloc(nobj * sizeof(body));
    if (bp == NULL)
	error("readsnapshot: not enuf memory for bodies\n");
    totalmass = 0.0;
    for (i = 0; i < nobj; i++) {
	Type(bp) = BODY;
	Mass(bp) = *mp++;
	SETV(Pos(bp), pp);
	pp += NDIM;
	SETV(Vel(bp), pp);
	pp += NDIM;
	bp++;
	totalmass += Mass(bp);
    }
    dprintf(0,"Total mass = %g\n", totalmass);
    free(mbuf);
    free(pbuf);
    *nobj_ptr = nobj;
}

real *dendata;			/* local density */

int n2btot, nbctot;		/* body-body, body-cell interactions */

real cputree, cpufcal;		/* CPU time to build tree, compute forces */

void dencalc()
{
    real hackden(), directden();
    real *pp, *work, rneib, newrneib, nudge;
    int neibnum;
    //double cputime(), cpubase;
    double cpubase;
      // , atof();
    // string *burststring(), *rminxstr;
    string *rminxstr;
    int xstrlen(), i, ibody;
    bodyptr bp;
    bool verbose;

    verbose=getbparam("verbose");
    rneib=getdparam("rneib");
    neibnum=getiparam("neib")+1;
    nudge = getdparam("nudge");
    if (nudge > 0) {
      set_xrandom(0);   /* should use seed= */
#if 0
      warning("Nudging all particles by +/-%s",nudge);
      for (bp=massdata; bp<massdata+nmass; bp++) {
	for (i=0; i<NDIM; i++)
	  Pos(bp)[i] += xrandom(-nudge,nudge);
      }
#endif
    }
#ifdef DEBUG
    dprintf(0,"neib=%d rneib=%f\n", neibnum, rneib);
#endif    
    rsize = getdparam("rsize");
    rminxstr = burststring(getparam("rmin"), ", ");
    if (xstrlen(rminxstr, sizeof(string)) < NDIM) {
	SETVS(rmin, - rsize / 2.0);
    } else
	for (i = 0; i < NDIM; i++)
	    rmin[i] = atof(rminxstr[i]);
    dprintf(0,"initial rsize: %8f    rmin: %8f  %8f  %8f\n",
	   rsize, rmin[0], rmin[1], rmin[2]);
    fcells = getdparam("fcells");
    dendata = pp = (real *) malloc(ntest * sizeof(real));
    work = (real *) malloc(ntest * sizeof(real));
    if (pp == NULL || work==NULL)
	error("forcecalc: not enuf memory for results");
    cpubase = cputime();
    maketree(massdata, nmass,nudge);
    cputree = cputime() - cpubase;
    dprintf(0,"  final rsize: %8f    rmin: %8f  %8f  %8f\n",
	   rsize, rmin[0], rmin[1], rmin[2]);
    cpubase = cputime();
    n2btot = nbctot = 0;
    ibody=0;
    for (bp = testdata; bp < testdata+ntest; bp++) {
	if (Qdirect) 
	  *pp = directden(bp, neibnum, rneib, work, testdata, ntest);
	else
	  *pp = hackden(bp, neibnum, rneib, &newrneib, work);
	rneib=0.95*rneib+0.05*newrneib;
	pp++;
	ibody++;
	if(verbose && ibody%100==0)dprintf(0," %d rn=%f\n", ibody,rneib);
    }
    cpufcal = cputime() - cpubase;
    if (norm==1) {
      dprintf(1,"Renormalizing %d densities\n",ntest);
      real factor = 1.0/ntest;
      if (Qdensity) factor *= totalmass;
      for (pp=dendata; pp < dendata+ntest; pp++) 
	*pp *= factor;
    }
}

stream outstr;

void outresult()
{
    string out;

    out = getparam("out");
    if (*out != 0) {
	outstr = stropen(out, "w");
	put_history(outstr);
	writesnapshot();			/* output testdata results */
	strclose(outstr);
    }
}

void writesnapshot()
{
    real *mbuf, *mp, *pspbuf, *pspp;
    bodyptr bp;
    int cs = CSCode(Cartesian, NDIM, 2);
    string options = getparam("options");
    bool scanopt();

    mbuf = mp = (real *) malloc(ntest * sizeof(real));
    pspbuf = pspp = (real *) malloc(ntest * 2 * NDIM * sizeof(real));
    if (mbuf == NULL || pspbuf == NULL)
	error("writesnapshot: not enuf memory for particle buffers");
    for (bp = testdata; bp < testdata+ntest; bp++) {
	*mp++ = Mass(bp);
	SETV(pspp, Pos(bp));
	pspp += NDIM;
	SETV(pspp, Vel(bp));
	pspp += NDIM;
    }
    put_set(outstr, SnapShotTag);
    put_set(outstr, ParametersTag);
    put_data(outstr, NobjTag, IntType, &ntest, 0);
    put_data(outstr, TimeTag, RealType, &tsnap, 0);
    put_tes(outstr, ParametersTag);
    put_set(outstr, ParticlesTag);
    put_data(outstr, CoordSystemTag, IntType, &cs, 0);
    if (scanopt(options, "mass"))
	put_data(outstr, MassTag, RealType, mbuf, ntest, 0);
    if (scanopt(options, "phase"))
	put_data(outstr, PhaseSpaceTag, RealType, pspbuf, ntest, 2, NDIM, 0);
    put_data(outstr, 
	     getbparam("write_at_phi") ? PotentialTag : DensityTag,
	     RealType, dendata, ntest, 0);
    put_tes(outstr, ParticlesTag);
    put_set(outstr, DiagnosticsTag);
    put_data(outstr, "cputree", RealType, &cputree, 0);
    put_data(outstr, "cpucal", RealType, &cpufcal, 0);
    put_tes(outstr, DiagnosticsTag);
    put_tes(outstr, SnapShotTag);
    free(mbuf);
    free(pspbuf);
}
