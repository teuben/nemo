/*
 * HACKFORCE.C: hackcode force calculation tool.
 *
 *   updates:
 *	7-jul-89  V1.1: added mass to options - interchanged order test,out PJT
 *     19-jan-90  V1.2: fixed typo, resulted in not enuf memory for Sparc
 *			also multi-snapshot in-file allowed (no times= yet)
 *			If test is blank, the time-looping is done in in=
 *			if test is not blank, time-looping applied to test=
 *	8-jul-90  V1.3  Using <snapshot> macros now			PJT
 *                      --- but abandoned again ---
 *     19-aug-90  V1.4  also transfers time as a parameter (Makino)     PJT
 *     29-oct-90  V1.4a -- but fixed bug created in previous version    PJT
 *     13-feb-91  V1.4b -- fixed yet another bug in reading test=       PJT
 *			also renamed 'time' to 'tsnap' to avoid time()  
 *     15-apr-92      c -- malloc() -> allocate()    PJT
 *     22-jul-93  V1.5  fixed bug when handling multisnapshots      PJT/JAN
 *	7-aug-94  V1.5a declaration of atof() fails on macro-versions (linux)
 *     20-sep-01      b NULL -> 0
 *     29-mar-04  V1.6  using 'global' macro to prevent mu;ltiple definitons
 */

#define global                                  /* don't default to extern  */
#include "code.h"
#include <getparam.h>
#include <filestruct.h>
#include <history.h>
#include <extstring.h>

#if 0
/*	Can't be done yet - Body etc. was already used in defs.h */
#include <snapshot/body.h>
#include <snapshot/snapshot.h>
#else
#include <archaic/snapshot.h>
#endif

string defv[] = {	
    "in=???\n       Input (snapshot) file with mass coordinates",
    "out=???\n      Output (snapshot) file with f.c. results",
    "test=\n        Input file with test coordinates (if blank, in= coords are used)",
    "tol=1.0\n      Cell subdivision tolerence",
    "eps=0.05\n     Standard softening parameter",
    "rsize=4.0\n    Side-length of initial box",
    "rmin=\n              Lower left corner of initial box [default is -rsize/2 (centered)",
    "options=mass,phase\n Output options: phase and/or mass",
    "fcells=0.75\n        Cell/body allocation ratio",
    "VERSION=1.6a\n       28-mar-05 PJT",
    NULL,
};

string usage="Add hackcode1 computed forces and potential to a snapshot";

extern double cputime(void);
extern string *burststring(string, string);
extern bool scanopt(string, string);


void nemo_main(void)
{
    while (input_data())	{			/* input mass and test data */
       force_calc();				/* find force at test pos */
       out_result();				/* write snap with results */
    }
}


static bodyptr massdata;	/* array of mass points */
static int nmass;		/* number of mass points */

static bodyptr testdata;	/* array of test points */
static int ntest;		/* number of test points */

static real tsnap;              /* some time that was obtained from input/test */

static stream instr=NULL;	/* input file for masses */
static stream tststr=NULL;	/* file for which force calc done (def: in-file */

int input_data(void)
{
    string input, test;
    int    i;

    if (instr == NULL) {        /* initialization only at first entry */
        input = getparam("in");             /* get filename */
        instr = stropen(input, "r");        /* open file */
        test = getparam("test");            /* get name for testfile */
        if (*test != 0)                     /* if testfile provided: */
	    tststr = stropen(test, "r");    /* ... open it */
    }

    if (tststr==NULL) {                     /* if no testfile, data from in */
        get_history(instr);
        i = read_snapshot(&massdata, &nmass, instr); /* read mass coord data */
        if (i==0) return 0;
	testdata = massdata;			/* use mass data for tests */
	ntest = nmass;
    } else {                                /* else data from test */
        if (massdata==NULL) {                   /* on first pass read masses */
            get_history(instr);            
            i=read_snapshot(&massdata, &nmass, instr);
            if (i==0) return 0;
        }
	get_history(tststr);                    /* read (next) testdata */
	i=read_snapshot(&testdata, &ntest, tststr);
        if (i==0) return(0);
    }
    return 1;
}

int read_snapshot(
		   bodyptr *btab_ptr,	     /* gets particle array */
		   int *nobj_ptr,	     /* gets number of bodies */
		   stream instr              /* stream to read from */
		   )
{
  int nobj, cs, i;
  real *mbuf, *mp, *pbuf, *pp;
  bodyptr bp;

  for(;;) {                        /* loop until done or proper snapshot */
    if (!get_tag_ok(instr,SnapShotTag))
        return 0;
    get_set(instr, SnapShotTag);
    if (!get_tag_ok(instr,ParametersTag)) {
    	get_tes(instr,SnapShotTag);
        continue;
    }
    get_set(instr, ParametersTag);
    get_data(instr, NobjTag, IntType, &nobj, 0);
    if (nobj < 1)
	error("read_snapshot: %s = %d  is absurd", NobjTag, nobj);
    if (get_tag_ok(instr,TimeTag))
        get_data(instr,TimeTag, RealType, &tsnap, 0);
    else {
        dprintf(0,"No time tag: time=0.0 assumed\n");
        tsnap = 0.0;
    }
    get_tes(instr, ParametersTag);
    if (!get_tag_ok(instr,ParticlesTag)) {
    	get_tes(instr,SnapShotTag);
        continue;
    }
    get_set(instr, ParticlesTag);
    get_data(instr, CoordSystemTag, IntType, &cs, 0);
    if (cs != CSCode(Cartesian, NDIM, 2))
	error("read_snapshot: cannot handle %s = %d", CoordSystemTag, cs);
    mbuf = mp = (real *) allocate(nobj * sizeof(real));
    pbuf = pp = (real *) allocate(nobj * 2 * NDIM * sizeof(real));
    get_data(instr, MassTag, RealType, mbuf, nobj, 0);
    get_data(instr, PhaseSpaceTag, RealType, pbuf, nobj, 2, NDIM, 0);
    get_tes(instr, ParticlesTag);
    get_tes(instr, SnapShotTag);
    *btab_ptr = bp = (bodyptr) allocate(nobj * sizeof(body));
    for (i = 0; i < nobj; i++) {
	Type(bp) = BODY;
	Mass(bp) = *mp++;
	SETV(Pos(bp), pp);
	pp += NDIM;
	SETV(Vel(bp), pp);
	pp += NDIM;
	bp++;
    }
    free(mbuf);
    free(pbuf);
    *nobj_ptr = nobj;
    return 1;
  }
}

real *phidata, *accdata;	/* per-test-pos potential, acceleration */

int n2btot, nbctot;		/* body-body, body-cell interactions */

real cputree, cpufcal;		/* CPU time to build tree, compute forces */

void force_calc(void)
{
    real *pp, *ap;
    double cpubase;
    string *rminxstr;
    int i;
    bodyptr bp;

    tol = getdparam("tol");
    eps = getdparam("eps");
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
    phidata = pp = (real *) allocate(ntest * sizeof(real));
    accdata = ap = (real *) allocate(ntest * NDIM * sizeof(real));
    cpubase = cputime();
    maketree(massdata, nmass);
    cputree = cputime() - cpubase;
    dprintf(0,"  final rsize: %8f    rmin: %8f  %8f  %8f\n",
	   rsize, rmin[0], rmin[1], rmin[2]);
    cpubase = cputime();
    n2btot = nbctot = 0;
    for (bp = testdata; bp < testdata+ntest; bp++) {
	hackgrav(bp);
	*pp++ = Phi(bp);
	SETV(ap, Acc(bp));
	ap += NDIM;
	n2btot += n2bterm;
	nbctot += nbcterm;
    }
    cpufcal = cputime() - cpubase;
}

stream outstr=NULL;

void out_result(void)
{
    string out;

    if (outstr==NULL) {
        out = getparam("out");
	outstr = stropen(out, "w");
	put_history(outstr);
    }
    write_snapshot();			/* output testdata results */
}

void write_snapshot(void)
{
    real *mbuf, *mp, *pspbuf, *pspp;
    bodyptr bp;
    int cs = CSCode(Cartesian, NDIM, 2);
    string options = getparam("options");

    mbuf = mp = (real *) allocate(ntest * sizeof(real));
    pspbuf = pspp = (real *) allocate(ntest * 2 * NDIM * sizeof(real));
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
    put_data(outstr, "tol", RealType, &tol, 0);
    put_data(outstr, "eps", RealType, &eps, 0);
    put_tes(outstr, ParametersTag);
    put_set(outstr, ParticlesTag);
    put_data(outstr, CoordSystemTag, IntType, &cs, 0);
    if (scanopt(options, "mass"))
	put_data(outstr, MassTag, RealType, mbuf, ntest, 0);
    if (scanopt(options, "phase"))
	put_data(outstr, PhaseSpaceTag, RealType, pspbuf, ntest, 2, NDIM, 0);
    put_data(outstr, PotentialTag, RealType, phidata, ntest, 0);
    put_data(outstr, AccelerationTag, RealType, accdata, ntest, NDIM, 0);
    put_tes(outstr, ParticlesTag);
    put_set(outstr, DiagnosticsTag);
    put_data(outstr, "n2btot", IntType, &n2btot, 0);
    put_data(outstr, "nbctot", IntType, &nbctot, 0);
    put_data(outstr, "cputree", RealType, &cputree, 0);
    put_data(outstr, "cpufcal", RealType, &cpufcal, 0);
    put_tes(outstr, DiagnosticsTag);
    put_tes(outstr, SnapShotTag);
    free(mbuf);
    free(pspbuf);
}
