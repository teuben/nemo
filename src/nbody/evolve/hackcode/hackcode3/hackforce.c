/*
 * HACKFORCE.C: hackcode force calculation tool.
 *
 *   updates:
 *	7-jul-89  V1.1: added mass to options - interchanged order test,out PJT
 *     19-jan-90  V1.2: fixed typo, resulted in not enuf memory for Sparc
 *			also multi-snapshot in-file allowed (no times= yet)
 *			If test is blank, the time-looping is done in in=
 *			if test is not blank, time-looping applied to test=
 */

#include "defs.h"
#include <getparam.h>
#include <filestruct.h>
#include <snapshot.h>

string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "in=???\nInput file with mass coordinates",
    "out=\nOutput file with f.c. results",
    "test=\nInput file with test coordinates (if blank, mass coords are used)",
    "tol=1.0\nCell subdivision tolerence",
    "eps=0.05\nStandard softening parameter",
    "rsize=4.0\nSide-length of initial box",
    "rmin=\nLower left corner of initial box [default is -rsize/2 (centered)",
    "options=mass,phase\nOutput options: phase and/or mass",
    "fcells=0.75\nCell/body allocation ratio",
    "VERSION=1.2\nPJT 19-jan-90",
    NULL,
};

main(argc, argv)
int argc;
string argv[];
{
    initparam(argv, defv);
    while (inputdata())	{			/* input mass and test data */
       forcecalc();				/* find force at test pos */
       outresult();				/* write snap with results */
    }
}

bodyptr massdata;	/* array of mass points */
int nmass;		/* number of mass points */

bodyptr testdata;	/* array of test points */
int ntest;		/* number of test points */

stream instr=NULL;	/* input file for masses */
stream tststr=NULL;	/* file for which force calc done (def: in-file */

inputdata()
{
    string input, test;
    int    i;

    if (instr == NULL) {        /* initialization only at first entry */
        input = getparam("in");             /* get filename */
        instr = stropen(input, "r");        /* open file */
        test = getparam("test");            /* get name for testfile */
        if (*test != NULL)                  /* if testfile provided: */
	    tststr = stropen(test, "r");    /* ... open it */
    }

    if (tststr==NULL) {                     /* if no testfile, data from in */
        get_history(instr);
        i = readsnapshot(&massdata, &nmass, instr); /* read mass coord data */
        if (i==0) return(0);
	testdata = massdata;			/* use mass data for tests */
	ntest = nmass;
    } else {                                /* else data from test */
        if (massdata==NULL) {                   /* on first pass read masses */
            get_history(instr);            
            readsnapshot(&massdata, &nmass, instr);
            if (i=0) return(0);
        }
	get_history(tststr);                    /* read (next) testdata */
	readsnapshot(&testdata, &ntest, tststr);
        if (i=0) return(0);
    }
    return(1);
}

readsnapshot(btab_ptr, nobj_ptr, instr)
bodyptr *btab_ptr;				/* gets particle array */
int *nobj_ptr;					/* gets number of bodies */
stream instr;					/* stream to read from */
{
  int nobj, cs, i;
  real *mbuf, *mp, *pbuf, *pp;
  bodyptr bp;

  for(;;) {                        /* loop until done or proper snapshot */
    if (!get_tag_ok(instr,SnapShotTag))
        return(0);
    get_set(instr, SnapShotTag);
    if (!get_tag_ok(instr,ParametersTag))
        continue;
    get_set(instr, ParametersTag);
    get_data(instr, NobjTag, IntType, &nobj, 0);
    if (nobj < 1)
	error("readsnapshot: %s = %d  is absurd\n", NobjTag, nobj);
    get_tes(instr, ParametersTag);
    if (!get_tag_ok(instr,ParticlesTag))
        continue;
    get_set(instr, ParticlesTag);
    get_data(instr, CoordSystemTag, IntType, &cs, 0);
    if (cs != CSCode(Cartesian, NDIM, 2))
	error("readsnapshot: cannot handle %s = %d\n", CoordSystemTag, cs);
    mbuf = mp = (real *) allocate(nobj * sizeof(real));
    pbuf = pp = (real *) allocate(nobj * 2 * NDIM * sizeof(real));
    if (mp == NULL || pp == NULL)
	error("readsnapshot: not enuf memory for buffers\n");
    get_data(instr, MassTag, RealType, mbuf, nobj, 0);
    get_data(instr, PhaseSpaceTag, RealType, pbuf, nobj, 2, NDIM, 0);
    get_tes(instr, ParticlesTag);
    get_tes(instr, SnapShotTag);
    *btab_ptr = bp = (bodyptr) allocate(nobj * sizeof(body));
    if (bp == NULL)
	error("readsnapshot: not enuf memory for bodies\n");
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
    return(1);
  }
}

real *phidata, *accdata;	/* per-test-pos potential, acceleration */

int n2btot, nbctot;		/* body-body, body-cell interactions */

real cputree, cpufcal;		/* CPU time to build tree, compute forces */

forcecalc()
{
    real *pp, *ap;
    double cputime(), cpubase, atof();
    string *burststring(), *rminxstr;
    int xstrlen(), i;
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
    printf("initial rsize: %8f    rmin: %8f  %8f  %8f\n",
	   rsize, rmin[0], rmin[1], rmin[2]);
    fcells = getdparam("fcells");
    phidata = pp = (real *) allocate(ntest * sizeof(real));
    accdata = ap = (real *) allocate(ntest * NDIM * sizeof(real));
    if (pp == NULL || ap == NULL)
	error("forcecalc: not enuf memory for results\n");
    cpubase = cputime();
    maketree(massdata, nmass);
    cputree = cputime() - cpubase;
    printf("  final rsize: %8f    rmin: %8f  %8f  %8f\n",
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

outresult()
{
    string out;

    if (outstr==NULL) {
        out = getparam("out");
	outstr = stropen(out, "w");
	put_history(outstr);
    }
    writesnapshot();			/* output testdata results */
}

writesnapshot()
{
    real *mbuf, *mp, *pspbuf, *pspp;
    bodyptr bp;
    int cs = CSCode(Cartesian, NDIM, 2);
    string options = getparam("options");
    bool scanopt();

    mbuf = mp = (real *) allocate(ntest * sizeof(real));
    pspbuf = pspp = (real *) allocate(ntest * 2 * NDIM * sizeof(real));
    if (mbuf == NULL || pspbuf == NULL)
	error("writesnapshot: not enuf memory for buffers\n");
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
