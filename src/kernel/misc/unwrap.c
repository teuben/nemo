/*
 * unwrap: routine that unwraps phases to a contiguous string
 *	22-jan-95	ansi prototypes
 */

#include <stdinc.h>
 
void unwrap(real w, int n, real *p, int idx, int mode, real fixval)
{
    int i, j;
    real d, x;

    if (n<2) return;
    x = p[0];
    for (i=0; i<n; i++) {
        d = (p[i]-x)/w;         /* diff with previous point */
        if (d > 0.5)
            for (j=i; j<n; j++) p[j] -= w;
        else if (d < -0.5)
            for (j=i; j<n; j++) p[j] += w;
	x = p[i];		/* remember last point */
    }
}

#ifdef TESTBED

#include <getparam.h>

string defv[] = {
    "in=???\n       Table with phases in column 1",
    "out=???\n      Output table with unwrapped phases",
    "period=360\n   Period of phases",
    "xcol=1\n       Column of the phases to unwrap",
    "mean=\n        Align around this as mean",
    "min=\n         Minimum phase should be above this",
    "max=\n         Maximum phase should be below this",
    "fix=1\n        Fix phase around this array element (1..n)",
    "wrap=f\n       Output wrap too?",	/* don't work */
    "nmax=10000\n   Max table allocation of lines",
    "VERSION=0.1b\n 9-nov-93 PJT",
    NULL,
};

string usage="unwrap phases";

nemo_main()
{
    real *coldat[1], *p, period, psum=0.0;
    int i, n, nmax, colnr[1];
    bool Qwrap;
    stream instr, outstr;

    Qwrap = getbparam("wrap");
    period = getdparam("period");
    instr = stropen(getparam("in"),"r");
    outstr = stropen(getparam("out"),"w");

    nmax = file_lines(getparam("in"),getiparam("nmax"));
    if (nmax<1) error("Problem reading input file");

        
    colnr[0] = getiparam("xcol");
    p = coldat[0] = (real *) allocate(nmax*sizeof(real));
    n = get_atable(instr,1,colnr,coldat, nmax);
    if (n<1) error("get_atable returns with %d data",n);

    unwrap(period, n, p, -1, 0, 0.0);

    for (i=0; i<n; i++)
	fprintf(outstr,"%g\n",p[i]);
}

#endif
