/* 
 * grid: aid for 1-D gridding
 *
 *      4-mar-94    ansi            pjt
 *     22-jan-95    fixed real=float problems
 *     20-jun-01    gcc3
 *     15-mar-05    proper real_proc for c++
 *
 */


#include <stdinc.h>
#include <grid.h>

void inil_grid(
    Grid *g,
    int n,
    real gmin,
    real gmax)
{
    if (n<1) error("inil_grid: n=%d",n);
    if (gmin==gmax) error("inil_grid: gmin=gmax=%g",gmin);

    g->mode = GRID_LINEAR;
    g->n = n;
    g->gmin = gmin;
    g->gmax = gmax;
    g->dg = (gmax-gmin)/n;
    g->up = gmin < gmax;
    g->g = NULL;
}

void inia_grid(
    Grid *g,
    int n,
    real *a)
{
    int i;

    if (n<2) error("inia_grid: n=%d",n);
    if (a[0]==a[1]) error("inia_grid: a(1)=a(2)=%g",a[0]);
    g->up = a[0] < a[1];
    for (i=1; i<n; i++) {
        if (a[i]==a[i-1]) error("inia_grid: a(%d)=a(%d)=%g",i,i+1,a[i]);
        if (g->up) {
            if (a[i-1] > a[i])
                error("inia_grid: array not monotonically increasing a(%d)",i+1);
        } else {
            if (a[i-1] < a[i])
	        error("inia_grid: array not monotonically decreasing a(%d)",i+1);
	}
    }
    g->mode = GRID_ARRAY;
    g->n = n;
    g->g = (real *) allocate(n*sizeof(real));
    for (i=0; i<n; i++) g->g[i] = a[i];
    g->gmin = a[0];
    g->gmax = a[n-1];
}

void inip_grid(
    Grid *g,
    int n,
    real fmin,
    real fmax,
    real_proc f)
{
    g->mode = GRID_PROC;
    g->n = n;
    g->gmin = fmin;
    g->gmax = fmax;
    g->dg = (fmax-fmin)/n;
    g->f = f;
    g->g = NULL;
}


int index_grid(
    Grid *g,
    real x)
{
    int imin, imax, idx = 0;

    switch (g->mode) {
      case GRID_PROC:
                        x = (*g->f)(x);
      case GRID_LINEAR:
      
                        if (x < g->gmin || x > g->gmax) return -1;
                        idx = (int) floor(  (x- g->gmin)/g->dg );
                        if (idx==g->n) idx--;   /* right edge exception */
                        break;
      case GRID_ARRAY:
      			if (g->up) {
	                    if (x < g->gmin || x > g->gmax) return -1;
	                } else {
	                    if (x < g->gmax || x > g->gmin) return -1;
			}
                        imin = 0;
                        imax = g->n;
                        error("grid_array not implemented yet");
                        break;
      default:
                        error("illegal grid mode %d",g->mode);
                        idx = -1;
                        break;
    }
    return idx;
}

real value_grid(
    Grid *g,
    real x)
{
    int idx;
    real f;
	
    switch (g->mode) {
      case GRID_PROC:
      case GRID_LINEAR:
            		return g->gmin + x * g->dg;
      case GRID_ARRAY:
			if (x<0) error("Cannot look extrapolate left");
			if (x> (real)g->n) error("Cannot look extrapolate right");
			idx = (int) floor(x);
			if (idx==g->n) return g->g[idx-1];	/* right edge */
#if 0
			f = remainder(x,1.0)/ABS(g->g[idx+1] - g->g[idx]);
#else
			f = 0.0;
			error("grid.c: need to implement remainder(real,real)");
#endif
			return f * g->g[idx] + (1-f)*g->g[idx+1];
      default:
                        error("illegal grid mode %d",g->mode);
    }
    return 0.0; /* never reached */
}



#ifdef TESTBED

#include <getparam.h>

string defv[] = {
	"linear=\n	N,xmin,xmax in case of linear grid",
	"array=\n	x1,x2,...xN in case of lookup table",
	"proc=\n	Procedure (not implemented)",
        "x=\n           Value(s) to return grid index(s) for",
	"i=\n		Index(s) to return grid value(s) for",
	"VERSION=1.0\n	5-nov-93 PJT",
	NULL,
};

#define MAXN  1000

nemo_main()
{
    real a[MAXN], x[MAXN];
    int n, i;
    Grid g;

    if (hasvalue("linear")) {
         if (nemoinpr(getparam("linear"),a,3) != 3) error("Linear Parsing");
         n = (int) a[0];
         inil_grid(&g,n,a[1],a[2]);
    } else if (hasvalue("array")) {
         n = nemoinpr(getparam("array"),a,MAXN);
         if (n<0) error("Array parsing");
         inia_grid(&g,n,a);
    } else if (hasvalue("proc"))
         error("proc= not implemented");
    else error("Need one of: linear=, array=, proc=");

    n = nemoinpr(getparam("x"),x,MAXN);
    for (i=0; i<n; i++)
    	printf("x=%g i=%d\n",x[i],index_grid(&g,x[i]));
    n = nemoinpr(getparam("i"),x,MAXN);
    for (i=0; i<n; i++)
    	printf("i=%g x=%g\n",x[i],value_grid(&g,x[i]));
}

#endif
