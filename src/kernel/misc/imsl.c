/*
 *	interface to simulate some old ISML spline routines
 *		ICSCCU, ICSEVU, DCSEVU
 *	Calls Forsyth routines, which we have implemented in NEMO 
 *	as spline(), seval() slighly differently -  (see libJ.a in NEMO) 
 *
 *	C calling conventions
 *	  originally 'ier' is returned in the arguments, now it is
 *	  by the function
 *
 *	4-jun-88	created, needed for anisot.c	      Peter Teuben	
 *	9-sep-90	now C calling convention - for mkexphot        PJT
 *     25-feb-92        happy gcc2.0                                   PJT
 *     21-sep-93        ansi
 */
 
#include <stdinc.h>
#include <spline.h>

/*
 *	ICSCCU:		cubic spline interpolation (setting it up)
 */
 

void icsccu (real *x, real *y, int nx, real *c, int ic, int ier)		
{
    if (nx < ic) error("icsccu IMSL simulator: needs NX>=IC %d %d\n",nx,ic);
	
    spline (c, x, y, nx);	/* call Forsyth et al */
}

/*
 *	ICSEVU:		evaluation of a cubic spline
 */


void icsevu (real *x, real *y, int nx, real *c, int ic,
             real *u, real *s, int m, int ier)
{
    int i;

    if (nx < ic) error("icsevu IMSL simulator: needs NX>=IC %d %d\n",nx,ic);

    for (i=0; i<m; i++)
        s[i] = seval(u[i],x,y,c,nx);
}


/*
 *	DCSEVU:		cubic spline first and second derivative evaluator
 */

void dcsevu (real *x, real *y, int nx, real *c, int ic,
             real *u, real *ds, int m1, real *dss, int m2, int ier)
{
    int i;

    if (nx < ic) error("dcsevu IMSL simulator: needs NX>=IC %d %d\n",nx,ic);

    for (i=0; i<m1; i++)
        ds[i] = spldif(u[i],x,y,c,nx);
	
    for (i=0; i<m2; i++)
        dss[i] = spldif2(u[i],x,y,c,nx);
}

