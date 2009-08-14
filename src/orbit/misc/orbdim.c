/*
 * CR:	 determine number of isolating integrals for an orbit
 *	 using phase space correllation function
 *		8-jun-88 initial version?			pjt
 *		3-dec-93 nemo v2.x installed?			pjt
 *	       11-jan-95 linux declaration fix			pjt
 *  Tancredi, G., Sanchez, A., and Roig, F. (2000). 
 *  "A Comparison between Methods to Compute Lyapunov Exponents", 
 *  Astron. J. 121, 1171-1179.
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <orbit.h>

string defv[] = {
	"in=???\n			input orbit",
	"lograd=-3:0:0.2\n		which radii to look at",
	"xnorm=1,1,1\n			normalization of phase-space",
	"vnorm=1,1,1\n			-",
	"nsteps=\n			number of step (def: all)",
	"VERSION=1.1a\n,		11-jan-95 pjt",
	NULL,
};
string usage = "Determine the number of isolating integrals in an orbit";

string	infile;				/* file names */
stream  instr;				/* file streams */

#define MRAD 1000
real	lograd[MRAD];
real	logcr[MRAD];
int     nrad;

real	xnorm[NDIM];
real	vnorm[NDIM];

orbitptr optr;
int      nsteps;


nemo_main ()
{
	int i;
	
	setparams();

	instr = stropen (infile,"r");
	
	optr=NULL;
	read_orbit(instr,&optr);
	if (hasvalue("nsteps"))
		nsteps=getiparam("nsteps");
	else
		nsteps=Nsteps(optr);
	normalize();
	make_cr (nrad, lograd, logcr, PhasePath(optr), nsteps , 2*Ndim(optr));
	strclose(instr);
}

setparams()
{
	int nx,nv,i;
	
	infile = getparam("in");
	nrad = nemoinpr(getparam("lograd"), lograd, MRAD);
	if (nrad<1) error ("Parsing lograd=");
	nx = nemoinpr(getparam("xnorm"),xnorm,NDIM);
	if (nx!=NDIM)
		error ("xnorm must have length NDIM");
	nv = nemoinpr(getparam("vnorm"),vnorm,NDIM);
	if (nv!=NDIM)
		error ("vnorm must have length NDIM");
}

/* normalize phase space coordinates */

normalize ()
{
	int i;
	real invsqx[NDIM], invsqv[NDIM];
	
	for (i=0; i<NDIM; i++)
		invsqx[i] = 1.0/(xnorm[i]);
	for (i=0; i<NDIM; i++)
		invsqv[i] = 1.0/(vnorm[i]);
	
	for (i=0; i<Nsteps(optr); i++) {
		Xorb(optr,i) *= invsqx[0];
		Yorb(optr,i) *= invsqx[1];
		Zorb(optr,i) *= invsqx[2];
		Uorb(optr,i) *= invsqv[0];
		Vorb(optr,i) *= invsqv[1];
		Worb(optr,i) *= invsqv[2];
	}

	printf ("Phase space is normalized by %f %f %f %f %f %f\n",
		xnorm[0],xnorm[1],xnorm[2],vnorm[0],vnorm[1],vnorm[2]);
}
		

/* 
 *  CR: compute the phase space correlation function for
 *     'n' values of the 'radius'. These radii are stored
 *     as their logs in an array 'r'. On output the
 *     array 'cr' contains the values of the correlation
 *     at those corresponding radii. (also in log)
 *
 *	On input x must contain an array x[nsteps][ndim] of
 *	(normalized) phase space coordinates.
 *		See: Carnevali & R...  - 1988 ...
 */
 
make_cr(n, r, cr, x, nsteps, ndim)
int  n, ndim, nsteps;
real r[], cr[];
real *x;
{
    int  i,j,k;
    real dsq, cadd, dcdr;

    for (i=0; i<n; i++) {
/*	printf ("logr=%f ",r[i]);		*/
        r[i] = pow(10.0,2*r[i]);	/* make proper r^2 */
/*	printf (" r^2=%f\n",r[i]);		*/
        cr[i]=0.0;
    }
    printf ("n=%d nsteps=%d ndim=%d\n",n,nsteps,ndim);
/*
    for (i=0; i<nsteps; i++) {
      for (k=0; k<ndim; k++)
        printf ("%f ",*(x+i*ndim+k));
      printf ("\n");
    }
*/

    for (i=0; i<nsteps; i++)
        for (j=0; j<=i; j++) {
            cadd = ((i==j)?1.0:2.0);	
            dsq=0.0;
            for (k=0; k<ndim; k++)
		dsq += sqr(*(x+i*ndim+k) - *(x+j*ndim+k));
	               /*    x[i][k]           x[j][k]    */
            for (k=0; k<n; k++)
                if (dsq<r[k]) {
                    cr[k] += cadd;
                }
        }

    for (k=0; k<n; k++) {
        cr[k] = log10(cr[k]/(nsteps*nsteps));	/* return log(cr) */
        r[k] = 0.5*log10(r[k]);		/* return log(r) */
	printf (" %f %f ",r[k],cr[k]);
	if (k>0) {
		dcdr=(cr[k]-cr[k-1])/(r[k]-r[k-1]);
		printf ("%f\n",dcdr);
	} else
		printf ("\n");
    }
}
