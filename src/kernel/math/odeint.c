/*
 *  ODEINT:  integrate o.d.e. 
 *      15-mar-90:   made it return 0 consistently PJT
 *	16-feb-97:   no more nexted extern's       PJT
 *      24-nov-03:   NEMOfied, prototype for gcc3  PJT
 *      23-jan-04:   prototypes, nemofied          pjt
 */

#include <stdinc.h>

#define MAXSTP  10000		/* maximum steps to be taken */
#define NMAX	10		/* maximum number of independant variables */
#define TWO     2.0
#define ZERO    0.0
#define TINY    1.E-30		/* to prevent underflow */

#define sign(x,s) (((s)>0) ? (x) : (-(x)))	/* fortran definition */
 
local  int     kmax, kount;			/* used for history record */
local  double  dxsav, xp[100], yp[10][100];	/* of the integration */
      /* kmax (and dxsav if kmax>0) must be initialized in calling routine */


int odeint (ystart, nvar, x1, x2, eps, h1, hmin, nok, nbad, derivs, rkqc)
double ystart[];			/* starting vector, will get result */
int    nvar;				/* length of vector */
double x1,x2,eps,h1,hmin;		/* some more input conditions */
int    *nok,*nbad;			/* return codes */
int    (*derivs)(), (*rkqc)();		/* pointers to functions */
{
      double yscal[NMAX], y[NMAX], dydx[NMAX];
      double x, h, xsav;
      int    i, nstp;
      double hdid, hnext;

      dprintf(2,"ODEINT: x1,x2,eps,h1,hmin=%f %f %f %f %f\n",x1,x2,eps,h1,hmin);
      x=x1;			/* start integrating here */
      h=sign(h1,x2-x1);		/* with this initial stepsize */
      dprintf(2," ...  starts at x=%f   with h=%f\n",x,h);
      *nok=0;			/* diagnostics: # good steps */
      *nbad=0;			/* .. and bad steps (stepsize decreased) */
      kount=0;			/* kounter in case to save intermediates */
      for (i=0; i<nvar; i++)
	y[i]=ystart[i];		/* copy initial vector in working area y[] */
      xsav=x-dxsav*TWO;		/* make sure first point is saved */
      for (nstp=0; nstp<MAXSTP; nstp++) {
	derivs(x,y,dydx);
	for (i=0; i<nvar; i++) 
          yscal[i]=fabs(y[i])+fabs(h*dydx[i])+TINY;	/* calibration */
	dprintf(2,"ODEINT: y,h,dydx,yscal=%f %f %f %f\n",y[0],h,dydx[0],yscal[0]);
        if (kmax>0) {				/* history record */
          if (fabs(x-xsav) > fabs(dxsav)) {	/* of KMAX steps */
            if (kount < kmax){			/* with steps DXSAV */
              xp[kount]=x;
              for (i=0; i<nvar; i++)
                yp[i][kount]=y[i];
	      kount++;
              xsav=x;				/* next x to be saved */
            }
          }
        }

        if ((x+h-x2)*(x+h-x1) > ZERO) 
           h=x2-x;			/* last step ? */
           
        dprintf(2,"STEP %d x=%f h=%f\n",nstp,x,h);
        rkqc (y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);    /* step */
        				/* note the x got increased inside */

        dprintf(2," hdid=%f hnext=%f\n",hdid,hnext);
        if (hdid==h) 
          *nok += 1;
        else
          *nbad += 1;
        
        if ((x-x2)*(x2-x1) >= ZERO) {		/* OK, enough integrated */
          for (i=0; i<nvar; i++)
            ystart[i]=y[i];

          if (kmax > 0) {		/* add last point to history record */
            xp[kount]=x;
            for (i=0; i<nvar; i++)
              yp[i][kount]=y[i];
            kount++;
          }
          return 0;
        }
        if (fabs(hnext) < hmin) 
           warning("ODEINT: Stepsize smaller than minimum.");
        h=hnext;
      }
      warning("ODEINT: Too many steps.\n");
      return 0;
}
