#include <stdinc.h>

/*
 *  RKQC:  integrate a rk4 step ..
 *
 *	20-jun-01	gcc3				PJT
 */
 

#define NMAX	10
#define FCOR    0.0666666667
#define ONE	1.0
#define SAFETY  0.9
#define ERRCON  0.0006


extern void rk4 (double *, double *, int, double, double, double *, iproc);



void rkqc(double *y, double *dydx, int n,
	double *x, double htry, double eps, double *yscal,
	double *hdid, double *hnext,
	int (*derivs)())
{
      double ytemp[NMAX], ysav[NMAX], dysav[NMAX];
      double pgrow, pshrnk, xsav, h, hh, errmax;
      int    i;
      
      pgrow = -0.20;		/* RK4 */
      pshrnk= -0.25;		/* RK4 */
      xsav = *x;
      for (i=0; i<n; i++) {	/* save vector and derivative */
        ysav[i]=y[i];
        dysav[i]=dydx[i];
      }

      h=htry;			/* take this step to start with */

restep:  		/*  restep label in case more accuracy needed */

      hh=0.5*h;					/* Take TWO half steps */
      rk4(ysav,dysav,n,xsav,hh,ytemp,derivs);
      derivs(xsav+hh,ytemp,dydx);
      rk4(ytemp,dydx,n,xsav+hh,hh,y,derivs);	/*  Result in y[] */
      if((xsav+h) == xsav)
         printf ("RKQC: Stepsize not significant in RKQC.");

      rk4(ysav,dysav,n,xsav,h,ytemp,derivs);	/*  Take ONE whole step */
						/*  Result in ytemp[] */
      errmax=0.0;
      for (i=0; i<n; i++) {
/* printf ("errmax loop @%f %f %f %f ",xsav,y[i],ytemp[i],yscal[i]);   */
        ytemp[i] = y[i] - ytemp[i];		/* note sign */
        errmax=MAX(errmax,fabs(ytemp[i]/yscal[i]));
/* printf ("errmax/eps %d => %f\n",i,errmax/eps);  */
      }
      errmax=errmax/eps;	/* largest relative error */

      if(errmax > ONE) {
        h = SAFETY*h*pow(errmax,pshrnk);
/*  printf ("RKQC RESTEP h=%f\n",h);	*/
        goto restep; 
      }
      else {
        *hdid = h;
        if (errmax > ERRCON)
          *hnext = SAFETY*h*pow(errmax,pgrow);
        else
          *hnext = 4.0*h;
      }

      for (i=0; i<n; i++)            /* add higher order correction  */
        y[i]=y[i]+ytemp[i]*FCOR;      /* to the integration vector */
      *x += h;				/* and also update x !  */
}
