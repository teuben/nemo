#define max(x,y)  ((x)>(y) ? (x) : (y) ) 

/*
 *  RKQS:  a better way (than rkqc) to integrate a rk4 step (NumRec2-p719)
 *   
 *  A firth-order Runge-Kutta step with monitoring of local truncation error
 * to ensure accuracy and adjust stepsize. Input are the dependant vector,
 * y[1..n], and its derivatives, dydx[1..n], at the starting value of the
 * independant variable x. Also input are the stepsize to be attempted
 * htry, the required accuracy eps, and the vector yscal[1..n] against
 * which the error is scaled. On output, y and x are replaced with their new
 * values, hdid is the stepsize that was actually accomplished, and hnext 
 * is the user-supplied routine that computes the right-hand side 
 * derivatives.
 */
 

#define NMAX	10  /* pjt */
#define PGROW   -0.2
#define PSHRNK  -0.25
#define ERRCON  1.89e-4
#define SAFETY  0.9


rkqc(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs) 
double  y[],dydx[],htry,eps,yscal[];
double  *x, *hdid, *hnext;				/* output (as well) */
int	n;
int	(*derivs)();
{
    double ytemp[NMAX], yerr[NMAX];
    double h, errmax, xnew;
    double pow(), fabs();
    int    rk4();
    int    i;
      

    h=htry;			/* take this step to start with */

    for(;;) {                   /* loop untilo happy */
      rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);        /* take a step */

      errmax=0.0;                           /* evaluate accuracy */
      for (i=0; i<n; i++)
        errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
      errmax /= eps;                  	/* largest relative error */

      if(errmax > 1.0) {             /* large truncation error, reduce step */
        h = SAFETY*h*pow(errmax,PSHRNK);
        if (h < 0.1*h) h *= 0.1;            /* no more than factor 10 */
        xnew = (*x) + h;
        if (xnew == *x) error("RKCK: stepsise underflow");
        continue;
      }
      else {                        /* step ok, compute size of next step */
        if (errmax > ERRCON)
          *hnext = SAFETY*h*pow(errmax,PGROW);
        else
          *hnext = 5.0*h;               /* no more than factor 5 increase */
        *x += (*hdid=h);
        for (i=0; i<n; i++)
          y[i] = ytemp[i];
        break;
      }
    }
}
