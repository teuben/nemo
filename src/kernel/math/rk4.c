/*
 *  RK4: take a single RK4 - step
 *
 *	20-jun-01	gcc3					pjt
 *
 *  Numerical Recipies, pp. xxx
 */
 
#define  NMAX   10


void rk4 (double *y, double *dydx, int n, double x, double h, double *yout, 
	int (*derivs)())
{
      double yt[NMAX], dyt[NMAX], dym[NMAX];
      double hh, h6, xh;
      int    i;
      
      hh=h*0.5;
      h6=h/6.0;
      xh=x+hh;
      for (i=0; i<n; i++)
        yt[i]=y[i]+hh*dydx[i];

      derivs (xh,yt,dyt);
      for (i=0; i<n; i++)
        yt[i]=y[i]+hh*dyt[i];

      derivs(xh,yt,dym);
      for (i=0; i<n; i++) {
        yt[i]=y[i]+h*dym[i];        
        dym[i]=dyt[i]+dym[i];
      }

      derivs(x+h,yt,dyt);
      for (i=0; i<n; i++)       
        yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.*dym[i]);

}
