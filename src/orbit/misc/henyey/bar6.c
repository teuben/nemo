/* BAR6: prolate ferrers k=2 for HENYEY */
extern double sqrt(double), log(double);         /* externals used in here */

double sqr(double x)
{ 
  return x*x; 
}


#include "real.h"

static real Grav_Const = 1.0;        /* this makes you have to enter GM instead */
static int  first = 1;


/*
 *    BAR6:    the k=2 Ferrers prolate bar
 *             this one is for 2D though (sep-91)
 *             used to return second derivatives... for HENYEY...
 *             input: elipm   mass of bar
 *                      a,b   axes of bar  (a>b)
 *                      x,y   position of point of interest
 *                    ax,ay   forces [-d(potential)/dr]
 *                      ept   potential
 *              axx,ayy,azz   second derivatives of potential - for testing
 */
void bar6 (real elipm, real a, real b, real x, real y,                      /* input */
	   real *ept, real *ax, real *ay, real *axx, real *axy, real *ayy)  /* output */
{
    real aa,bb,xx,yy,zz,rr,gi,side,em,ee,e,lne;
    real kappa,i,a1,a3,a11,a13,a33,a111,a113,a133,a333,aconst;
    real x4,y4,z4,x6,y6,z6;

    aa=a*a;
    bb=b*b;
    xx=x*x;
    yy=y*y;
    rr=xx+yy;
    side=xx/aa+yy/bb-1.0;
    gi=1.0;
    if (side>0.0) {
        kappa=(rr-aa-bb+sqrt(sqr(rr-aa-bb) + 4.0*aa*bb*side)) * 0.5;
        aa += kappa;
        bb += kappa;
        gi=a*b*b/(sqrt(aa)*bb);
    }
    em=bb/aa;
    ee=1.0-em;
    e=sqrt(ee);
    lne=log((1.0+e)/(1.0-e));
                                /* index symbols */
    a3=gi*(1.0-em*lne/e*0.5)/ee;
    a1=2.0*(gi-a3);
    i=aa*a1+2.0*bb*a3;
    
    a33 = gi*(4*ee/em-6+3*em*lne/e)/(8*aa*sqr(ee));
    a13 = 2*(gi/bb-2*a33);
    a11 = 2*(gi/aa-a13)/3;
    
    a333 = gi*(30-20*ee/em+16*sqr(ee)/sqr(em)-15*em*lne/e) /
                (48*ee*sqr(ee)*sqr(aa));
    a133 = 2*(gi/sqr(bb)-3*a333);
    a113 = 2*(gi/(aa*bb)-2*a133)/3;
    a111 = 2*(gi/sqr(aa)-a113)*0.2;
                                /* other useful constants */
    x6=a111*xx + 3*a113*yy;
    y6=3*a133*xx + a333*yy;
    z6=3*a133*xx + 3*a333*yy;
    
    x4=a11*xx+a13*yy;
    y4=a13*xx+a33*yy;
    z4=y4;
                                /* and the potential, derivatives etc. */
    *ept = i-3*(a1*xx+a3*(yy))+3*(xx*x4+yy*y4)
         -(x6*sqr(xx)+y6*sqr(yy));

    *ax = -6*x*a1 + 12*x*x4 - 2*x*(2*xx*x6)
         -2*x*(sqr(xx)*a111+3*sqr(yy)*a133);

    *ay = -6*y*a3+  12*y*y4 - 2*y*(2*yy*y6)
         -2*y*(3*sqr(xx)*a113+sqr(yy)*a333);

    *axx = -6*a1 + 12*x4 + 24*a11*xx - 8*xx*(2*a111*xx+x6)
         -2*(2*x6*xx+a111*xx*xx+3*a133*yy*yy);

    *ayy = -6*a3 + 12*y4 + 24*a33*yy - 8*yy*(2*a333*yy+y6)
	 -2*(2*y6*yy+3*a113*xx*xx+a333*yy*yy);

#if defined(TESTBED)
    *axy = -6*a3 + 12*z4 - 12*a133*xx*yy        /* used for Poisson test */
         -2*(3*a113*sqr(xx)+3*a333*sqr(yy));
#else
    *axy = 8*x*y*(3*a13 - 2*a113*xx - 3*a133*yy);	/* ??? very old bug ??? */
    *axy = 24*x*y*(a13 - a113*xx - a133*yy);		/* the correct one */
#endif

    aconst = elipm*Grav_Const*35/(32*a*b*b);
    *ept *= aconst;            /* note ept > 0 for henyey */
    *ax *= aconst;              /* and forces < 0 */
    *ay *= aconst;
    *axx *= aconst;
    *axy *= aconst;
    *ayy *= aconst;

    if (first) {
      printf("# bar6: M,a,b=%g %g %g   aconst=%g\n", elipm, a, b, aconst);
      first = 0;
    }

    
}



void barx_ (real *gm,real *a,real *b,real *x,real *y,
	    real *ept,real *ax,real *ay,real *axx,real *axy,real *ayy)  /* called by FORTRAN as: CALL BARX() */
{
  bar6(*gm,*a,*b,*x,*y,ept,ax,ay,axx,axy,ayy);	/* Call the C routine */
}

#if defined(TESTBED)

main(ac,av)                         /* main to test Possoin */
int ac;
char *av[];
{
  real x,y;
  double atof();
  real ax,ay,ept,axx,ayy,azz;
  real r,norden,pden,a,b;

  x = atof(av[1]);
  y = atof(av[2]);
  a = 1.0;
  b = 0.5;
  bar6(1.0,a,b,x,y,&ept,&ax,&ay,&axx,&azz,&ayy);   /* needs azz hidden in axy !! */
  r = 1.0-x*x/a/a-y*y/b/b;
  norden = (r<0 ? 0 : r*r);
  pden = (axx+ayy+azz);
    
  printf("x=%f y=%f: pot=%f den=%f norden=%f norden/pden=%f\n",
	 x,y,ept,pden,norden,norden/pden);
}

#endif
