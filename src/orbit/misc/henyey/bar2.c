/* BAR2: prolate ferrers k=0 (homegeneous) for HENYEY - derived from BAR6 */
extern double sqrt(double), log(double);         /* externals used in here */

double sqr(double x)
{ 
  return x*x; 
}

double Grav_Const = 1.0;        /* this makes you have to enter GM instead */

int PotBit = (1<<0);
int FxBit  = (1<<1);
int FyBit  = (1<<2);
int FzBit  = (1<<3);
int FxxBit = (1<<4);
int FyyBit = (1<<5);
int FzzBit = (1<<6);
int FxyBit = (1<<7);
int FxzBit = (1<<8);
int FyzBit = (1<<9);
#if 0
int F0Bit = PotBit;
int F1Bit = (FxBit|FyBit|FzBit);
int F2Bit = (FxxBit|FyyBit|FzzBit|FxyBit|FxzBit|FyzBit);
#endif

/*
 *    BAR2:    the k=0 Ferrers prolate bar
 *             this one is for 2D though (25-sep-91)
 *             used to return second derivatives... for HENYEY...
 *             input: elipm   mass of bar
 *                      a,b   axes of bar  (a>b)
 *                      x,y   position of point of interest
 *                    ax,ay   forces [-d(potential)/dr]
 *                      ept   potential
 *              axx,ayy,azz   second derivatives of potential - for testing
 */
bar2 (elipm,a,b,x,y,ept,ax,ay,axx,axy,ayy)
double elipm,a,b,x,y;                          /* input */
double *ax,*ay,*ept,*axx,*axy,*ayy;            /* output */
{
    double aa,bb,xx,yy,zz,rr,gi,side,em,ee,e,lne;
    double kappa,i,a1,a3,a11,a13,a33,a111,a113,a133,a333,aconst;
    double x4,y4,z4,x6,y6,z6;

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
                                /* and the potential, derivatives etc. */
    *ept = i-(a1*xx+a3*(yy));

    *ax = -6*x*a1;

    *ay = -6*y*a3;

    *axx = -6*a1;

    *ayy = -6*a3;

#if defined(TESTBED)
    *axy = -6*a3;
#else
    *axy = 0.0; 	/* can't be good */
#endif

    aconst = elipm*Grav_Const*35/(32*a*b*b);
    *ept *= -aconst;            /* note ept > 0 for henyey */
    *ax *= aconst;              /* and forces < 0 */
    *ay *= aconst;
    *axx *= aconst;
    *axy *= aconst;
    *ayy *= aconst;
}



barx_ (gm,a,b,x,y,ept,ax,ay,axx,axy,ayy)  /* called by FORTRAN as: CALL BARX() */
double *gm,*a,*b,*x,*y;                   /* input */
double *ax,*ay,*ept,*axx,*axy,*ayy;       /* output */
{
  bar2(*gm,*a,*b,*x,*y,ept,ax,ay,axx,axy,ayy);	/* Call the C routine */
}

#if defined(TESTBED)

main(ac,av)                         /* main to test Possoin */
int ac;
char *av[];
{
  double x,y;
  double atof();
  double ax,ay,ept,axx,ayy,azz;
  double r,norden,pden,a,b;

  x = atof(av[1]);
  y = atof(av[2]);
  a = 1.0;
  b = 0.5;
  bar2(1.0,a,b,x,y,&ept,&ax,&ay,&axx,&azz,&ayy);   /* needs azz hidden in axy !! */
  r = 1.0-x*x/a/a-y*y/b/b;
  norden = (r<0 ? 0 : r*r);
  pden = (axx+ayy+azz);
    
  printf("x=%f y=%f: den=%f norden=%f norden/pden=%f\n",
			x,y,pden,norden,norden/pden);
}

#endif

