/*
 * bar83.c: procedures for intializing and calculating the forces and
 *	       potential of a model as described in Teuben et al. (1983)
 *	Note:  only valid in the z=0 plane!
 *
 *	July 1987 - Peter Teuben @ IAS, Princeton, NJ
 *	March 88  - finally according to TS85 paper (sort of)
 *	Feb 90 - added time parameter
 */

/*CTEX
 *  {\bf potname=bar83
 *	 potpars={\it $\Omega,f_m,f_x,{c\over a}$}}
 *
 *  Barred potential as described by Teuben and Sanders (1983),
 *  see also {\bf teusan83}.
 *
 *	Note:  the potential only valid in the z=0 plane!
 */
#include <stdinc.h>
#include <potential_float.h>

		/* composite model parameters */

static double omega = 0.0;		/* pattern speed */	     
static double fm    = 1.334697416;	/* mass ratio bar/disk */
static double fx    = 8.485281374;	/* length scale ratio bar/disk */
static double ca    = 0.2;		/* bar axial ratio */


static double    M_core, M_h, M_b, M_c, A_b, B_b, A_c, A_h;  /* handy variables */
static double    Grav_Const;

static void prol6 (double elipm,
		   double a,double b,double c, 
		   double x,double y,double z,
		   double *ax,double *ay,double *az,
		   double *ept);


void inipotential (int *npar, double *par, char *name)
{
    if (*npar>0)
        omega = par[0];
    if (*npar>1)
        fm=par[1];
    if (*npar>2)
    	fx=par[2];
    if (*npar>3)
        ca=par[3];

    M_core = 4.6724720e10;		/* mass of 'core' (bar and bulge) */
    M_h = 4.93 * M_core;		/* mass of halo */
    M_b = M_core / (1.0 + fm)*fm;	/* eventual mass of bar */
    M_c = M_core / (1.0 + fm);		/* eventual mass of bulge/disk */
    A_b = 9000.0;			/* length of bar along major axis */
    B_b = ca * A_b;			/*       "             minor  "   */
    A_c = A_b / fx;			/* length scale of bulge/disk */
    A_h = 1.57 * A_b;			/* length scale of halo */
    Grav_Const = 1.0/232.59;    
    dprintf (1,"INI_POTENTIAL Bar83  name=%s\n",name);
    dprintf (1,"Parameters : \nPattern Speed = %f \n",omega);
    dprintf (1,"fm=%f   fx=%f  ca=%f\n\n",fm,fx,ca);
}
    
void potential_double (int *ndim,double *pos,double *acc,double *pot,double *time)
{
	double rr, r;
	double q_h, q_c, pot_h, pot_c, acc_h, acc_c;

	prol6 (M_b,A_b,B_b,B_b,pos[0],pos[1],pos[2],
					&acc[0],&acc[1],&acc[2],pot);

	rr = sqr(pos[0]) + sqr(pos[1]) + sqr(pos[2]);
        r = sqrt(rr);

	q_h = rr + sqr(A_h);
	pot_h = M_h/sqrt(q_h);
	*pot -= pot_h;			/* add halo potential */

	q_c = rr + sqr(A_c);
	pot_c = M_c/sqrt(q_c);
	*pot -= pot_c;			/* add core potential */

	*pot *= Grav_Const;		/* rescale */

	if (r==0.0)			/* no forces at r=0 */
	  return;			

	acc_h = pot_h/(q_h*r);
	acc[0] -= acc_h * pos[0];
	acc[1] -= acc_h * pos[1];
	acc[2] -= acc_h * pos[2];

	acc_c = pot_c/(q_c*r);
	acc[0] -= acc_c * pos[0];
	acc[1] -= acc_c * pos[1];
	acc[2] -= acc_c * pos[2];

	acc[0] *= Grav_Const;			/* rescale */
	acc[1] *= Grav_Const;
	acc[2] *= Grav_Const;


}

/*
 *    PROL6:	the funny bar
 *		input:	a,b,c	axes of bar (b is dummy)
 *			x,y,z	position of point of interest
 *		     ax,ay,az	forces at point
 *		          ept   potential at point
 */

static void prol6 (double elipm,
		   double a,double b,double c, 
		   double x,double y,double z,
		   double *ax,double *ay,double *az,
		   double *ept)
     
{
    double aa,cc,xx,yy,zz,rr,gi,side,em,ee,e,lne;
    double kappa,i,a1,a3,a11,a13,a33,a111,a113,a133,a333,aconst;
    double x4,y4,z4,x6,y6,z6;
    
    aa=sqr(a);
    cc=sqr(c);
    xx=sqr(x);
    yy=sqr(y);
    zz=sqr(z);
    rr=xx+yy+zz;
    side=xx/aa+(yy+zz)/cc-1.0;
    gi=1.0;
    if (side>0.0) {
        kappa=(rr-aa-cc+sqrt(sqr(rr-aa-cc) + 4.0*aa*cc*side)) * 0.5;
        aa=aa+kappa;
        cc=cc+kappa;
        gi=a*c*c/(sqrt(aa)*cc);
    }
    em=cc/aa;
    ee=1.0-em;
    e=sqrt(ee);
    lne=log((1.0+e)/(1.0-e));
                                /* index symbols */
    a3=gi*(1.0-em*lne/e*0.5)/ee;
    a1=2.0*(gi-a3);
    i=aa*a1+2.0*cc*a3;
    
    a33=1;
    a13=1;
    a11=1;
	a33 = gi*(4*ee/em-6+3*em*lne/e)/(8*aa*sqr(ee));
	a13 = 2*(gi/cc-2*a33);
	a11 = 2*(gi/aa-a13)/3;
    
    a333=1;
    a333=1;



	a333 = gi*(30-20*ee/em+16*sqr(ee)/sqr(em)-15*em*lne/e) /
		(48*ee*sqr(ee)*sqr(aa));
	a133 = 2*(gi/sqr(cc)-3*a333);
	a113 = 2*(gi/(aa*cc)-2*a133)/3;
	a111 = 2*(gi/sqr(aa)-a113)*0.2;
    
    x6=a111*xx + 3*a113*(yy+zz);
    y6=3*a133*xx + a333*(yy+3*zz);
    z6=3*a133*xx + a333*(3*yy+zz);
    



    x4=a11*xx+a13*(yy+zz);
    y4=a13*xx+a33*(yy+zz);
    z4=y4;    

    *ept = i-3*(a1*xx+a3*(yy+zz))+3*(xx*x4+yy*y4+zz*z4)
         -(x6*sqr(xx)+y6*sqr(yy)+z6*sqr(zz)+6*a133*xx*yy*zz);
    *ax = -6*x*a1 + 12*x*x4 - 2*x*(2*xx*x6+6*a133*yy*zz)
         -2*x*(sqr(xx)*a111+3*sqr(yy)*a133+3*sqr(zz)*a133);

    *ay = -6*y*a3+  12*y*y4 - 2*y*(2*yy*y6+6*a133*xx*zz)
	 -2*y*(3*sqr(xx)*a113+sqr(yy)*a333+3*sqr(zz)*a333);
    

    *az = -6*z*a3 + 12*z*z4 - 2*z*(2*zz*z6+6*a133*xx*yy)
	-2*z*(3*sqr(xx)*a113+3*sqr(yy)*a333+sqr(zz)*a333);

    aconst=elipm*35/(32*a*c*c);		/*  G=1 */
    *ept *= -aconst;
    *ax *= aconst;
    *ay *= aconst;
    *az *= aconst;
}


#ifdef TESTBED

#include <stdio.h>

main()
{

	double x,y,z,ax,ay,az,ept;

	printf ("Enter x y z: ");
	scanf("%lf %lf %lf",&x,&y,&z);
	printf (" x y z = %f %f %f\n",x,y,z);

	prol6(1.0,1.0,0.2,0.2,x,y,z,&ax,&ay,&az,&ept);

	printf ("ax ay az = %f %f %f\n",ax,ay,az);
	printf ("ept = %f\n",ept);
}

#endif
