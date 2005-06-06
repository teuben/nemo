/*
 * teusan85.c: procedures for intializing and calculating the forces and
 *             potential of a model as described in Teuben & Sanders (1985)
 *	       This bar is oriented along the X axis.
 *      >> This is the 2D version for forces, if a correct version which
 *	   returns Z forces too (potlist needs them if densities are
 *	   required) - compiled with -DNEED3D or so
 *
 *      March 88  - finally according to TS85 paper (sort of)
 *      Feb 90 - added time parameter for new potential-interface
 *      June 91 - resurected and tested for Sellwood & Wilkinson
 *                1) Original 3D:       17.76" CPU
 *                2) 3d->2d version:    14.08"
 *                3) save inside var's:  9.03"
 *                   -O compilation:     8.82"
 *	June 91 - compile option to get full ... 3D
 *	March 92 - happy gcc2.0					pjt
 *	Oct 93  - get_pattern
 *	Aug 96  - changed order of pars to improve readability  Dave Shone/pjt
 *                prototyped and #ifdeffed the NEED3D a bit too.
 *	Jun 01  - dprintf level 0->1
 *      Sep 04  - double/float
 */

/*CTEX
 *  {\bf potname=teusan85}
 *
 * This potential is that of a barred galaxy model as 
 * described in Teuben \& Sanders (1985)
 * This bar is oriented along the X axis.
 * This is the 2D version for forces.  This version should give (near)
 * identical results to {\bf bar83} and very simlar to {\bf athan92}.
 */
 
#include <stdinc.h>
#include <potential_float.h>

/* Official model parameters as they go into inipotential(): */

local double omega = 0.0381;      /* pattern speed (km/s/pc) */
local double fm    = 1.334697416; /* mass ratio bar/disk */
local double fx    = 6.0;         /* a/rmx - usually not varied ... */
local double ca    = 0.2;         /* bar axial ratio */

/* Handy model parameters to keep around: */

local double    M_h, M_b, M_c, A_b, B_b, A_c, A_h;
local double    g_c, g_h;
local double    Grav_Const;
local double    aconst, A2_h, A2_c, A2_b, B2_b;

#if defined(NEED3D)
local void prol6(double,double,double,double,double,double,double,
                 double*,double*,double*,double*);
#else
local void prol6_nr_2d (double,double,
                        double*, double*, double*);
#endif            
/*-------------------------------------------------------------------------
 * INI_POTENTIAL: initializes the potential.
 *-------------------------------------------------------------------------
 */

void inipotential (int *npar, double *par, string name)
{
    double vmx, rmx;

    if (*npar>0) omega = par[0];	/* note; get_pattern */
    if (*npar>1) fm=par[1];
    if (*npar>2) fx=par[2];
    if (*npar>3) ca=par[3];

    Grav_Const = 1.0/232.59;            /* Units: km/s, pc,  */

    A_b = 9000.0;                       /* major axis of bar  */
    vmx = 270.1097818;                  /* max rotcur scaling factor */

    B_b = A_b * ca;                     /* minor axis of bar */
    rmx = A_b / fx;                     /* unperturbed max rotcur in core */
    A_c = rmx/sqrt(2.0);                /* length scale in core */
    A2_c = sqr(A_c);
    g_c = -sqr(vmx)*rmx*pow(1.5,1.5);   /* GM of core */
    M_c = -g_c/Grav_Const;              /* Mass of core */

    A_h = 14140.0;                      /* halo */
    A2_h = sqr(A_h);
    g_h = -9.9e8;
    M_h = -g_h/Grav_Const;

    g_c /= (1+fm);                      /* redistribute mass to bar/core */
    M_b = M_c*fm/(1+fm);
    M_c /= (1+fm);

#if defined(NEED3D)
    dprintf (1,"INI_POTENTIAL TeuSan85 [Full 3D] \n");
#else
    dprintf (1,"INI_POTENTIAL TeuSan85 [optimzed 2D] \n");
#endif
    dprintf (1,"Parameters : \nPattern Speed = %f \n",omega);
    dprintf (1,"fm=%f   fx=%f  ca=%f\n\n",fm,fx,ca);
    dprintf (2,"Derived parameters:\n");
    dprintf (2,"Halo: M_h=%g A_h=%g\n",M_h, A_h);
    dprintf (2,"Core: M_c=%g A_c=%g\n",M_c, A_c);
    dprintf (2,"Bar:  M_b=%g A_b=%g B_b=%g\n",M_b, A_b, B_b);

/*  save some other variables for this optimized version */

    aconst = Grav_Const*M_b*35/(32*A_b*B_b*B_b);
    A2_b = sqr(A_b);
    B2_b = sqr(B_b);

    par[0] = omega;
}

/*------------------------------------------------------------------------------
 *  POTENTIAL: the worker routine: gets forces and pot at (x,y) (z=dummy in 2D)
 *------------------------------------------------------------------------------
 */
void potential_double (int *ndim,double *pos,double *acc,double *pot,double *time)
{
    double rr;
    double dn_h, dn_c, ftr;

#if defined(NEED3D)
    prol6(M_b, A_b, B_b, B_b, pos[0],  pos[1],  pos[2],
                             &acc[0], &acc[1], &acc[2], pot);
    rr = pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
#else
    prol6_nr_2d(pos[0], pos[1], &acc[0], &acc[1], pot);
    dprintf(3,"NPROL6_2d: fac=%g accx,accy,pot=%g %g %g\n",           
        aconst, acc[0], acc[1], *pot);
    acc[2] = 0.0;                            /* 2D : set z-force 0 */
    rr = pos[0]*pos[0] + pos[1]*pos[1];
#endif

    dn_h = 1.0/sqrt(A2_h+rr);
    dn_c = 1.0/sqrt(A2_c+rr);
    *pot += g_c*dn_c + g_h*dn_h;

    if (rr==0.0)                    /* no forces at r=0 from spheres */
        return;                     /* can return now */

    ftr = g_c * dn_c*dn_c*dn_c + g_h * dn_h*dn_h*dn_h;  /* radial force */

    acc[0] += ftr*pos[0];
    acc[1] += ftr*pos[1];
#if defined(NEED3D)
    acc[2] += ftr*pos[2];
#endif
}

/*
 *    PROL6:   the k=2 Ferrers prolate bar (2D non-re-usable!!)
 *                        x,y   position of point of interest (I)
 *                      ax,ay   forces at point               (O)
 *                        ept   potential at point            (O)
 */

#if !defined(NEED3D)

local int first=1;
local double i_i, a1_i, a3_i, a11_i, a13_i, a33_i, 
              a111_i, a113_i, a133_i, a333_i;

local void
prol6_nr_2d (double x,double y,                         /* input */
             double *ax, double *ay, double *ept)       /* output */
{
    double aa,cc,xx,yy,rr,side,gi, x4,y4,x6,y6; /* temp variables */
    double i,a1,a3,a11,a13,a33,a111,a113,a133,a333; /* index symbols */

    aa=A2_b;
    cc=B2_b;
    xx=x*x;
    yy=y*y;
    rr=xx+yy;
    side=xx/aa+yy/cc-1.0;
    gi=1.0;
    if (side>0.0) {
        double kappa;               /* temp variable if outside bar */

        kappa=(rr-aa-cc+sqrt(sqr(rr-aa-cc) + 4.0*aa*cc*side)) * 0.5;
        aa += kappa;
        cc += kappa;
        gi=A_b*B2_b/(sqrt(aa)*cc);
    } 
    if (side>0.0 || first) {        /* make new index symbols */
       double em, ee, e, lne;           /* temporary variables */

       em=cc/aa;
       ee=1.0-em;
       e=sqrt(ee);
       lne=log((1.0+e)/(1.0-e));

       a3=gi*(1.0-em*lne/e*0.5)/ee;
       a1=2.0*(gi-a3);
       i=aa*a1+2.0*cc*a3;
    
       a33 = gi*(4*ee/em-6+3*em*lne/e)/(8*aa*ee*ee);
       a13 = 2*(gi/cc-2*a33);
       a11 = 2*(gi/aa-a13)/3;
    
       a333 = gi*(30-20*ee/em+16*ee*ee/(em*em)-15*em*lne/e) /
                (48*ee*ee*ee*aa*aa);
       a133 = 2*(gi/(cc*cc)-3*a333);
       a113 = 2*(gi/(aa*cc)-2*a133)/3;
       a111 = 2*(gi/(aa*aa)-a113)*0.2;
    } else {                        /* retrieve inside coeffs */
       i = i_i;   
       a1 = a1_i;  a3 = a3_i;
       a11 = a11_i; a13 = a13_i; a33 = a33_i;
       a111 = a111_i; a113 = a113_i; a133 = a133_i; a333 = a333_i;
    }
    if (first && side<=0.0) {           /* save inside coefs */
       dprintf(1,"Saving inside coeffs\n");
       first = 0;
       i_i = i;   
       a1_i = a1;  a3_i = a3;
       a11_i = a11; a13_i = a13; a33_i = a33;
       a111_i = a111; a113_i = a113; a133_i = a133; a333_i = a333;
    }
    
    x4=a11*xx+a13*yy;
    y4=a13*xx+a33*yy;

    x6=a111*xx + 3*a113*yy;
    y6=3*a133*xx + a333*yy;

    *ept = i-3*(a1*xx+a3*yy)+3*(xx*x4+yy*y4)
         -(x6*xx*xx+y6*yy*yy);
    *ax = -6*x*a1 + 12*x*x4 - 2*x*(2*xx*x6)
         -2*x*(xx*xx*a111+3*yy*yy*a133);
    *ay = -6*y*a3 + 12*y*y4 - 2*y*(2*yy*y6)
         -2*y*(3*xx*xx*a113+yy*yy*a333);

    *ept *= -aconst;        /* rescale to physical units */
    *ax *= aconst;
    *ay *= aconst;
}
#else

/*
 *    PROL6:    the k=2 Ferrers prolate bar - more or less original code
 *              input:  a,b,c   axes of bar (b is dummy since prolate)
 *                      x,y,z   position of point of interest
 *                   ax,ay,az   forces at point
 *                        ept   potential at point
 *
 *	See also:
 *	S. Chandrasekhar (1969) 'Ellipsoidal Figures of Equilibrium' par.21
 */

local void prol6 (double elipm,double a,double b,double c,      /* input */
                               double x,double y,double z,
                  double *ax,double *ay,double *az,double *ept) /* output */
{
    double aa,cc,xx,yy,zz,rr,gi,side,em,ee,e,lne;
    double kappa,i,a1,a3,a11,a13,a33,a111,a113,a133,a333,aconst;
    double x4,y4,z4,x6,y6,z6;

    aa=a*a;
    cc=c*c;
    xx=x*x;
    yy=y*y;
    zz=z*z;
    rr=xx+yy+zz;
    side=xx/aa+(yy+zz)/cc-1.0;
    gi=1.0;
    if (side>0.0) {
        kappa=(rr-aa-cc+sqrt(sqr(rr-aa-cc) + 4.0*aa*cc*side)) * 0.5;
        aa += kappa;
        cc += kappa;
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
    
    a33 = gi*(4*ee/em-6+3*em*lne/e)/(8*aa*sqr(ee));
    a13 = 2*(gi/cc-2*a33);
    a11 = 2*(gi/aa-a13)/3;
    
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

    aconst = elipm*Grav_Const*35/(32*a*c*c);
    *ept *= -aconst;
    *ax *= aconst;
    *ay *= aconst;
    *az *= aconst;
}
#endif

#if defined(TESTBED)
main()
{
  double dum[1], pos[3], acc[3], pot, ptime;
  int ndim=3, npar=0;

  inipotential(&npar,dum,"");
  potential(&ndim,pos,acc,&pot,&ptime);
}

double sqr(x)
double x;
{ return (x*x); }
#endif
