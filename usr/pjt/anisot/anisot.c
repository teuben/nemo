/*
 * Calculate anisotropic distribution function tables for various models
 *
 *	Based on a fortran program by S. Casertano
 *	converted to C with the help of FortrixC (what a pain that was)
 *
 *	Interfaces to IMSL cubic spline and RK integrator have been defined
 *
 *	xx-xxx-8x	V0.0 original program in Fortran  Stefano Casertano
 *	 4-jun-88	V1.0 Translated to C for Nemo          Peter Teuben
 *      28-oct-88       V1.1 in= keyword added + Plummer bug removed    PJT
 *	 9-sep-90	V1.2 the new C-type imsl emulator		PJT
 *	 7-may-92	V1.3 imsl in library - more to be fixed
 *				 fixed up with local odeint/rkqc
 *				 drange -> nemoinpd
 *         nov-93       adding some CTEX comments
 *         jun-97       prototypes to make it compile for NEMO 2.4
 *         nov-01       minor cleanup to make it compile for NEMO 3.x
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#include <math.h>

string defv[] = {
   "out=???\n          name of table output (binary)",
   "model=plummer\n    model (plummer, king, devauc, jaffe, hernquist)",
   "rmax=100\n         Plummer/King : model cuttoff radius",
   "b=0\n              Plummer/King : anisotropy parameter b=1/ra",
   "nrad=512\n         number of radii in integration",
   "stride=1\n         stride from integration table to output table[1=full]",
   "ngauss=48\n        interval subsection for gauss-legendre integration",
   "sigma=0.53192304\n Plummer : vel. disp.?",
   "r0=0.58904862\n    Plummer : scale length",
   "w0=1\n             King : dimensionless central potential",
   "emtot=-1\n         King :",
   "rc=1\n             King : core radius",
   "in=\n              optional (ascii) data file with r,rho stored",
   "radcol=1\n         (in=) column # where to get radii",
   "denscol=2\n        (in=) column # where to get densities",
   "masscol=0\n        (in=) column # where to get cum. masses (optional)",
   "VERSION=1.3b\n     7-nov-01 PJT",
   NULL,
};

string usage = "compute anisotropic distribution function tables";

#define DYNAM	1		/* static gave too many local variables */
#if !defined(DYNAM)
#define NPOINT	100
#endif

int nrad;			/* size of tables */

#define	NGAUTO	2000

int     kmax, kount;                    /* in case a history record of */
double  dxsav, xp[100], yp[10][100];    /* integration is needed */
                                        /* used in NumRec functions */
static double gravconst = 1.0;



nemo_main()
{
#if defined(DYNAM)      
        double *poten, *radius, *distr, *dens, *dens1, *emint, *distr2, *value;
        double *thdis, *aux1, *aux2, *aux3, *cdens1, *cvalue, *qplus;
        int *index;
#else
        double poten[NPOINT], radius[NPOINT], distr[NPOINT], 
                dens[NPOINT], dens1[NPOINT], emint[NPOINT], 
                distr2[NPOINT], value[NPOINT], thdis[NPOINT];
        double aux1[NPOINT], aux2[NPOINT], aux3[NPOINT], 
                cdens1[3][NPOINT], cvalue[3][NPOINT], qplus[NPOINT];
        int index[NPOINT];
#endif
	double abscis[NGAUTO], weight[NGAUTO];
	double b, width, ra, fmin;
	int idens, idens1, ipoten, ivalue, iaux1, iaux2, iaux3;
	int njump, ngauss, i, ier, idistr, j, npri, nsec;
	int ndistr;
	int gauleg(), icsccu(), dcsevu(), plummer(), king(), devauc(), jaffe();
	stream outstr, fstr;
        string fname, model;

        model = getparam("model");          /* model */
        b = getdparam("b");
        nrad = getiparam("nrad");           /* integration constants */
        njump = getiparam("stride");        
        ngauss = getiparam("ngauss");
        fname = getparam("in");             /* optional filename */
        /* REST OF PARAMETERS ARE ASKED IN THEIR RESPECTIVE ROUTINES */

#if !defined (DYNAM)
	nrad = NPOINT;
#endif	

	idens = nrad;
	idens1 = nrad;
	ipoten = nrad;
	ivalue = nrad;
	iaux1 = nrad;
	iaux2 = nrad;
	iaux3 = nrad;

#if defined(DYNAM)
	poten  = (double *) allocate(nrad*sizeof(double));
	radius = (double *) allocate(nrad*sizeof(double));
	distr  = (double *) allocate(nrad*sizeof(double));
	dens   = (double *) allocate(nrad*sizeof(double));
	dens1  = (double *) allocate(nrad*sizeof(double));
	emint  = (double *) allocate(nrad*sizeof(double));
	distr2 = (double *) allocate(nrad*sizeof(double));
	value  = (double *) allocate(nrad*sizeof(double));
	thdis  = (double *) allocate(nrad*sizeof(double));
	aux1   = (double *) allocate(nrad*sizeof(double));
	aux2   = (double *) allocate(nrad*sizeof(double));
	aux3   = (double *) allocate(nrad*sizeof(double));
	aux1   = (double *) allocate(nrad*sizeof(double));
	cdens1 = (double *) allocate(nrad*sizeof(double) * 3);
	cvalue = (double *) allocate(nrad*sizeof(double) * 3);
	qplus  = (double *) allocate(nrad*sizeof(double) * 3);
	index  = (int *)    allocate(nrad*sizeof(int));
#endif
/*CTEX
 * 
 *         This program is calculating the distribution function
 *      for a given spherical density distribution.
 *      The density distribution will be
 *      given numerically, in the form of a table.  The potential will al
 *      have to be computed numerically.  The models are supposed to be 
 *      regular in the center, and finite.
 * 
 *         Once density and potential are known, the distribution
 *      function is be computed with the aid of the formula
 * $$
 *            f(E) = {\sqrt{2} \over 4\pi} {d \over dE} \int_E^0
 *         {d \rho_1 \over d U} { d U \over \sqrt{U-E} } \> ,
 * $$
 *      $\rho_1$  being the `corrected' (in Merritt's sense) density.
 * 
 *         It may be possible to integrate by parts.  In this case, one
 *      can avoid either numerical differentiation or integration over a
 *      singularity, possibly both.  We will indeed use the alternative me
 *      of integrating once by parts to calculate  DISTR2, as a check.
 *      The first tests seem to indicate that DISTR2 is LESS accurate than
 * 
 *         The vectors RADIUS, POTEN and DENS will contain the potential a
 *      density as a function of radius.  They will be filled in a separat
 *      subroutine, which is also expected to fill in  EMINT,  the mass in
 *      to a given radius, and  THDIS,  the theoretical distribution
 *      function (when known at least for the isotropic problem).
 * 
 */

    gauleg(abscis, weight, ngauss);

/*CTEX
 * 
 *         The subroutine  $GAULEG$  will calculate abscissae $X_i$ and
 *      weights $W_i$ for the Gauss-Legendre integration with $N = 2*NGAUSS$
 *      (from Press et al, {\it Numerical Recipes}, p. 110).
 *      They have been tested - for $NGAUSS = 48$ - against the values
 *      in A\&S, Table 25.4, pp. 916-919.  $GAULEG$ returns the values:
 * $$
 *      ABSCIS(i) = {X_i}^2
 * $$
 *	and
 * $$
 *      WEIGHT(i) = 2 * W_i
 * $$
 *      which are appropriate for integrals of the form:
 * $$
 *      I = \int_a^b {f(t) \, dt \over \sqrt(t-a)} \, ;
 * $$
 *      in fact we have
 * $$
 *      I \approx \sqrt{b-a} \cdot \sum_{i=1)^n WEIGHT(i) \, f(t_i) \, ,}
 * $$
 *      with  $ t_i = a + ABSCIS(i) * (b-a)$
 * 
 *      (derived from AS 25.4.36, with a simple change of variable).
 * 
 * -----------------------------------------------------------------------
 * 
 */
    if (*fname) {
        fstr = stropen(fname,"r");
        read_file (fstr, radius, dens, poten, emint, thdis, nrad, b);
        strclose(fstr);
    } else {
        if(scanopt(model,"king"))
            king        (radius, dens, poten, emint, thdis, nrad, b);
        else if (scanopt(model,"plummer"))
            plummer     (radius, dens, poten, emint, thdis, nrad, b);
        else if (scanopt(model,"devauc"))
            devauc      (radius, dens, poten, emint, thdis, nrad, b);
        else if (scanopt(model,"jaffe"))
            jaffe       (radius, dens, poten, emint, thdis, nrad, b);
        else if (scanopt(model,"hernquist"))
            hernquist   (radius, dens, poten, emint, thdis, nrad, b);
        else
            error("Unknown model [king,plummer,devauc,jaffe,hernquist]\n");
    }
    ra = (b>0.0 ? 1.0/b : -1.0);		/* anisotropy radius */
    printf (" model complete\n");

/*CTEX
 *         Define the `corrected' density  $DENS1$  by
 * $$ 
 *                DENS1 = DENS * (1 + (B*R)**2)
 * $$ 
 *      $B$  is the anisotropy parameter $= 1/RA$ (in Merritt's notation).
 *      $B=0$ for isotropic models
 * 
 */

    for (i=0; i<nrad; i++)
        dens1[i] = dens[i] * (1.0 + sqr(b * radius[i]));

/*
 * 
 *         Now prepare the interpolation subroutine for the density as
 *      a function of the potential
 * 
 */
    icsccu(poten, dens1, nrad, cdens1, idens1, ier);
/*CTEX
 * 
 *         In order to calculate the distribution function we will need
 *      first the value of the integral
 * $$ 
 *      \int_Q^0 {dU \over \sqrt{U-Q)} {d \rho1 \over dU}}
 * $$ 
 *      This we will calculate for values of  Q  coinciding
 *      with tabulated values of the potential.  The integral will
 *      be performed by Gauss's method (see beginning).
 * 
 */
    for (i=0, idistr=0; i<nrad; i += njump, idistr++) {
        index[idistr] = i;
        qplus[idistr] = poten[i];
        width =  - qplus[idistr];
/*
 * 
 *         lower limit of integration and width of the interval
 * 
 */
        value[idistr] = 0.0;
        distr2[idistr] = 0.0;
        if(i == (nrad-1))
            break;			/* DONE !! */
	for (j=0; j<ngauss; j++)
            aux1[j] = qplus[idistr] + abscis[j] * width;
/*
 * 
 *      AUX1  now contains the abscissae of the points where the function
 *      has to be evaluated;  AUX2  will contain the values (of d\rho/dU);
 *      AUX3  those of the second derivative.
 * 
 */
        npri = ngauss;
        nsec = ngauss;
        dcsevu(poten, dens1, nrad, cdens1, idens1, 
                		aux1, aux2, npri, aux3, nsec, ier);
/*
 * 
 *         values available; sum up for the integral
 * 
 */
        for (j=0; j<ngauss; j++) {
            value[idistr] += aux2[j] * weight[j];
            distr2[idistr] += aux3[j] * weight[j];
        }
/*CTEX
 *         Finally, we need to 
 *         correct for the scale and for the factor
 *         $ {\sqrt{2} \over 4 \pi}$.
 * 
 */
        value[idistr] *= sqrt(2*width) / (4 * PI*PI);
        distr2[idistr] *= sqrt(2*width) / (4 * PI*PI);  /* PJT TEST */
/*
 * 
 *         fill in the vectors reporting density, mass and theoretical
 *      distr funct at selected points
 * 
 */
#if 0 
 	printf ("(%d %d): %f %f\n",i,idistr,value[idistr],distr2[idistr]);
#endif
    } /* for (i,idistr) */
    ndistr = idistr+1;		/* FORTRAN AND C differ one here */
    printf (" distribution function computed with ndistr=%d\n",ndistr);

/*CTEX
 * 
 *      Now, the distribution function  $f$  is given by
 * $$
 *      f = {\sqrt 2 \over 4 \pi} {d \over d Q_+} {the previous integral}
 * $$
 *      Again, splines will do the job for the differentiation.
 * 
 */
    icsccu(qplus, value, ndistr, cvalue, ivalue, ier);
    for (idistr=0; idistr < ndistr; idistr++)
        aux1[idistr] = qplus[idistr];

    npri = ndistr;
    nsec = 0;
    dcsevu(qplus, value, ndistr, cvalue, ivalue, 
        			aux1, distr, npri, aux2, nsec, ier);
	
    outstr = stropen(getparam("out"),"w");
	        /* BUG: THIS IS NOT CORRECT WHEN NJUMP > 1 */
    put_set(outstr,  "OsipkovMerrittModel");
    put_data(outstr, "AnisoRadius", RealType, &ra, 0);
    put_data(outstr, "Ntab", IntType, &nrad, 0);
    put_data(outstr, "Radius", RealType, radius, nrad, 0);
    put_data(outstr, "Density", RealType, dens, nrad, 0);
    put_data(outstr, "Mass", RealType, emint, nrad, 0);
    put_data(outstr, "Potential", RealType, qplus, nrad, 0);
    put_data(outstr, "DistFunc", RealType, distr2, nrad, 0);    /* distr2 better ? */
    put_tes(outstr, "OsipkovMerrittModel");
    strclose(outstr);

    dprintf (2,"First table:\n");
    for (idistr=0; idistr < ndistr; idistr++) {
        i = index[idistr];
        dprintf (2,"%d: %f %f %f %f %f %f %f %f\n",
		i,radius[i], dens[i], emint[i], qplus[idistr], value[idistr], 
		distr[idistr], distr2[idistr], thdis[i]);
    }
    nsec = 0;       /* count negative DF's */
    fmin = 0.0;
    for (i=0; i<ndistr; i++)
        if (distr2[i] < fmin) {
            fmin = distr2[i];
            nsec++;
        }
    if (nsec>0)
        error("DF < 0 in %d out of %d points; DF_min\n",nsec, ndistr, fmin);
}

/*
 * GAULEG: Finds the abscissae and weights for the Gauss-Legendre
 *         integration scheme.  To be used in DISTR.
 * 
 */

gauleg( abscis, weight, ngauss)
double abscis[];
double weight[];
int ngauss;
{
    double x[NGAUTO], w[NGAUTO], z, p1, p2, p3, pp, z1;
    int n, m, i, j;

    dprintf(0,"gauleg ... creating abscis\n");
    n = 2 * ngauss;

    m = (n + 1) / 2;

    for (i=0; i<m; i++) {
        z =cos(PI * (i + 0.75) / (n + 0.50));
        do {
            p1 = 1.0;
            p2 = 0.0;
	    for (j=0; j<n; j++) {
                p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j+1.0);
            }
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;
        } while (ABS(z-z1) > 5.e-13);		/* WAS: -16 */
        x[i] =  - z;
        x[n - i] = z;		/* THIS DOESN"T SEEM TO BE USED */
        w[i] = 2.0 / ((1.0 - z * z) * pp * pp);
        w[n - i] = w[i];
        abscis[(ngauss - i)] = z * z;
        weight[(ngauss - i)] = 2 * w[i];
    } /* for(i) */
                    /* check accuracy of integrations by this method */
    p1 = p2 = 0;
    for (i=0; i<ngauss; i++) {
        p1 += weight[i];                /* p1 = \int_0^1{dt \over \sqrt(t)} */
        p2 += weight[i] * abscis[i];    /* p2 = \int_0^1{t dt \over \sqrt(t)} */
    }
    dprintf (1,"GAULEG: sum(wi)=%f [exact: 2] sum(wi.ai)=%f [exact: 2/3]\n",
                p1, p2);
}

int plummer(radius, dens, poten, emint, thdis, nrad, b)
double radius[];
double dens[];
double poten[];
double emint[];
double thdis[];
int nrad;
double b;
{
    double rmax, sigma, r0, rstep, rhocen, thdis0, qtrue, qqq;
    int i, idens = nrad;
    static int nexp = 2;
    int calcpot_mass();

    printf("Creating tables for Plummer model\n");

    rmax = getdparam("rmax");       /* cutoff radius */
    sigma = getdparam("sigma");     /* vel.disp? */
    r0 = getdparam("r0");           /* scale length */

    rstep = rmax / pow((double)(nrad - 1), (double)nexp);
    rhocen = 9.0 * sqr(sigma) / (TWO_PI * gravconst * sqr(r0));
    printf ("Creating radius and dens for Plummer\n");
    for (i=0; i<nrad; i++) {
        radius[i] = pow((double)i, (double)nexp) * rstep;
        dens[i] = rhocen / 
                  pow(1.0 + sqr(radius[i]) / sqr(r0), 2.5);
        dprintf (3,"(%d) radius,dens=%f %f\n",i,radius[i],dens[i]);
    }
/*
 * 
 *         calculate gravitational potential (in CALCPOT)
 * 
 */
    calcpot_mass(radius, dens, poten, emint, nrad);

/*
 * 
 *         fill in the theoretical distribution function
 * 
 */
    thdis0 = sqrt(2.0) / (378 * pow(PI, 3.0) * gravconst * sqr(r0) * sigma);
    for (i=0; i<nrad; i++) {
        qtrue =  - 6.0 * sqr(sigma) / 
                   sqrt(1.0 + sqr((radius[i] / r0)));
        qqq =  - qtrue / sqr(sigma);
        thdis[i] = thdis0 * (pow(qqq,3.5) * (1.0 - sqr(b * r0)) + 
       	             pow(qqq,1.5) * 63.0 / 4.0 * sqr(b * r0) );
    }

    printf("Plummer model: sigma r0 rho0=%f %f %f\n", sigma, r0, rhocen);
    printf("               rmax totmas= %f %f\n", rmax, emint[nrad - 1]);
    printf("               b = %f\n",b);
}

int devauc(radius, dens, poten, emint, thdis, nrad, b)
double radius[];
double dens[];
double poten[];
double emint[];
double thdis[];
int nrad;
double b;
{
    error("De Vaucouleurs models not implemented yet");
}

int jaffe(radius, dens, poten, emint, thdis, nrad, b)
double radius[];
double dens[];
double poten[];
double emint[];
double thdis[];
int nrad;
double b;
{
    error("Jaffe models not implemented");
}

int hernquist(radius, dens, poten, emint, thdis, nrad, b)
double radius[];
double dens[];
double poten[];
double emint[];
double thdis[];
int nrad;
double b;
{
    error("Hernquist models not implemented");
}


/*CTEX
 *
 * {\bf CALCPOT\_MASS}: calculate potential (and its minimum) from
 *  density profile
 * 
 *  	input:  radius,dens	(both arrays of length nrad)
 *
 *	output:	poten,emint	
 *
 *         This function will compute the potential and the mass inside
 *      r given the density.  It will assume that the density, as given
 *      in the vector  $DENS$  and spline-interpolated, is exact.
 *      We will define two auxiliary functions,  $EM$  and  $TI$, being the
 *      integrals of  $RHO*R*R$  and of  $RHO*R$:
 * $$
 *         EM(r) = 4 \pi \int_0^r \rho(s) s^2 ds
 * $$
 *	and
 * $$
 *         TI(r) = 4 \pi \int_0^r \rho(s) s ds
 * $$
 *	and
 * $$
 *         U(r) = U(0) + G * (TI(r) - EM(r)/r)
 * $$
 *      (the last is obtained from the definition of the potential, after
 *      exchanging the order of integration)
 * 
 *      The gravitational potential is defined to be zero at the outer
 *      boundary (supposed finite = $RADIUS(NRAD)$), and negative inside
 *      (a physicist's definition).
 * 
 */
 
int calcpot_mass(radius, dens, poten, emint, nrad)
double radius[];
double dens[];
double poten[];
double emint[];
int nrad;
{
#if defined(DYNAM)
    double *cdens;
#else
    double cdens[3][NPOINT];
#endif	
    double em, ti, s, c0, c1, c2, c3, r;
    double a0, a1, a2, a3, a4, b0, b1, b2, b3;
    double b4, b5, rout, emtot, poten0;
    int ier, i;
    int icsccu();
    static int idens = 0;

    if (idens < nrad) {
        if (idens == 0) {   /* probably first time around */
#if defined(DYNAM)
	    cdens = (double *) allocate(nrad*3*sizeof(double));
#else
            idens = nrad;
#endif		
        } else {
#if defined(DYNAM)
	    cdens = (double *) reallocate(cdens, 3*nrad*sizeof(double));
            if (cdens==NULL)
            idens = nrad;
#else
            error("Cannot reallocate\n");
#endif		
        }
    }
    /* set up spline coefficients */
    icsccu(radius, dens, nrad, cdens, idens, ier);

/*
 *         Now we know the coefficients of the spline interpolation
 *      routine.  The density is described as a polynomial inside each
 *      subinterval.  Therefore the integration -- to calculate the
 *      quantities  EM  and  TI  -- is straightforward, implying only
 *      recalculating the relevant coefficients.
 *      In particular, we will define the coefficients  a0-a4 and b0-b5,
 *      of the polynomial representation of  \rho \cdot r and of
 *      \rho \cdot r^2  respectively.
 * 
 */
    em = 0.0;
    ti = 0.0;
    poten[0] = 0.0;
    emint[0] = em;
/*
 *         assumes  RADIUS[0] = 0  -- necessary for the interpolation
 *      routines to work correctly
 *		PJT: This has to be fixed, e.g. models with rho(0) = infinity
 *		     would have some problems here, won't they?
 * 
 */
    for (i=1; i<nrad; i++) {
        s = radius[i] - radius[i-1];
        c0 = dens[i - 1];
#if defined(DYNAM)
        c1 = cdens[i-1];
    	c2 = cdens[nrad + i-1];
	c3 = cdens[nrad*2 + i-1];
#else                
        c1 = cdens[0][i - 1];
        c2 = cdens[1][i - 1];
        c3 = cdens[2][i - 1];
#endif                
        r = radius[i-1];
        a0 = c0 * r;
        a1 = c1 * r + c0;
        a2 = c2 * r + c1;
        a3 = c3 * r + c2;
        a4 = c3;
        b0 = (c0 * r) * r;
        b1 = (c1 * r + 2.0 * c0) * r;
        b2 = (c2 * r + 2.0 * c1) * r + c0;
        b3 = (c3 * r + 2.0 * c2) * r + c1;
        b4 = (2.0 * c3) * r + c2;
        b5 = c3;
        ti +=  FOUR_PI * s * (a0 + s * (a1 / 2.0 + s * 
            	(a2 / 3.0 + s * (a3 / 4.0 + s * a4 / 5.0))));
        em +=  FOUR_PI * s * (b0 + s * (b1 / 2.0 + s * 
            	(b2 / 3.0 + s * (b3 / 4.0 + s * 
            		(b4 / 5.0 + s * b5 / 6.0)))));
        poten[i] = gravconst * (ti - em / radius[i]);
        emint[i] = em;
    } /* for (i=1;) */
    rout = radius[nrad - 1];
    emtot = emint[nrad - 1];
    poten0 =  - poten[nrad - 1];
    printf ("Last potential (for renorm) = %f\n",poten0);
/*
 * 
 *         Redefine the zero of the potential
 * 
 */
    for (i=0; i<nrad; i++) {
        poten[i] = poten[i] + poten0;
        dprintf (3,"calcpot_mass(%d): R,PSI,M= %f %f %f\n",
            i,radius[i],poten[i],emint[i]);
   }
}

/*
 * CALCPOT_DENS: calculate potential (and its minimum) from `
 *               cumulative mass profile
 *  	input:  radius,emint	(both arrays of length nrad)
 *	output:	poten,dens	
 */
 
int calcpot_dens(radius, dens, poten, emint, nrad)
double radius[];
double dens[];
double poten[];
double emint[];
int nrad;
{
    error("calcpot_dens: not implemented yet");
}

double w0;
double rho0;

int king (radius, dens, poten, emint, thdis, nrad, b)
double radius[];
double dens[];
double poten[];
double emint[];
double thdis[];
int nrad;
double b;
{
    double yy[2], wk[9][2], comm[24], emtot, rc, wstep, tol;
    double w, wend, sigma, rhocen, thdis0;
    double func(), rho();
    int i, neqs, ind, nwk, ie;

    w0 = getdparam("w0");           /* dimensionless central potential */
    emtot = getdparam("emtot");     /* ?? */
    rc = getdparam("rc");           /* core radius */
     
    w0 = -ABS(w0);   
    wstep = -(w0) / (nrad - 1);
    rho0 = rho(w0);
    for (i=0; i<nrad; i++)
        poten[i] = w0 + wstep * i;
/*
 * 
 *         first solve for the unscaled model (depending on  W0  only)
 * 
 */
        neqs = 2;
        yy[0] = 0.0;
        yy[1] = 2.0 / 3.0;
        tol = 1.e-6;
        ind = 1;
        nwk = 2;
        w = w0;
        radius[0] = 0.0;
        dens[0] = 1.0;
        emint[0] = 0.0;

        for (i=1; i<nrad; i++) {
                wend = poten[i];
		dverk (&neqs, func, &w, yy, &wend, &tol, &ind, comm,
			&nwk, wk, &ie);
                radius[i] = sqrt(yy[0]);
                dens[i] = rho(w) / rho0;
                emint[i] = 2.0 * yy[0] * sqrt(yy[0]) / yy[1];
        }
/*
 * 
 *      now find the scaling factors (from the core radius, given, and the
 *      total mass) and calculate the theoretical ISOTROPIC distr func
 * 
 */
        sigma = sqrt( gravconst * emtot / rc / emint[((nrad) - (1))] );
        rhocen = (9.0 * sqr(sigma)) / (FOUR_PI * gravconst * sqr(rc));
        thdis0 = 9.0 / (FOUR_PI * gravconst * sqr(rc) * 
        		sigma * FOUR_PI * sqrt(2.0) * rho0);

	printf ("King model w0=%f M rc sigma rhocen= %f %f %f %f b=%f\n",
		w0, emtot, rc, sigma, rhocen, b);

	for (i=0; i<nrad; i++)
        {
                radius[i] = radius[i] * rc;
                emint[i] = emint[i] * 
                	sqr(sigma) * rc / gravconst;
                dens[i] = dens[i] * rhocen;
                thdis[i] = thdis0 * (exp(- poten[i] ) - 1);
                poten[i] = poten[i] * sqr(sigma);
        }
}

/* NOTE that the parameter 'n' has been left out in this C-version */

/*CTEX
 * 
 *         The differential equation to be solved is
 * $$
 *      x" - (3/2) (x')^2 + (9/4) (rho(w)/rho(w0) (x')^3 = 0
 * $$
 *      where             $x = radius^2$
 * 
 *         The two components of the vector  $Y$  are  $X$, $X'$ respectively.
 * 
 */

double func(w, y, yprime)
double w;
double y[];
double yprime[];
{
    double rho();

    yprime[0] = y[1];
    if (y[0] > 1.e-8)
        yprime[1] = 1.50 * y[1] * y[1] * (1.0 - 1.5 * y[1] * 
                	rho(w) / rho0) / y[0];
    if (y[0] <= 1.e-8)
        yprime[1] = 0.40 * (1.0 + pow(-w0, 1.5) / rho0);
/*CTEX
 * 
 *         Despite the apparent singularity in the first form, the
 *      differential equation - with the right initial conditions - has a
 *      regular solution.  The second statement amounts to saying that
 * $$
 *                       X" = - (2/5) rho'(w0)/rho(w0)
 * $$ 
 */
    return(0.0);
}



/*CTEX
 * $$
 *      COEFF  =  \sqrt{\pi} / 2  -- to convert  ERF  into the integral
 * $$
 * 
 *         The function  RHO  calculates the integral
 * $$ 
 *      rho = [\exp(y^2) \int_0^y \exp(-v^2) dv] - y - 2 y^3 / 3 
 *                                                       (y = \sqrt{-w})
 * $$
 * 
 *      In terms of this, the density is defined as
 * $$ 
 *      dens = 4 \sqrt(2) \pi k \sigma^3 \exp(w0)
 * $$ 
 *      where  $k$  is the constants in the definition of the distribution
 *      function, which can be expressed in terms of the other quantities
 *      by using the condition
 * $$ 
 *         dens(w0) = {9 \sigma^2 \over 4 \pi g rc^2},
 * $$ 
 *      so that
 * $$ 
 *         k = {9 \over {(4 \pi)}^2 g rc^2 \sigma \sqrt{2} rho(w0) \exp{w0}}
 * $$ 
 */

double rho(w)
double w;
{
    double Frho, y;
    static double coeff = 0.88622692550;

    Frho = 0.0;
    if (w > 0.0)
       return(Frho);

    y = sqrt(-w);
    Frho = exp(w) * coeff * erf(y) - y - 2.0 * pow(y, 3.0) / 3.0;
    return(Frho);
}



#define MAXLIN 256
#define MAXCOL 64

int read_file (instr, radius, dens, poten, emint, thdis, nrad, b)
stream instr;
double radius[];
double dens[];
double poten[];
double emint[];
double thdis[];
int nrad;
double b;
{
    int i, n, nlines, nval, radcol, denscol, masscol;
    char line[MAXLIN];
    double dval[MAXCOL];

    radcol = getiparam("radcol");
    denscol = getiparam("denscol");
    masscol = getiparam("masscol");
    if (radcol==0)
        error ("Must have radius as input");
    if (denscol==0 && masscol==0)
        error ("Must have either density or mass as input");
    if (denscol>0 && masscol>0)
        error ("Cannot have both density and mass as input");

    n = nrad;          /* maximum number to read from file */
    nlines=0;           /* line counter in file */
    for(;;) {
        if (!get_line(instr, line)) {
            nrad = nlines;
            break;
        } else if (nlines >= n) {
            dprintf(0,"Warning: declared (nrad=%d) space exhausted\n",n);
            break;
        }
	nval = nemoinpd(line,dval,MAXCOL);
        if (nval>MAXCOL || nval<0)
            error("Too many numbers on a line or error in parsing line: %s",line);

        if (nval<radcol)
            error("radius column referenced outside range");
        radius[nlines] = dval[radcol-1];

        if (nval<denscol)
            error("density column referenced outside range");
        if (denscol>0)
           dens[nlines] = dval[denscol-1];
        if (masscol>0)
           emint[nlines] = dval[masscol-1];

        nlines++;
    }
    
    if (masscol==0)
        calcpot_mass(radius, dens, poten, emint, nrad);
    if (denscol==0)
        calcpot_dens(radius, dens, poten, emint, nrad);
    
    for (i=0; i<n; i++)
        thdis[i] = 0.0;     /* cannot be known from this general procedure */
        
    printf ("Model read from file:  b=%f\n",b);
    printf ("    read %d lines from file\n",nrad);
    printf ("    radcol=%d denscol=%d masscol=%d\n",radcol, denscol, masscol);
}


/*
 * DVERK: differential euqation solver - calling odeint in (double) NumRec
 *
 */
 
#define NMAX 10

int dverk (n, fcn, x, y, xend, tol, ind, c, nw, w, ier)
double *x, *y, *xend, *tol, *c, *w;
int *n, *ind, *nw, *ier;
int (*fcn)();
{
    double ystart[NMAX], x1, x2, eps, h1, hmin, nok, nbad;
    int rkqc();
    int i, nvar;

    if (*ind != 1)
		error("dverk: IMSL simulator requires IND=1\n");
    if (*nw != *n)
		error("dverk: IMSL simulator requires NW=N?\n");
    if (*n > NMAX)
		error("dverk: IMSL simulator, recompile with NMAX = %d\n",*n); 

    nvar = *n;
    x1 = *x;
    x2 = *xend;
    for (i=0; i<*n; i++)
		ystart[i] = y[i];
    eps = *tol;
    h1 = (x2-x1)/20.0;		/* try 20 steps */
    hmin = 0.0;

    odeint(ystart,n,x1,x2,eps,h1,&nok,&nbad,fcn,rkqc);
    *x = *xend;

    *ier = 0;
    return(0);
}
