/*
 * MKHERNQUIST: 
 *
 * Generate a two-component, isotropic, spherical galaxy/halo model using a
 * superposition of two Hernquist models.  The galaxy is assumed to be a
 * unit Hernquist model with M=a=1 while the halo mass and scale radius are
 * free parameters specified on the command line.  The distribution
 * function is calculated on the fly using Eddington's formula and then
 * sampled to generate a N-body realization of a galaxy.  The galaxy and
 * halo particles are dumped in succession -- first half are galaxy
 * particles.  
 *
 *                  original version                John Dubinski - 1999
 *
 *      9-mar-99    Adopted for NEMO                Peter Teuben
 *                   - removed sorting and fcut stuff (snapsort,unbind)
 *                   - float->real (including the NR stuff; no query/ran1)
 *      9-sep-01    gsl/xrandom
 *     27-mar-03    fixed double/float usage or qromb (it never worked before)
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <history.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>  
#include <snapshot/body.h>
#include <snapshot/put_snap.c>

string  defv[] = {     
    "out=???\n                Output file name",
    "nbody=???\n              Total number of particles for both components",
    "m=1\n                    Mass of Halo (Galaxy=1)",
    "a=1\n                    Scale length of Halo (Galaxy=1)",
    "seed=0\n                 Standard seed",
    "addphi=f\n               Add potentials to snapshot",
    "zerocm=t\n               Centrate snapshot (t/f)?",
    "headline=\n              Optional verbiage",
    "VERSION=1.1\n            27-mar-03 PJT",
    NULL,
};

string usage="two-component, isotropic, spherical galaxy+halo Hernquist model";

#define CONST (1.0/sqrt(32.)/(M_PI*M_PI*M_PI))

int galaxy;
int nobj;
int *indx;
Body *bp, *btab;

int n; 
double dE;
double a, b, mu, E; 

/* local functions */

double rad(double eta), f(double E2), drho2_d(double u), pot(double r),
       radpsi(double p), radmass(double m), rho1(double r), rho2(double r);
void   cofm(void);
float  drho2(float u);

extern float qromb(float (*func)(float), float a, float b);

extern double xrandom(double,double);
extern int set_xrandom(int);



void nemo_main()
{
	int i, seed, bits;
	real mg, mh, tsnap;
	double vx, vy, vz;
	double fmax, f0, f1, v2, vmax, vmax2;
        bool Qcenter = getbparam("zerocm");
        bool Qphi = getbparam("addphi");
        stream outstr = stropen(getparam("out"),"w");

        mu = getdparam("m");
        a = getdparam("a");
        seed = init_xrandom(getparam("seed"));
        nobj = getiparam("nbody");
	if (nobj%2) 
	    warning("Total number of particles reset to %d\n",2*((nobj+1)/2));
        nobj = ((nobj+1)/2);
        if (hasvalue("headline"))
            set_headline(getparam("headline"));
        put_history(outstr);

	b = a - 1;

	btab = (Body *) allocate(2*nobj*sizeof(Body));

	for(i=0, bp=btab; i<2*nobj; i++, bp++) {
		double eta, radius, cth, sth, phi;
		double psi0;

		eta = (double) xrandom(0.0,1.0);
		radius = rad(eta);
		if( i >= nobj ) {
			if( mu == 0.0 ) break;
			radius *= a;
		}
		if( i<nobj ) {
			galaxy = 1;
		} else {
			galaxy = 0;
		}

		phi = xrandom(0.0,2*M_PI);
		cth = xrandom(-1.0,1.0);
		sth = sqrt(1.0 - cth*cth);
		Pos(bp)[0] = (real) (radius*sth*cos(phi));
		Pos(bp)[1] = (real) (radius*sth*sin(phi));
		Pos(bp)[2] = (real) (radius*cth);

		psi0 = -pot(radius);
                if (Qphi) Phi(bp) = psi0;
		vmax2 = 2.0*psi0;
		vmax = sqrt(vmax2);
		fmax = f(psi0);
		f0 = fmax; f1 = 1.1*fmax;        /* just some initial values */
		while( f1 > f0 ) {

                        /* pick a velocity */
			v2 = 2.0*vmax2;
			while( v2 > vmax2 ) {    /* rejection technique */
				vx = vmax*xrandom(-1.0,1.0);
				vy = vmax*xrandom(-1.0,1.0);
				vz = vmax*xrandom(-1.0,1.0);
				v2 = vx*vx + vy*vy + vz*vz;
			}

			f0 = f((psi0 - 0.5*v2));
			f1 = fmax*xrandom(0.0,1.0);

		}
		Vel(bp)[0] = vx;
		Vel(bp)[1] = vy;
		Vel(bp)[2] = vz;
        } 
	mg = 1.0/nobj;
	mh = mu/nobj;
	for(i=0, bp=btab; i<2*nobj; i++, bp++) {
            if (i<nobj)
                Mass(bp) = mg;
            else
                Mass(bp) = mh;
	}

	if (Qcenter)
	    cofm();

        bits = MassBit | PhaseSpaceBit;
        if (Qphi)  bits |= PotentialBit;
        nobj *= 2;
        tsnap = 0.0;
        put_snap(outstr, &btab, &nobj, &tsnap, &bits);
        strclose(outstr);
}
		
double rad(double eta)
{
	double sqeta;
	sqeta = sqrt(eta);
	return(sqeta/(1-sqeta));
}

void cofm()
{
	int i;
	real xcm, ycm, zcm;
	real vxcm, vycm, vzcm;
	real mtot;
	Body *bp;

	xcm = ycm = zcm = vxcm = vycm = vzcm = mtot = 0;
	for(i=0, bp=btab; i<2*nobj; i++, bp++) {
		xcm += Mass(bp)*Pos(bp)[0];
		ycm += Mass(bp)*Pos(bp)[1];
		zcm += Mass(bp)*Pos(bp)[2];
		vxcm += Mass(bp)*Vel(bp)[0];
		vycm += Mass(bp)*Vel(bp)[1];
		vzcm += Mass(bp)*Vel(bp)[2];
		mtot += Mass(bp);
	}

	xcm /= mtot;
	ycm /= mtot;
	zcm /= mtot;
	vxcm /= mtot;
	vycm /= mtot;
	vzcm /= mtot;
	dprintf(1,"Centering snapshot (%g,%g,%g) (%g,%g,%g)\n",
		xcm,ycm,zcm,vxcm,vycm,vzcm);


	for(i=0, bp=btab; i<2*nobj; i++, bp++) {
		Pos(bp)[0] -= xcm;
		Pos(bp)[1] -= ycm;
		Pos(bp)[2] -= zcm;
		Vel(bp)[0] -= vxcm;
		Vel(bp)[1] -= vycm;
		Vel(bp)[2] -= vzcm;
	}
}

/* Distribution function of the 2-component Hernquist model */

double f(double E2)
{
	double ans;

	E = E2;
	ans = qromb(drho2,0.0,2.0*sqrt(E))*CONST;

	return ans;
}

float drho2(float u)
{
  return (float)drho2_d( (double) u);
}

double drho2_d(double u)
{
	double ans;
	double r;
	double c0, c1, c02, c03, c12, c13, c04, c14;
	double r2, r3;
	double p;
	double dp1, dp2, drho1r, drho2r;

	p = E - 0.25*u*u;
	if( p <= 0 )
		return 0;
	else
		r = radpsi(p);

	r2 = r*r; r3 = r2*r;

	c0 = 1.0 + r;  c1 = a + r;
	c02 = c0*c0; c12 = c1*c1;
	c03 = c02*c0; c13 = c12*c1; c04 = c02*c02; c14 = c12*c12;
	dp1 = -1/c02 - mu/c12;
	dp2 = 2/c03 + 2*mu/c13;

	if( galaxy ) {
		drho1r = -(1+4*r)/(r2*c04);
		drho2r =  2.0*(10*r2 + 5*r + 1)/(r3*c04*c0);
	}
	else { /* its the halo */
		drho1r = - mu*a*(a + 4*r)/(r2*c14);
		drho2r = 2*mu*a*(10*r2 + 5*a*r + a*a)/(r3*c14*c1);
	}

	ans = drho2r/(dp1*dp1) - drho1r/(dp1*dp1*dp1)*dp2;

	return ans;
}

double pot(double r)
{
	double p1, p2;

	p1 = -1.0/(1.0 + r);
	p2 = -mu/(a + r);

	return p1+p2;
}


double radpsi(double p)
{
	double b0, c0;
	double p1;
	double ans;

	p1 = 1/p;
	b0 = 0.5*((1+a) - (1+mu)*p1);
	c0 = a*(1-p1) - mu*p1;

	ans = -b0 + sqrt(b0*b0 - c0);
	return ans;
}

double radmass(double m)
{
	double eps;
	double rguess, r0, r1;
	double dr;
	double c0, c1, c02, c12; 

	eps = 1.0e-6; dr = 1;
	rguess = 0;
	r0 = 0.001;
	while( dr > eps ) {
		c0 = 1 + r0; c1 = a + r0;
		c02 = c0*c0; c12 = c1*c1;
		r1 = r0 - (1/c02 + mu/c12 - m/(r0*r0))/
			(-2/(c0*c02) -2*mu/(c1*c12) + 2*m/(r0*r0*r0));
		dr = fabs(r1-r0);
		r0 = r1;
	}
	return r0;
}
		
double rho1(double r)
{
	return 1.0/r/pow(1+r,3.0);
}

double rho2(double r)
{
	return mu*a/r/pow(a+r,3.0);
}

