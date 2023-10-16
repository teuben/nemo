/*
 * NKING.C: construct an n-component King model.
 *
 *	27-oct-88	??          - Josh Barnes ?
 *	 7-may-92  PJT	new NEMO V2.
 */

#include <stdinc.h>
#include <getparam.h>

string defv[] = {
    "rho=1.0\n              density scale parameter(s)",
    "sig=1.0\n              velocity scale parameter(s)",
    "psi_0=1.0\n            scaled central potential",
    "step=0.05\n            step in terms of core radius",
    "tab=-\n                Output table",
    "VERSION=1.1\n          7-may-92 PJT",
    NULL,
};

string usage="construct an n-component King model.";

#define MAXCOMP  4			/* max. no. of components */

int ncomp;				/* actual no. of components */

real rhotab[MAXCOMP];			/* density scale parameters */
real sigtab[MAXCOMP];			/* velocity scale parameters */

real psi_0;				/* scaled central potential */

real step;				/* step in terms of core radius */

#define MAXVALS  512			/* max. no. of tabulated values */

int nvals;				/* number of tabulated values */

real radv[MAXVALS];			/* array of tabulated radii */
real psiv[MAXVALS];			/* array of relative potentials */
real massv[MAXCOMP][MAXVALS];		/* arrays of cumulative mass */

real r_core, m_core;			/* overall core radius and mass */

real Mtot, Ttot, Utot;			/* total mass and binding energy */

real m_comp[MAXCOMP];			/* total mass of each component */
real v_comp[MAXCOMP];			/* rms velocity of each component */

stream ostr;

void diffstep();

nemo_main()
{
    setparams();
    kingmodel();
}

setparams()
{
    string *burststring(), *rptr, *sptr;
    real atof();

    rptr = burststring(getparam("rho"), ", ");
    sptr = burststring(getparam("sig"), ", ");
    if (*rptr == NULL || *sptr == NULL)
	error("nking: rho or sig missing");
    ncomp = 0;
    while (*rptr != NULL || *sptr != NULL) {
	if (ncomp >= MAXCOMP)
	    error("nking: too many components; MAXCOMP=%d",MAXCOMP);
	rhotab[ncomp] = (*rptr != NULL ? atof(*rptr++) : rhotab[ncomp-1]);
	sigtab[ncomp] = (*sptr != NULL ? atof(*sptr++) : sigtab[ncomp-1]);
	ncomp++;
    }
    psi_0 = getdparam("psi_0");
    step = getdparam("step");
    ostr = stropen(getparam("tab"),"w");
}

kingmodel()
{
    real tmp, rho_0, sigsq, rho(), sqrt();
    int i;

    rho_0 = sigsq = 0.0;
    dprintf(2,"Model\trho\tsig\tpsi\n");
    for (i = 0; i < ncomp; i++) {
        tmp = rho(psi_0, i);
	rho_0 += tmp;
	sigsq += sigtab[i]*sigtab[i] / ncomp;
	dprintf(2,"%d\t%g\t%g\t%g\n",i+1,rhotab[i],sigtab[i],tmp);
    }
    dprintf(2,"\n");
    r_core = sqrt(9 * sigsq / (FOUR_PI * rho_0));
    m_core = FRTHRD_PI * r_core*r_core*r_core * rho_0;

    dprintf(1,"rho_0:  %12.5f\tsigsq:  %12.5f\n", rho_0, sigsq);
    dprintf(1,"r_core: %12.5f\tm_core: %12.5f\n", r_core, m_core);
    dprintf(1,"\n");

    kingsolve();
}

#define MAXVARS  (3 + 2*MAXCOMP)	/* max. no. of integ. vars. */

#define Rad(vv)     (vv[0])		/* radius (indp var) */
#define Psi(vv)     (vv[1])		/* relative potential */
#define Epot(vv)    (vv[2])		/* cumulative P.E. */
#define Mass(vv,i)  (vv[3+2*i])		/* mass of ith component */
#define Ekin(vv,i)  (vv[4+2*i])		/* K.E. of ith component */

#define DMASSMIN  1.0E-4		/* mass increment in terms of m_core */

kingsolve()
{
    real vv[MAXVARS], mlast, f;
    int i;
    void diffvv();

    Rad(vv) = 0.0;
    Psi(vv) = psi_0;
    Epot(vv) = 0.0;
    for (i = 0; i < ncomp; i++) {
	Mass(vv, i) = 0.0;
	Ekin(vv, i) = 0.0;
    }
    nvals = 0;
    mlast = -1.0;
    while (Psi(vv) > 0.0) {
	Mtot = 0.0;
	for (i = 0; i < ncomp; i++)
	    Mtot += Mass(vv, i);
	if (Mtot - mlast > DMASSMIN * m_core) {
	    if (nvals >= MAXVALS)
		error("nking: too many steps");
	    radv[nvals] = Rad(vv);
	    psiv[nvals] = Psi(vv);
	    for (i = 0; i < ncomp; i++) 
		massv[i][nvals] = Mass(vv, i);
            if (ostr) {
                fprintf(ostr,"%g %g",radv[nvals],psiv[nvals]);
              	for (i = 0; i < ncomp; i++) 
                    fprintf(ostr," %g",massv[i][nvals]);
                fprintf(ostr,"\n");
            }
	    nvals++;
	    mlast = Mtot;
	}
	diffstep(vv, 3 + 2*ncomp, diffvv, step * r_core);
    }
    f = psiv[nvals-1] / (psiv[nvals-1] - Psi(vv));

    dprintf(1,"nvals:%d\tf: %f\n", nvals, f);
    dprintf(1,"radv[%d] = %f (before)\n", nvals-1, radv[nvals-1]);
    radv[nvals-1] = radv[nvals-1] + f * (Rad(vv) - radv[nvals-1]);
    dprintf(1,"radv[%d] = %f (after)\n", nvals-1, radv[nvals-1]);
    psiv[nvals-1] = 0.0;
    Mtot = Ttot = 0.0;
    for (i = 0; i < ncomp; i++) {
	Mtot += massv[i][nvals-1];
	Ttot += Ekin(vv, i);
    }
    Utot = Epot(vv) + 0.5 * Mtot*Mtot / radv[nvals-1];
    dprintf(1,"Mtot: %f\tTtot: %f\tUtot: %f\tEpot: %f\n",
	   Mtot, Ttot, Utot, Epot(vv));
}    

void diffvv(dvv, vv)
real dvv[], vv[];
{
    real r, psi, rho(), eps(), M, dM;
    int i;

    r = Rad(vv);
    psi = Psi(vv);
    for (i = 0; i < ncomp; i++) {
	Mass(dvv, i) = FOUR_PI * r*r * rho(psi, i);
	Ekin(dvv, i) = FOUR_PI * r*r * eps(psi, i);
    }
    M = dM = 0.0;
    for (i = 0; i < ncomp; i++) {
	M += Mass(vv, i);
	dM += Mass(dvv, i);
    }
    Rad(dvv) = 1.0;
    Psi(dvv) = (M > 0.0 ? - M / r*r : 0.0);
    Epot(dvv) = 0.5 * dM * psi;
}

real rho(psi, i)
real psi;
int i;
{
    real sigsq, sqrt(), Sn();

    sigsq = sigtab[i] * sigtab[i];
    return sqrt(2 / PI) * rhotab[i] * Sn(psi / sigsq, 2);
}

real eps(psi, i)
real psi;
int i;
{
    real sigsq, sqrt(), Sn();

    sigsq = sigtab[i] * sigtab[i];
    return 1.5 * sqrt(2 / PI) * rhotab[i] * sigsq * Sn(psi / sigsq, 3);
}

real Sn(w, n)
real w;
int n;
{
    real sum, sqrt(), s;
    int i;

    sum = 0.0;
    if (w > 0.0) {
	s = sqrt(2 * w);
	i = 0;
	while (ABS(s) > 1.0e-7 * ABS(sum)) {
	    if (i >= n)
		sum += s;
	    i++;
	    s = 2 * w * s / (2*i + 1);
	}
    }
    return sum;
}

void diffstep(x0, nx, diff, delta)
real x0[];
int nx;
proc diff;
real delta;
{
    real dxi, x1[MAXVARS], dx0[MAXVARS], dx1[MAXVARS];
    int i;

    (*diff)(dx1, x0);
    for (i = 0; i < nx; i++) {
	dxi = 0.5 * delta * dx1[i];
	dx0[i] = dxi;
	x1[i] = x0[i] + dxi;
    }
    (*diff)(dx1, x1);
    for (i = 0; i < nx; i++) {
	dxi = delta * dx1[i];
	dx0[i] = dx0[i] + dxi;
	x1[i] = x0[i] + 0.5 * dxi;
    }
    (*diff)(dx1, x1);
    for (i = 0; i < nx; i++) {
	dxi = delta * dx1[i];
	dx0[i] = dx0[i] + dxi;
	x1[i] = x0[i] + dxi;
    }
    (*diff)(dx1, x1);
    for (i = 0; i < nx; i++)
	x0[i] = x0[i] + (dx0[i] + delta * dx1[i]) / 3;
}

