/* 
 * VRT: flow potential : derives vx and vy velocities given vr and vt.
 *
 *	18-nov-03  cloned off vrt.c for the M51 project     Rahul & Peter
 *      19-nov-03  added phase equations to account for vr and vt being
 *                 constant along logarithmic arms.
 */

#include <stdinc.h>
#include <filestruct.h>
#include <image.h>
#include <table.h>

#define VERSION "flowcode:vrtm51 V1.1 20-nov-03"

local double omega = 0.0;		/*   pattern speed  */
local double pitch = 10.0;              /*    pitch angle   */
local double rref = 1.0;               /* reference radius */
local double thetaref = 0.0;           /*  reference angle */

local int nrad, nmax;
local real *theta, *vr, *vt, *den, *coef_vr, *coef_vt, *coef_den;
local int entries=0;


local double   vscale = 1.0;		/* rescale velocity unit */
local double   rscale = 1.0;		/* rescale radius unit */

local stream   potstr = NULL;
local int      nr, np;
local real     *rads, *phis;
local real     dphi;
local real     tanp;


#define MAXCOL  4

extern int nemo_file_lines(string,int);
extern void spline (double *coef, double *x, double *y, int n);
extern double seval(double x0, double *x, double *y, double *coef, int n);


void inipotential (int *npar, double *par, string name)
{
    int n, colnr[MAXCOL];
    real *coldat[MAXCOL];
    stream instr;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) pitch = par[1];
    if (n>2) rref = par[2];
    if (n>3) thetaref = par[3];
    if (n>4) warning("inipotential(flowrt): npar=%d only 4 parameters accepted (ome,pitch,rref,tref)",n);
    
    if (entries>0) {
        warning("Re-entering inipotential(5NEMO): removed previous tables");
        free(theta);
        free(vr);
        free(vt);
	free(den);
    }
    entries++;

    tanp = tan(pitch*PI/180.0);

    dprintf (1,"INIPOTENTIAL vrt potential %s\n",name);
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  Parameters : Pitch Angle = %f\n",pitch);
    dprintf (1,"  Parameters : Reference Radius = %f\n",rref);
    dprintf (1,"  Parameters : Reference Angle = %f\n",thetaref);
    dprintf (1,"  Table = %s\n",name);

    if (pitch == 0) error("inipotential: Need a non-zero pitch angle");

    nmax = nemo_file_lines(name,0);
    dprintf (1,"  Nmax = %d\n",nmax);
    if (nmax<=0) error("file_lines returned %d lines in %s\n",
             nmax,name);
    dprintf (1,"  Nmax = %d\n",nmax);
    
    theta = (real *) allocate(nmax * sizeof(real));
    vr = (real *) allocate(nmax * sizeof(real));
    vt = (real *) allocate(nmax * sizeof(real));
    den = (real *) allocate(nmax * sizeof(real));

    coldat[0] = theta;        colnr[0] = 1;
    coldat[1] = vr;           colnr[1] = 2;
    coldat[2] = vt;           colnr[2] = 3;
    coldat[3] = den;          colnr[3] = 4;

    instr = stropen(name,"r");
    nrad = get_atable(instr,MAXCOL,colnr,coldat,nmax);
    dprintf(0,"table has %d entries\n",nrad);
    strclose(instr);
    if (nrad==0) error("No lines (%d)read from %s",nrad,name);
    if (nrad<0) {
	nrad = -nrad;
	warning("only read part of table");
    }
    coef_vr = (real *) allocate(nrad*4*sizeof(real));
    coef_vt = (real *) allocate(nrad*4*sizeof(real));
    coef_den = (real *) allocate(nrad*4*sizeof(real));

    spline(coef_vr, theta, vr, nrad);
    spline(coef_vt, theta, vt, nrad);
    spline(coef_den, theta, den, nrad);

    dprintf(2,"vrtm51[1]: %g %g %g %g\n",theta[0],vr[0],vt[0],den[0]);
    dprintf(2,"vrtm51[%d]: %g %g %g %g\n",nrad,theta[nrad-1],vr[nrad-1],vt[nrad-1],den[nrad-1]);
}

#define  UNWRAP(phi)    ((phi) - 360 * floor((phi)/ 360 ))


void potential(int *ndim,double *pos,double *acc,double *pot,double *time)
{
    real rad, phi, phase;
    real x, y, vrad, vtan;
    bool mirror;

    x = pos[0];
    y = pos[1];
    rad = sqrt(x*x + y*y);
    phi = atan2(y,x)*180./PI;
    phi -= log(rad/rref)/tanp + thetaref;
    phase = UNWRAP(phi);
      
    vrad = seval(phase,theta,vr,coef_vr,nrad);
    vtan = seval(phase,theta,vt,coef_vt,nrad);
    *pot = seval(phase,theta,den,coef_den,nrad);
    dprintf(1,"x,y,rad,phi,phase,vt,vr,den: %g %g   %g %g %g   %g %g %g\n", x,y,rad,phi,phase,vtan,vrad,*pot);

    if (rad > 0) {
	vtan -= omega * rad;
        acc[0] = (vrad * x - vtan * y) / (rad * vscale);
        acc[1] = (vrad * y + vtan * x) / (rad * vscale);
    } else {
        acc[0] = acc[1] = 0.0;
    }
    dprintf(1,"r,t=%g %g  vr,vt=%g %g vx,vy=%g %g\n",
            rad,phi,vrad,vtan,acc[0],acc[1]);
    acc[2] = 0.0;
}
