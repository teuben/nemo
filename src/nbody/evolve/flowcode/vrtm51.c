/* 
 * VRT: flow potential : derives vx and vy velocities given vr and vt.
 *
 *  Example use:
 *
 *   potname=vrtm51 
 *   potfile=tab1,tab2,...,tabN
 *   potpars=omega,pitch_angle,radius_ref,theta_ref,r0,r1,...,rN
 *
 *   potfiles are ascii tables with <angle(deg),vr,vt,den
 *
 *
 *      ** ONLY WORKS FOR LOG SPIRALS **
 *
 *	18-nov-03  cloned off vrt.c for the M51 project     Rahul & Peter
 *      19-nov-03  added phase equations to account for vr and vt being
 *                 constant along logarithmic arms.
 *      24-nov-03  allow multiple input files w/ added error checking
 *      25-nov-03  ieck, another emberassing reference frame error
 *      26-nov-03  indexing error in binsearch; sign error in thetaref
 *       1-jan-04  option to use constant vt in ring after omega*r
 *      29-dec-04  cleanup a valgrind complaint
 */

#include <stdinc.h>
#include <filestruct.h>
#include <image.h>
#include <table.h>
#include <extstring.h>
#include <spline.h>

#define VERSION "flowcode:vrtm51 V1.8b 21-jul-04"

local double omega = 0.0;		/*   pattern speed  */
local double pitch = 10.0;              /*    pitch angle   */
local double rref = 1.0;                /* reference radius */
local double thetaref = 0.0;            /*  reference angle */

#define MAXTAB  64

local int nrad[MAXTAB], nmax;
local real *theta[MAXTAB], *vr[MAXTAB], *vt[MAXTAB], *den[MAXTAB], 
  *coef_vr[MAXTAB], *coef_vt[MAXTAB], *coef_den[MAXTAB];
local real rings[MAXTAB+1];    /* rmin(0)... rmax(ntab+1), where rmin(i)=rmax(i-1) */
local int ntab=0;              /* number of tables used */
local int entries=0;           /* counter how often this routine was called */

local real     tanp;

local bool Qstick = TRUE;      /* set v=0 when hitting the inner or outer edge */
local bool Qconst = FALSE;     /* testing */


#define MAXCOL  4

extern int nemo_file_lines(string,int);
extern string *burststring(string,string);

/* 
 * binary search into an monotonically increasing array x or length n
 * returns 0 .. n, where 0 and n are outside the array
 *
 *
 *  x[0]   x[1]  ....   x[n-1]       <-- boundaries
 *    |      |             |
 *
 *  0 |   1  |   ....  n-1 |   n     <-- returned index
 */

local int binsearch(real u, real *x, int n)
{
  int i, j, k;
  
  if (u < x[0])                   /* check below left edge */
    return 0;
  else if (u >= x[n-1])           /* and above right edge */
    return n;
  else {
    i = 0;
    k = n;
    while (i+1 < k) {
      j = (i+k)/2;
      if (x[j] <= u)
	i = j;
      else
	k = j;
    }
    return i+1;
  }
}

void inipotential (int *npar, double *par, string name)
{
    int i, n, colnr[MAXCOL];
    real *coldat[MAXCOL];
    string *fnames;
    stream instr;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) pitch = par[1];
    if (n>2) rref = par[2];
    if (n>3) thetaref = par[3];

    if (entries>0) {
        warning("Re-entering inipotential(5NEMO): removed previous tables");
	for (i=0; i<ntab; i++) {
	  free(theta[i]);
	  free(vr[i]);       free(coef_vr[i]);
	  free(vt[i]);       free(coef_vt[i]);
	  free(den[i]);      free(coef_den[i]);
	}
    }
    entries++;

    tanp = tan(pitch*PI/180.0);

    dprintf (1,"INIPOTENTIAL vrt potential %s\n",name);
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  Parameters : Pitch Angle = %f\n",pitch);
    dprintf (1,"  Parameters : Reference Radius = %f\n",rref);
    dprintf (1,"  Parameters : Reference Angle = %f\n",thetaref);
    dprintf (1,"  Constant_vt: %s (hardcoded)\n",Qconst ? "TRUE" : "FALSE");
    dprintf (1,"  Sticky     : %s (hardcoded)\n",Qstick ? "TRUE" : "FALSE");
    dprintf (1,"  Table = %s\n",name);

    if (pitch == 0) error("inipotential: Need a non-zero pitch angle");

    fnames = burststring(name,",");
    ntab = xstrlen(fnames,sizeof(string))-1;
    if (ntab < 1) error("potpars= needs at least one filename");
    if (ntab > MAXTAB) 
      error("Not enough slots for input tables: %d > %d=MAXTAB",ntab,MAXTAB);
    if (ntab > 1) {
      if (ntab+5 != n) 
	error("found %d/%d , need potpars=ome,pitch,rs,ts,r_0,r_1,....r_%d\n",
	      n,ntab+5,ntab);
      for(i=0; i<=ntab; i++) {
	rings[i] = par[i+4];
	if (i>0 && rings[i] < rings[i-1])
	  error("ring radii in potpars= need to be sorted (%g > %g)",rings[i-1],rings[i-1]);
      }
    }

    for (i=0; i<ntab; i++) {      /* loop over all input files */

      nmax = nemo_file_lines(fnames[i],0);
      dprintf (1,"  Nmax = %d\n",nmax);
      if (nmax<=0) error("file_lines returned %d lines in %s\n",
			 nmax,fnames[i]);
      dprintf (1,"  Nmax = %d\n",nmax);
    
      coldat[0] = theta[i] = (real *) allocate(nmax * sizeof(real));
      coldat[1] = vr[i]    = (real *) allocate(nmax * sizeof(real));
      coldat[2] = vt[i]    = (real *) allocate(nmax * sizeof(real));
      coldat[3] = den[i]   = (real *) allocate(nmax * sizeof(real));
      
      colnr[0] = 1;
      colnr[1] = 2;
      colnr[2] = 3;
      colnr[3] = 4;

      instr = stropen(fnames[i],"r");
      nrad[i] = get_atable(instr,MAXCOL,colnr,coldat,nmax);
      dprintf(0,"table has %d entries\n",nrad[i]);
      strclose(instr);

      if (nrad[i] == 0) error("No lines (%d)read from %s",nrad[i],fnames[i]);
      if (nrad[i] < 0) {
	nrad[i] = -nrad[i];
	warning("only read part of table");
      }
      coef_vr[i]  = (real *) allocate(nrad[i]*4*sizeof(real));
      coef_vt[i]  = (real *) allocate(nrad[i]*4*sizeof(real));
      coef_den[i] = (real *) allocate(nrad[i]*4*sizeof(real));

      spline(coef_vr[i], theta[i], vr[i], nrad[i]);
      spline(coef_vt[i], theta[i], vt[i], nrad[i]);
      spline(coef_den[i], theta[i], den[i], nrad[i]);
      
      dprintf(2,"vrtm51[1]: %g %g %g %g\n",theta[i][0],vr[i][0],vt[i][0],den[i][0]);
      dprintf(2,"vrtm51[%d]: %g %g %g %g\n",
	      nrad[i],   theta[i][nrad[i]-1],
	                 vr[i][nrad[i]-1],
	                 vt[i][nrad[i]-1],
	                 den[i][nrad[i]-1]);
    }
}

#define  UNWRAP(p)    ((p) - 360 * floor((p)/ 360 ))

void potential(int *ndim,double *pos,double *acc,double *pot,double *time)
{
    real rad, rad0, rad1, phi, phase;
    real x, y, vrad, vtan, tmp;
    int i;

    x = pos[0];
    y = pos[1];
    rad = sqrt(x*x + y*y);

    phi = atan2(y,x)*180./PI;
    tmp = log(rad/rref)/tanp * 180/PI;
    phase = phi + tmp - thetaref;
    phase = UNWRAP(phase);        /* make sure 0..360 */

    if (ntab > 1) {
      i = binsearch(rad,rings,ntab+1);
      if (Qstick) {
	if (i==0 || i==ntab+1) {
	  acc[0] = acc[1] = acc[2] = 0.0;
	  *pot = 1.0;
	  return;
	}
      }
      if (i>0) i--;   /* binsearch gives one more than we want */
      if (i==ntab) i--;  /* outside outer ring */
      /* inside inner ring always point to inner ring */
      /* outside outer ring always point to outer ring */
    } else {
      /* if only one ring, all radii are valid, even if bad */
      i = 0;
    }

#if 0
    /* hmmm.... valgrind is actually failing the next stmt */
    rad1 = Qconst ? rad0 : rad ;
#else
    rad1 = rad;
#endif
      
    vrad = seval(phase,theta[i],vr[i],coef_vr[i],nrad[i]);
    vtan = seval(phase,theta[i],vt[i],coef_vt[i],nrad[i]) - omega*rad1;
    *pot = seval(phase,theta[i],den[i],coef_den[i],nrad[i]);
    dprintf(1,"x,y,rad,iring,phi,phase,DELTA,vt,vr,den: %g %g   %g %d %g %g [%g]  %g %g %g\n", 
	    x,y,rad,i,phi,phase,tmp,vtan,vrad,*pot);

    if (rad > 0) {
        acc[0] = (vrad * x - vtan * y) / rad;
        acc[1] = (vrad * y + vtan * x) / rad;
    } else {
        acc[0] = acc[1] = 0.0;
    }
    dprintf(1,"r,t=%g %g  vr,vt=%g %g vx,vy=%g %g\n",
            rad,phi,vrad,vtan,acc[0],acc[1]);
    acc[2] = 0.0;
}


