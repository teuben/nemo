/*
 * VELMAP: 
 *
 *	VELMAP performs various mapping functions on a projected
 *      galactic velocity field.
 *      
 *    3feb11  pjt   New task, cloned off velfit, for Nurur
 *
 */

#include <nemo.h>
#include <image.h>

string defv[] = {
  "in=???\n       input velocity field",
  "out=???\n      output map",
  "radii=\n       radii of the ring boundaries (Nring+1)",
  "pa=0\n         position angle of disk",
  "inc=45\n       inclination angle of disk",
  "center=\n      rotation center (mapcenter if left blank, 0,0=lower left)",
  "vsys=0\n       systemic velocity",
  "den=\n         input density image, unity if left blank",
  "frang=0\n      free angle around minor axis (2*frang is the total)",
  "blank=0.0\n    Value of the blank pixel to be ignored",
  "coswt=1\n      power of cos(theta) weighting",
  "VERSION=0.1\n  3-feb-2011 PJT",
  NULL,
};

string usage="various mapping functions on velocity fields";

string cvsid="$Id$";


#define MAXRING    1000


bool Qden = FALSE;
bool Qout = FALSE;
imageptr denptr = NULL, velptr = NULL, outptr = NULL;

real rad[MAXRING];
int nrad;

real pixe[MAXRING], flux[MAXRING], vsum[MAXRING], vsqu[MAXRING], wsum[MAXRING];
real vrot[MAXRING];

real pa, inc, vsys, xpos, ypos;
real undf;

int ring_index(int n, real *r, real rad)
{
  int i;
  if (rad < r[0]) return -1;
  if (rad > r[n-1]) return -2;
  for (i=0;i<n;i++)
    if (rad >= r[i] && rad < r[i+1]) return i;
  error("ring_index: should never gotten here %g in [%g : %g]",
	rad,r[0],r[n-1]);
}

nemo_main()
{
  stream denstr, velstr, outstr;
  real center[2], cospa, sinpa, cosi, sini, cost, costmin, x, y, xt, yt, r, den;
  real vr, wt, frang, dx, dy, xmin, ymin, rmin, rmax, fsum, ave, tmp, rms;
  real sincosi, cos2i, tga, dmin, dmax, dval;
  int i, j, k, nx, ny, ir, nring, nundf, nout, nang, nsum, coswt;

  velstr = stropen(getparam("in"),"r");

  read_image(velstr,&velptr);
  nx = Nx(velptr);
  ny = Ny(velptr);

  if (hasvalue("den")) {
    Qden = TRUE;
    denstr = stropen(getparam("den"),"r");
    read_image(denstr,&denptr);
  } 

  if (hasvalue("out")) {
    warning("New out= option not well tested yet");
    Qout = TRUE;
    outstr = stropen(getparam("out"),"w");
    copy_image(velptr,&outptr);
  } 

  nrad = nemoinpd(getparam("radii"),rad,MAXRING);
  if (nrad < 2) error("got no rings, use radii=");
  nring = nrad-1;
  if (hasvalue("center")) {
    if (nemoinpd(getparam("center"),center,2) != 2)
      error("not enuf for center=");
    xpos = center[0];
    ypos = center[1];
  } else {
    xpos = (Nx(velptr)-1.0)/2.0;
    ypos = (Ny(velptr)-1.0)/2.0;
  }
  pa = getdparam("pa");
  inc = getdparam("inc");
  vsys = getdparam("vsys");
  undf = getdparam("blank");
  frang = getdparam("frang");
  coswt = getiparam("coswt");

  cospa = cos(pa*PI/180.0);
  sinpa = sin(pa*PI/180.0);
  sini = sin(inc*PI/180.0);
  cosi = cos(inc*PI/180.0);
  costmin = sin(frang*PI/180.0);
  sincosi = sini*cosi;
  cos2i = cosi*cosi;
    
  for (i=0; i<nring; i++)
    pixe[i] = flux[i] = vsum[i] = vsqu[i] = wsum[i] = 0.0;
  nundf = nout = nang = 0;

  ymin = Ymin(velptr);
  xmin = Xmin(velptr);
  dx = Dx(velptr);       dx = ABS(dx);    dx = -dx;
  dy = Dy(velptr);       dy = ABS(dy);
  rmin = nx*dx*10.0;
  rmax = 0.0;

  dmin = dmax = vsys;

  /* loop over the map, accumulating data for fitting process */

  for (j=0; j<ny; j++) {
    y = (j-ypos)*dy;
    for (i=0; i<nx; i++) {
      if (MapValue(velptr,i,j) == undf) {
	nundf++;
	if (Qout) MapValue(outptr,i,j) = undf;
	continue;
      }
      x = (i-xpos)*dx;
      yt = x*sinpa + y*cospa;
      xt = x*cospa - y*sinpa;
#if 0
      if (Qout) {
	if (xt==0.0) {
	  MapValue(outptr,i,j) = undf;
	} else {
	  tga = xt/yt; /* X and Y are reversed  */
	  dval = MapValue(velptr,i,j)-vsys;
	  if (dval < 0) dval = -dval;
	  dval = vsys + dval*sqrt(cos2i+tga*tga)/sincosi;
	  MapValue(outptr,i,j) = dval;
	  if (dval > dmax) dmax = dval;
	  if (dval < dmin) dmin = dval;
	}
      }
#endif
      xt /= cosi;
      r  = sqrt(xt*xt+yt*yt);
      rmin = MIN(r,rmin);
      rmax = MAX(r,rmax);
      ir = ring_index(nrad,rad,r);
      dprintf(2,"r=%g ir=%d  (x,y)=%g,%g  (xt,yt)=%g,%g\n",
	      r,ir,x,y,xt,yt);
      if (ir < 0) {
	nout++;
	continue;
      }
      cost = yt/r;
      vr = (MapValue(velptr,i,j)-vsys)/cost/sini;

      dval = vr/r;
      MapValue(outptr,i,j) = dval;
      if (dval > dmax) dmax = dval;
      if (dval < dmin) dmin = dval;

      den = Qden ? MapValue(denptr,i,j) : 1.0;
      wt = den;
      for (k=0; k<coswt; k++)  wt *= cost;
      wt = ABS(wt);
      if (ABS(cost) > costmin) {
	pixe[ir] += 1.0;
	flux[ir] += den;
	vsum[ir] += wt*vr;
	vsqu[ir] += wt*vr*vr;
	wsum[ir] += wt;
      } else
	nang++;
    }
  }

  /* write output map(s), if needed */
  if (Qout) {
    dprintf(0,"Data min/max = %g %g\n",dmin,dmax);    
    MapMin(outptr) = dmin;
    MapMax(outptr) = dmax;
    write_image(outstr,outptr);
  }


  /* find the rotation curve */


  for (i=0; i<nring; i++) {
    if (pixe[i] < 1) continue;
    if(wsum[i] != 0)
      vrot[i] = vsum[i]/wsum[i];
    else
      vrot[i] = 0.0;
    /* reset arrays for next round of accumulations */
    wsum[i] = vsum[i] = vsqu[i] = 0.0;
  }


  /* loop over the map again, computing the residuals */

  for (j=0; j<ny; j++) {
    y = (j-ypos)*dy;
    for (i=0; i<nx; i++) {
      if (MapValue(velptr,i,j) == undf) continue;
      x = (i-xpos)*dx;
      yt =  x*sinpa + y*cospa;
      xt = (x*cospa - y*sinpa)/cosi;
      r  = sqrt(xt*xt+yt*yt);
      rmin = MIN(r,rmin);
      rmax = MAX(r,rmax);
      ir = ring_index(nrad,rad,r);
      if (ir < 0) continue;
      cost = yt/r;
      den = Qden ? MapValue(denptr,i,j) : 1.0;
      wt = den;
      for (k=0; k<coswt; k++)  wt *= cost;
      wt = ABS(wt);
      if (ABS(cost) > costmin) {
	vr = MapValue(velptr,i,j)-vsys-vrot[ir]*cost*sini;
	wsum[ir] += wt;
	vsum[ir] += wt*vr;
	vsqu[ir] += wt*vr*vr;
	dprintf(1,"%d %g %g %g %g\n",ir,xt,yt,vr,wt);
      }
    }
  }

  /* report */

  fsum = 0.0;
  nsum = 0;
  for (i=0; i<nring; i++) {
    if (wsum[i] == 0.0) continue;
    nsum++;
    r = 0.5*(rad[i] + rad[i+1]);
    tmp = flux[i]/pixe[i];
    ave = vsum[i]/wsum[i];
    rms = vsqu[i]/wsum[i]-ave*ave;
    if (rms <  0) rms=0.0;
    rms = sqrt(rms);
    fsum += rms*rms;

    dprintf(1,"%g %g %g %g %g %g ;; %g %g %g\n",
	   r,vrot[i],rms,pixe[i],tmp,sqrt(fsum/nsum),    ave,vsum[i],wsum[i]);
  }

  dprintf(0,"Nundf=%d/%d Nout=%d Nang=%d (sum=%d)\n",
	  nundf,nx*ny,nout,nang,nout+nundf+nang);
  dprintf(0,"Rmin/max = %g %g\n",rmin,rmax);

}

