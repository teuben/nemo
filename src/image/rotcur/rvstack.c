/*
 * RVSTACK:  Exact Radius-Velocity diagrams for later stacking
 *
 *  See also:     https://github.com/jbjolly/LineStacker
 *                jbjolly/LineStacker     Jolyy JV et al 2020
 *
 *    5-may-2021  PJT  Drafted
 *
 */

#include <nemo.h>
#include <image.h>
#include <spline.h>

string defv[] = {
  "in=???\n       input velocity field",
  "pa=0\n         position angle of receiding side of disk",
  "inc=45\n       inclination angle of disk",
  "vsys=0\n       systemic velocity to subtract",
  "center=\n      rotation center (mapcenter if left blank, 0,0=lower left)",
  "pan=90\n       position angle of near side of disk",
  "angle=10\n     (small) angle around major or minor axis",
  "blank=0.0\n    Value of the blank pixel to be ignored",
  "rscale=1\n     Scaling factor for radius",
  "gscale=f\n     Scale radius and velocity by the appropriate geometric sin/cos factor",
  "mode=rot\n     Rotation or Outflow",
  "VERSION=0.1\n  5-may-2021 PJT",
  NULL,
};

string usage="Exact Radius-Velocity diagrams for later stacking";

string cvsid="$Id$";


#define MAXRING    4096


bool Qwwb73;
bool Qden = FALSE;
bool Qout = FALSE;
bool Qtab = FALSE;
imageptr denptr = NULL, velptr = NULL, outptr = NULL;

real rad[MAXRING];
int nrad;

real pixe[MAXRING], flux[MAXRING], vsum[MAXRING], vsqu[MAXRING], wsum[MAXRING];
real vrot[MAXRING], radius[MAXRING], coeff[3*MAXRING];

real pa, pan, inc, vsys, xpos, ypos, zpos;
real undf;


void nemo_main()
{
  stream denstr, velstr, outstr, tabstr;
  real center[2], cospa, sinpa, cosi, sini, sint, cost, costmin, 
    x, y, v, xt, yt, r, den, vrot_s, dvdr_s, vmul, phi, dphi;
  real vr, wt, frang, dx, dy, dz, xmin, ymin, rmin, rmax, fsum, ave, tmp, rms;
  real angle = 0.5*getdparam("angle");
  real sincosi, cos2i, tga, dmin, dmax, dval, vmod;
  int i, j, k, nx, ny, nz, ir, nring, nundf, nout, nang, nsum, coswt;
  string outmode;
  int mode = -1;
  bool Qrot = TRUE;

  velstr = stropen(getparam("in"),"r");
  read_image(velstr,&velptr);

  nx = Nx(velptr);
  ny = Ny(velptr);
  nz = Ny(velptr);

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
  pan = getdparam("pan");
  inc = getdparam("inc");
  vsys = getdparam("vsys");
  undf = getdparam("blank");

  cospa   = cos(pa*PI/180.0);
  sinpa   = sin(pa*PI/180.0);
  sini    = sin(inc*PI/180.0);
  cosi    = cos(inc*PI/180.0);
  costmin = sin(frang*PI/180.0);
  sincosi = sini*cosi;
  cos2i   = cosi*cosi;
    
  nundf = nout = nang = 0;

  ymin = Ymin(velptr);
  xmin = Xmin(velptr);
  dx = Dx(velptr);       dx = ABS(dx);    dx = -dx;
  dy = Dy(velptr);       dy = ABS(dy);
  dz = Dz(velptr);

  dmin = dmax = vsys;
  zpos = nz/2.0;       // figure out where vsys is

  
  

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
      
      phi = atan2(y,x) * 180/PI;

      if (Qrot) {
	dphi = phi - pa;
	dphi = ABS(phi);
	if (dphi < angle)
	  vmul = 1.0;
	else {
	  continue;  // for now
	}
      } else {
	dphi = phi - pan;
	dphi = ABS(phi);
	if (dphi < angle)
	  vmul = 1.0;
	else {
	  continue; // for now
	}
      }

      r = sqrt(x*x + y*y);
      

      for (k=0; k<nz; k++) {
	v = (k-zpos)*dz - vsys;
	v *= vmul;
	printf("%g %g  %g    %d %d %g\n", r,v,CubeValue(velptr,i,j,k),i,j,dphi);
      } // k
    } // i
  } // j

}

