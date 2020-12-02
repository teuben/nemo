/*
 * CCDELLINT:   integrate properties in elliptical rings
 *              cloned off VELMAP
 *      
 *    30-nov-2020   0.1    New task, cloned off velmap     PJT
 *
 *
 */

#include <nemo.h>
#include <image.h>

string defv[] = {
  "in=???\n       input velocity field",
  "radii=\n       radii of the ring boundaries (Nring+1)",
  "pa=0\n         position angle of disk",
  "inc=45\n       inclination angle of disk",
  "vsys=0\n       systemic velocity",
  "center=\n      rotation center (mapcenter if left blank, 0,0=lower left)",
  "blank=0.0\n    Value of the blank pixel to be ignored",
  "norm=t\n       Normalize RV image to number of pixels in ring",
  "out=\n         RV image",
  "tab=\n         Optional output table",
  "VERSION=0.2\n  1-dec-2020 PJT",
  NULL,
};


string usage="integrate map/cube in elliptical rings";

string cvsid="$Id$";


#ifndef MAXRING
#define MAXRING    2048
#endif

bool Qtab = FALSE;
imageptr denptr = NULL, velptr = NULL, outptr = NULL;

real rad[MAXRING];
int nrad;

int  pixe[MAXRING];
real vsum[MAXRING], vsqu[MAXRING], wsum[MAXRING];
real vrot[MAXRING];

real pa, inc, vsys, xpos, ypos;
real undf;

int string_index(string options, string s)
{
  int i=0;
  string *sa = burststring(options,",");

  while (sa[i]) {
    if (streq(sa[i],s))  {
      freestrings(sa);
      return i;
    }
    i++;
  }
  freestrings(sa);	   
  return -1;
}

int ring_index(int n, real *r, real rad)
{
  int i;
  if (rad < r[0]) return -1;
  if (rad > r[n-1]) return -2;
  for (i=0;i<n;i++)
    if (rad >= r[i] && rad < r[i+1]) return i;
  error("ring_index: should never gotten here %g in [%g : %g]",	rad,r[0],r[n-1]);
  return -1;
}

void nemo_main(void)
{
  stream denstr, velstr, outstr, tabstr;
  real center[2], cospa, sinpa, cosi, sini, x, y, xt, yt, r, den;
  real vr, wt, dx, dy, xmin, ymin, rmin, rmax, ave, tmp, rms;
  real sincosi, cos2i, tga, dmin, dmax, dval, dr, area, fsum1, fsum2;
  int i, j, k, nx, ny, nz, ir, nring, nundf, nout, nang, nsum;
  bool Qnorm = getbparam("norm");

  velstr = stropen(getparam("in"),"r");
  
  read_image(velstr,&velptr);
  nx = Nx(velptr);
  ny = Ny(velptr);
  nz = Nz(velptr);
  dprintf(0,"Image %d x %d x %d pixels = %g x %g x %g size\n",
	  nx,ny,nz, nx*Dx(velptr), ny*Dy(velptr), nz*Dz(velptr));


  if (hasvalue("tab")) {
    Qtab = TRUE;
    tabstr = stropen(getparam("tab"),"w");
  } 

  nrad = nemoinpd(getparam("radii"),rad,MAXRING);
  if (nrad < 2) error("got no rings (%d), use radii=",nrad);
  nring = nrad-1;

  create_image(&outptr, nring, nz);
  Xref(outptr) = 0.0;
  Xmin(outptr) = rad[0];
  Dx(outptr)   = rad[1] - rad[0];
  Yref(outptr) = Zref(velptr);
  Ymin(outptr) = Zmin(velptr);
  Dy(outptr)   = Dz(velptr);
  if (Namez(velptr)) {
    Namex(outptr)= strdup("deg");
    Namey(outptr)= strdup(Namez(velptr));
  } else {
    Namex(outptr)= strdup("R");
    Namey(outptr)= strdup("Z");
  }
    
  outstr = stropen(getparam("out"), "w");
    
  if (hasvalue("center")) {
    if (nemoinpd(getparam("center"),center,2) != 2)
      error("not enuf values for center=, need 2");
    xpos = center[0];
    ypos = center[1];
  } else {
    xpos = (Nx(velptr)-1.0)/2.0;
    ypos = (Ny(velptr)-1.0)/2.0;
  }
  pa    = getdparam("pa");
  inc   = getdparam("inc");
  vsys  = getdparam("vsys");
  undf  = getdparam("blank");

  cospa   = cos(pa*PI/180.0);
  sinpa   = sin(pa*PI/180.0);
  sini    = sin(inc*PI/180.0);
  cosi    = cos(inc*PI/180.0);
  sincosi = sini*cosi;
  cos2i   = cosi*cosi;
    
  for (i=0; i<nring; i++) 
    pixe[i] = 0.0;
  nundf = nout = nang = 0;

  ymin = Ymin(velptr);
  xmin = Xmin(velptr);
  dx = Dx(velptr);       dx = ABS(dx);    dx = -dx;
  dy = Dy(velptr);       dy = ABS(dy);
  rmin = -nx*dx*10.0;
  rmax = 0.0;
  dprintf(0,"Map %d x %d pixels, size %g x %g\n",
	  Nx(velptr), Ny(velptr), -dx*Nx(velptr), dy*Ny(velptr));
  dprintf(0,"Pixel size: %g x %g\n",-dx, dy);

  dmin = dmax = vsys;

  /* loop over the map, accumulating data */

  for (j=0; j<ny; j++) {
    y = (j-ypos)*dy;
    for (i=0; i<nx; i++) {
      x = (i-xpos)*dx;
      yt = x*sinpa + y*cospa;      /* major axis now along Y  */
      xt = x*cospa - y*sinpa;      /* minor axis along X      */
      xt /= cosi;                  /* deproject to the circle */
      r  = sqrt(xt*xt+yt*yt);      /* radius in the disk      */
      rmin = MIN(r,rmin);
      rmax = MAX(r,rmax);
      ir = ring_index(nrad,rad,r);
      dprintf(2,"r=%g ir=%d  (x,y)=%g,%g  (xt,yt)=%g,%g\n",
	      r,ir,x,y,xt,yt);
      if (ir < 0) {
	nout++;
	continue;
      }
      for (k=0; k<nz; k++) {
	if (CubeValue(velptr,i,j,k) == undf) {
	  nundf++;
	  continue;
	}
	dval = CubeValue(velptr,i,j,k);
	MapValue(outptr,ir,k) = MapValue(outptr,ir,k) + dval;
	pixe[ir]++;
      } /* k */
    } /* i */
  } /* j */

  if (Qnorm)
    for (i=0; i<nring; i++)
      if (pixe[i] > 0)
	for (k=0; k<nz; k++)
	  MapValue(outptr,i,k) = MapValue(outptr,i,k) / pixe[i];
  

  /* output R-V image */
  write_image(outstr, outptr);

  /* report on the rings */

  if (Qtab) {
    dprintf(0,"tab= not implemented yet");
  }
  dprintf(0,"Nundf=%d/%d Nout=%d Nang=%d (sum=%d)\n",
	  nundf,nx*ny,nout,nang,nout+nundf+nang);
  dprintf(0,"Rmin/max = %g %g\n",rmin,rmax);

}

