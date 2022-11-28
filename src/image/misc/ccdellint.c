/*
 * CCDELLINT:   integrate properties in elliptical rings
 *              cloned off VELMAP
 *      
 *    30-nov-2020   0.1    New task, cloned off velmap, for PPV -> RV     PJT
 *    23-nov-2022   0.4    Mods for 2D maps                               PJT
 *
 *
 */

#include <nemo.h>
#include <image.h>
#include <moment.h>

string defv[] = {
  "in=???\n       input velocity field",
  "radii=\n       radii of the ring boundaries (Nring+1)",
  "pa=0\n         position angle of disk",
  "inc=45\n       inclination angle of disk",
  "center=\n      rotation center (mapcenter if left blank, 0,0=lower left)",
  "vsys=0\n       systemic velocity (if PPV)",
  "blank=0.0\n    Value of the blank pixel to be ignored",
  "norm=f\n       Normalize RV image to number of pixels in ring",
  "out=\n         RV image",
  "tab=\n         Optional output table",
  "rscale=1\n     Scale applied ot radii",
  "iscale=1\n     Scale applied to intensities",
  "VERSION=0.5\n  25-nov-2022 PJT",
  NULL,
};


string usage="integrate map/cube in elliptical rings";

#ifndef MAXRING
#define MAXRING    2048
#endif

bool Qtab = FALSE;
bool Qout = FALSE;
imageptr velptr = NULL, outptr = NULL;

real rad[MAXRING];
int nrad;

int  pixe[MAXRING];

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
  stream velstr, outstr, tabstr;
  real center[2], cospa, sinpa, cosi, sini, x, y, xt, yt, r;
  real dx, dy, xmin, ymin, rmin, rmax, sum;
  real sincosi, cos2i, dmin, dmax, dval, dr, nppb;
  int i, j, k, nx, ny, nz, ir, nring, nundf, nout, nang;
  bool Qnorm = getbparam("norm");
  real iscale = getrparam("iscale");
  real rscale = getrparam("rscale");
  Moment *mp;

  velstr = stropen(getparam("in"),"r");
  
  read_image(velstr,&velptr);
  nx = Nx(velptr);
  ny = Ny(velptr);
  nz = Nz(velptr);
  dprintf(0,"Image %d x %d x %d pixels = %g x %g x %g size\n",
	  nx,ny,nz, nx*Dx(velptr), ny*Dy(velptr), nz*Dz(velptr));

  // Number of points per pixel = pi/(4*ln(2)) * (beam/pixel)^2
  nppb = 1.13309 * Beamx(velptr) * Beamy(velptr) / (Dx(velptr)*Dy(velptr));
  nppb = ABS(nppb);
  if (nppb==0) nppb=1;
  dprintf(0,"nppb = %g\n",nppb);


  if (hasvalue("tab")) {
    Qtab = TRUE;
    tabstr = stropen(getparam("tab"),"w");
    fprintf(tabstr,"#  nppb=%g rscale=%g iscale=%g\n",nppb,rscale,iscale);
    fprintf(tabstr,"# r npix int rms sum sumcum\n");
    dprintf(0,"#  nppb=%g rscale=%g iscale=%g\n",nppb,rscale,iscale);    
  }

  if (hasvalue("out")) {
    Qout = TRUE;
    outstr = stropen(getparam("out"), "w");
  }

  if (!Qtab && !Qout)
    warning("No output (out= or tab=) selected");

  nrad = nemoinpd(getparam("radii"),rad,MAXRING);
  if (nrad < 2) error("got no rings (%d), use radii=",nrad);
  nring = nrad-1;

  mp = (Moment *) allocate(nring * sizeof(Moment));
  for (i=0; i<nring; i++)
    ini_moment(&mp[i], 2, 0);

  create_image(&outptr, nring, nz);
  Xref(outptr) = 0.0;
  Xmin(outptr) = rad[0];
  Dx(outptr)   = dr = rad[1] - rad[0];
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

    
  if (hasvalue("center")) {
    if (nemoinpd(getparam("center"),center,2) != 2)
      error("not enuf values for center=%s, need 2",getparam("center"));
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
	if (nz==1) {    // for 2D maps keep track of RMS
	  accum_moment(&mp[ir], dval, 1.0);
	}
      } /* k */
    } /* i */
  } /* j */

  if (Qnorm) {
    for (i=0; i<nring; i++)
      if (pixe[i] > 0)
	for (k=0; k<nz; k++) {
	  MapValue(outptr,i,k) = MapValue(outptr,i,k) / pixe[i];
	  //printf("%d %d  %g\n",i,k,MapValue(outptr,i,k));
	}
  } 

  /* output R-V image */
  if (Qout)
    !write_image(outstr, outptr);

  /* report on the rings */

  if (Qtab) {
    if (Qnorm)      // un-normalize back if it has been
      for (i=0; i<nring; i++)
	if (pixe[i] > 0)
	  MapValue(outptr,i,0) = MapValue(outptr,i,0) * pixe[i];

    if (nz > 1) error("Cannot print table for 3D cube ellint");
    sum = 0.0;
    for (i=0; i<nring; i++) {
      r =  rad[i+1];
      sum += MapValue(outptr,i,0);
      // r npix int rms sum sumcum
      fprintf(tabstr,"%g %d %g %g %g %g\n",
	      rscale*r,
	      n_moment(&mp[i]),
	      iscale*mean_moment(&mp[i]),
	      iscale*sigma_moment(&mp[i]),
	      iscale*sum_moment(&mp[i]),
	      iscale*sum);
    }
    strclose(tabstr);
  }
  dprintf(0,"Nundf=%d/%d Nout=%d Nang=%d (sum=%d)\n",
	  nundf,nx*ny,nout,nang,nout+nundf+nang);
  dprintf(0,"Rmin/max = %g %g\n",rmin,rmax);

}

