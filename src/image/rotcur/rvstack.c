/*
 * RVSTACK:  exact Radius-Velocity diagrams for later stacking
 *
 *  See also:     https://github.com/jbjolly/LineStacker      Jolyy JV et al 2020
 *
 *    5-may-2021  PJT  Drafted, first working version
 *
 */

#include <nemo.h>
#include <image.h>

string defv[] = {
  "in=???\n       input cube",
  "out=???\n      output R-V-map",
  "pa=0\n         position angle of disk",
  "inc=45\n       inclination angle of disk",
  "vsys=0\n       systemic velocity to subtract",
  "center=\n      rotation center (mapcenter if left blank, 0,0=lower left)",
  "angle=10\n     (small) angle around major or minor axis",
  "blank=0.0\n    Value of the blank pixel to be ignored",
  "rscale=1\n     Scaling factor for radius",
  "vscale=1\n     Scaling factor for velocity",
  "gscale=f\n     Scale radius and velocity by the appropriate geometric sin/cos factor",
  "mode=rot\n     Rotation (r) or Outflow (o) ",
  "side=0\n       Both (0), or positive (1) or negative (-1) side",
  "tab=f\n        Write a test table?",
  "VERSION=0.3\n  6-may-2021 PJT",
  NULL,
};

string usage="exact Radius-Velocity diagrams for later stacking";

string cvsid="$Id$";


imageptr velptr = NULL, outptr = NULL, sumptr = NULL;

real pa, pan, inc, vsys, xpos, ypos, zpos;
real undf;


void nemo_main()
{
  stream denstr, velstr, outstr, tabstr;
  real center[2], cospa, sinpa, cosi, sini, sint, cost,
    x, y, v, xt, yt, r, vmul, phi, dphi;
  real vr, wt, frang, dx, dy, dz;
  real angle ;
  real rscale = getrparam("rscale");
  real vscale = getrparam("vscale");
  int  side   = getiparam("side");
  real dmin, dmax;
  int i, j, k, nx, ny, nz, nundf, nout, nang;
  int ir, iv;
  string mode;
  bool Qrot;
  bool gscale = getbparam("gscale");
  bool Qtab = getbparam("tab");
  bool Qran = FALSE;    // hack tom wiggle cube pixels to avoid blank map pixels when gscale=True

  velstr = stropen(getparam("in"),"r");
  read_image(velstr,&velptr);
  nx = Nx(velptr);
  ny = Ny(velptr);
  nz = Nz(velptr);

  outstr = stropen(getparam("out"),"w");

  mode = getparam("mode");
  Qrot = (*mode == 'r');

  if (hasvalue("center")) {
    if (nemoinpd(getparam("center"),center,2) != 2)
      error("not enuf for center=");
    xpos = center[0];
    ypos = center[1];
  } else {
    xpos = Xref(velptr);   // reference pixel
    ypos = Yref(velptr);
  }
  pa = getdparam("pa") * PI/180;
  inc = getdparam("inc") * PI/180;
  angle = getdparam("angle") * PI/180 * 0.5;
  vsys = getdparam("vsys");
  undf = getdparam("blank");

  cospa = cos(pa);
  sinpa = sin(pa);
  sini  = sin(inc);
  cosi  = cos(inc);
    
  nundf = nout = nang = 0;

  dx = Dx(velptr);    
  dy = Dy(velptr);    
  dz = Dz(velptr);

  zpos = (vsys - Zmin(velptr)) / Dz(velptr) + Zref(velptr);

  dprintf(0,"cube:   %d,%d,%d\n",nx,ny,nz);
  dprintf(0,"center: %g,%g,%g\n",xpos,ypos,zpos);

  create_image(&outptr, nx, nz);
  Xmin(outptr) = 0.0;
  Xref(outptr) = 0.0;
  Dx(outptr)   = ABS(Dx(velptr)) * rscale;
  Namex(outptr)= "Radius";
  Namey(outptr)= "Velocity";
  Ymin(outptr) = 0.0;
  Yref(outptr) = Ny(outptr)/2.0;
  Dy(outptr)   = ABS(Dz(velptr)) * vscale;
  create_image(&sumptr, nx, nz);
  dprintf(0,"Output map: %d x %d    %g %g %g x %g %g %g\n",
	  Nx(outptr), Ny(outptr),
	  Xmin(outptr), Dx(outptr), Xref(outptr),
	  Ymin(outptr), Dy(outptr), Yref(outptr));
	 

  for (j=0; j<ny; j++) {       // loop over all points in the map
    y = (j-ypos)*dy;
    if (Qran) y += grandom(0,dy);
    for (i=0; i<nx; i++) {
      x = (i-xpos)*dx;
      if (Qran) x += grandom(0,dy);      
      
      phi = atan2(-x,y);       // astronomical convention of PA

      if (Qrot) {              // rotation
	dphi = phi - pa;
	dphi = SYM_ANGLE(dphi);
	if (ABS(dphi) < angle && (side >=0) ) 
	  vmul = 1.0;
	else {
	  dphi += PI;
	  dphi = SYM_ANGLE(dphi);	  
	  if (ABS(dphi) < angle && (side <=0))
	    vmul = -1.0;
	  else
	    continue;
	}
      } else {                  // outflow
	dphi = phi - pa;
	dphi = SYM_ANGLE(dphi);	
	if (ABS(dphi) < angle && (side >=0)) 
	  vmul = 1.0;
	else {
	  dphi += PI;
	  dphi = SYM_ANGLE(dphi);	  
	  if (ABS(dphi) < angle && (side <= 0))
	    vmul = -1.0;
	  else
	    continue;
	}
      }
      
      for (k=0; k<nz; k++) {                  // loop over spectral point
	r = sqrt(x*x + y*y) * rscale;
	v = (k-zpos)*dz;
	if (Qran) v+= grandom(0,dz);
	v *= vmul;
	v *= vscale;
	if (gscale) {
	  if (Qrot)
	    v /= sini;
	  else {
	    r /= cosi;
	    v /= cosi;
	  }
	}
	ir = (int) (r/Dx(outptr));
	iv = (int) (((v-Ymin(outptr))/Dy(outptr)) + Yref(outptr));
	
	if (Qtab && (k == (int)zpos))
	  printf("%d   %g %g  %g    %d %d %g   %d %d\n",nang, r,v,CubeValue(velptr,i,j,k),i,j,dphi,ir,iv);

	if (ir >= Nx(outptr)) continue;
	if (iv >= Ny(outptr)) continue;
	if (ir < 0) continue;	
	if (iv < 0) continue;

	nang++;
	    
	MapValue(outptr,ir,iv) += CubeValue(velptr,i,j,k);
	MapValue(sumptr,ir,iv) += 1.0;
      } // k
    } // i
  } // j

  dmin = MapMax(velptr);
  dmax = MapMin(velptr);
  dprintf(0,"old Cube minmax: %g %g\n", dmax, dmin);
  
  for (j=0; j<Ny(outptr); j++)
    for (i=0; i<Nx(outptr); i++)
      if (MapValue(sumptr,i,j) > 0) {
	MapValue(outptr,i,j) /= MapValue(sumptr,i,j);
	if (MapValue(outptr,i,j) > dmax) dmax = MapValue(outptr,i,j);
	if (MapValue(outptr,i,j) < dmin) dmin = MapValue(outptr,i,j);
      }
  MapMin(outptr) = dmin;
  MapMax(outptr) = dmax;
  dprintf(0,"new Map  minmax: %g %g\n", dmin, dmax);  
  dprintf(0,"%d points selected\n", nang);
  write_image(outstr,outptr);

}

