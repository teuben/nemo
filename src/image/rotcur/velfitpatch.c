/*
 * VELFITPATH: 
 *
 *	VELFITPATH fits a local velocity field in a small area
 *	around a set of points taking along a (smooth) curve.
 *
 *   24-jul-2015   pjt    Written out specs 
 *
 */

#include <nemo.h>
#include <image.h>
#include <spline.h>

string defv[] = {
  "in=???\n       input velocity field",
  "curve=???\n    Table with X and Y positions (in image coordinates) for curve",
  "spline=n\n     Spline fit along curve?",
  "step=1\n       Steps along curve",
  "box=10,10\n    Boxsize in X and Y (X along curve)",
  "VERSION=0.1\n  24-jul-2015 PJT",
  NULL,
};

string usage="fit local velocity field along a curve";

string cvsid="$Id$";


#define MAXRING    4096


imageptr denptr = NULL, velptr = NULL, outptr = NULL;


nemo_main()
{
  stream curstr, velstr, outstr, tabstr;
  real center[2], cospa, sinpa, cosi, sini, sint, cost, costmin, 
    x, y, xt, yt, r, den, vrot_s, dvdr_s;
  real vr, wt, frang, dx, dy, xmin, ymin, rmin, rmax, fsum, ave, tmp, rms;
  real sincosi, cos2i, tga, dmin, dmax, dval, vmod;
  int i, j, k, nx, ny, ir, nring, nundf, nout, nang, nsum, coswt;
  string outmode;
  int mode = -1;

  velstr = stropen(getparam("in"),"r");
  curstr = stropen(getparam("curve"),"r");

  /* read velocity field */

  read_image(velstr,&velptr);
  nx = Nx(velptr);
  ny = Ny(velptr);

  /* read in curve */

  /* optionally: make spline through curve */


  /* walk along the line in steps S */


  /* find the points in the box LX,LY around this point */

  /* fit a plane:
          v(x,y) = v0 + a*x + b*y
     where (x,y) are the coordinates w.r.t. the center of the box
  */


  /* print out results :
     Sx = center in X along curve
     Sy = center in Y along curve
     Sl = length along curve from starting point
     S_phi = PA of curve 
     S_vel = fitted PA of center velocity
     V0  = fitted center velocity
     dVdS = fitted velocity gradient
  */



}

