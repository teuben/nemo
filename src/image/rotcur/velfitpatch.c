/*
 * VELFITPATCH: 
 *
 *	VELFITPATCH fits a local velocity field in a small area
 *	around a set of points taking along a (smooth) curve.
 *
 *   24-jul-2015   pjt    Written out specs 
 *   26-jul-2015   pjt    implemented - realized ccdtrace should probably do the curve part
 *
 */

#include <nemo.h>
#include <image.h>
#include <spline.h>

string defv[] = {
  "in=???\n       input velocity field",
  "out=\n         basename for 3 output files (vc,ox,oy) or (vc,om,pa)",
  "patch=3\n      Patch (size of box will be 2*patch+1)",
  "scale=1\n      Scale factor for gradients",
  "blank=0.0\n    Blank Value in velocity field",
  "mode=xy\n      Output mode (xy derivates, vs. gp (gradient/positionangle)",
  "tab=\n         If given, tabular output file name",
  "VERSION=0.2\n  26-jul-2015 PJT",
  NULL,
};

string usage="fit local linear velocity field (full map, or along a curve)";

string cvsid="$Id$";


#define MAXRING    4096

void fit_patch(int npt, real *x, real *y, real *v, real *sol);
void xy2rt(real *sol);

imageptr denptr = NULL, velptr = NULL, outptr = NULL;


nemo_main()
{
  bool Qsample, Qtab, Qxy;
  stream curstr, velstr, outstr, tabstr;
  real center[2], cospa, sinpa, cosi, sini, sint, cost, costmin, 
    xt, yt, r, den, vrot_s, dvdr_s;
  real vr, wt, frang, dx, dy, xmin, ymin, rmin, rmax, fsum, ave, tmp, rms;
  real sincosi, cos2i, tga, dmin, dmax, dval, vmod;
  real *x, *y, *v, sol[3];
  int i, j, id, jd, nx, ny, nd, d, d2;
  string outmode = getparam("mode");
  int mode = -1;
  int npt = 0;

  velstr = stropen(getparam("in"),"r");

  /* read velocity field */

  read_image(velstr,&velptr);
  nx = Nx(velptr);
  ny = Ny(velptr);

  /* needed arrays */

  d = getiparam("patch");
  d2 = 2*d + 1;
  x = (real *) allocate(d2*d2*sizeof(real));
  y = (real *) allocate(d2*d2*sizeof(real));
  v = (real *) allocate(d2*d2*sizeof(real));

  /* output mode */

  Qtab = hasvalue("tab");
  Qxy  = (*outmode == 'x' || *outmode == 'X');

  for (i=d; i<nx-d; i++) {
    for (j=d; j<ny-d; j++) {
      nd = 0;
      for (id=-d; id<=d; id++) {
	for (jd=-d; jd<=d; jd++) {
	  x[nd] = id;
	  y[nd] = jd;
	  v[nd] = MapValue(velptr, i+id, j+jd);
	  if (v[nd] != 0.0) nd++;
	}
      }
      fit_patch(nd, x, y, v, sol);
      if (!Qxy) xy2rt(sol);
      if (Qtab) printf("%d %d   %g %g %g\n",i,j,sol[0],sol[1],sol[2]);
      npt++;
    }
  }
  dprintf(0,"Processed %d points with patch area %d x %d\n",npt,d2,d2);


}

/*
 * input:  v,dx,dy
 * output: v,dr,phi
 */

void xy2rt(real *sol)
{
  real dx = sol[1];
  real dy = sol[2];
  sol[1] = sqrt(dx*dx+dy*dy);
  sol[2] = atan2(dy,dx) * 180.0 / PI;
}


#define MAXCOL 4

void fit_patch(int npt, real *x, real *y, real *v, real *sol)
{ 
  real mat[(MAXCOL+1)*(MAXCOL+1)], vec[MAXCOL+1], a[MAXCOL+2];
  bool Qpoly = FALSE;
  int i, order=2;

  lsq_zero(order+1, mat, vec);
  for (i=0; i<npt; i++) {
    a[0] = 1.0;
    a[1] = x[i];
    a[2] = y[i];
    a[3] = v[i];
    lsq_accum(order+1,mat,vec,a,1.0);
  }
  lsq_solve(order+1,mat,vec,sol);
}

