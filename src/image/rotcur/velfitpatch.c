/*
 * VELFITPATCH: 
 *
 *	VELFITPATCH fits a local velocity field in a small area
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
  "out=\n         basename for 3 output files (.vc,.ox,.oy) or (.vc,.om,.pa)",
  "patch=3\n      Patch (size of box will be 2*patch+1)",
  "scale=1\n      Extra scale factor for gradients",
  "blank=0.0\n    Blank Value in velocity field to skip",
  "mode=xy\n      Output mode (xy derivates, vs. gp (gradient/position angle)",
  "tab=\n         If given, tabular output file name",
  "VERSION=0.4\n  26-jul-2015 PJT",
  NULL,
};

string usage="fit local linear velocity field";

string cvsid="$Id$";


void fit_patch(int npt, real *x, real *y, real *v, real *sol, real scale);
void xy2rt(real *sol);

imageptr velptr = NULL, outptr0 = NULL, outptr1 = NULL, outptr2 = NULL;

nemo_main()
{
  bool Qsample, Qtab, Qout, Qxy;
  stream velstr, tabstr, outstr0, outstr1, outstr2;
  real scale;
  real *x, *y, *v, sol[3];
  int i, j, id, jd, nx, ny, nd, d, d2;
  string outmode = getparam("mode");
  string bname, oname;
  int mode = -1;
  int npt = 0;

  /* output mode */

  Qtab = hasvalue("tab");
  Qxy  = (*outmode == 'x' || *outmode == 'X');

  scale = getdparam("scale");

  /* read velocity field */

  velstr = stropen(getparam("in"),"r");
  read_image(velstr,&velptr);
  nx = Nx(velptr);
  ny = Ny(velptr);

  Qout = hasvalue("out");
  if (Qout) {
    copy_image(velptr, &outptr0);
    copy_image(velptr, &outptr1);
    copy_image(velptr, &outptr2);
    bname = getparam("out");
    oname = (string) allocate(strlen(bname) + 5);
    sprintf(oname,"%s.%s", bname, "vc");
    outstr0 = stropen(oname,"w");
    sprintf(oname,"%s.%s", bname, Qxy ? "ox" : "om");
    outstr1 = stropen(oname,"w");
    sprintf(oname,"%s.%s", bname, Qxy ? "oy" : "pa");
    outstr2 = stropen(oname,"w");
  }

  /* needed arrays */

  d = getiparam("patch");
  d2 = 2*d + 1;
  x = (real *) allocate(d2*d2*sizeof(real));
  y = (real *) allocate(d2*d2*sizeof(real));
  v = (real *) allocate(d2*d2*sizeof(real));

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
      fit_patch(nd, x, y, v, sol, scale);
      if (!Qxy) xy2rt(sol);
      if (Qtab) printf("%d %d   %g %g %g\n",i,j,sol[0],sol[1],sol[2]);
      npt++;
      if (Qout) {
	MapValue(outptr0,i,j) = sol[0];
	MapValue(outptr1,i,j) = sol[1];
	MapValue(outptr2,i,j) = sol[2];
      }
    }
  }
  dprintf(0,"Processed %d points with patch area %d x %d\n",npt,d2,d2);
  if (Qout) {
    write_image(outstr0,outptr0);
    write_image(outstr1,outptr1);
    write_image(outstr2,outptr2);
  }
}

/*
 * input:  v,dx,dy
 * output: v,dr,pa
 */

void xy2rt(real *sol)
{
  real dx = sol[1];
  real dy = sol[2];
  sol[1] = sqrt(dx*dx+dy*dy);
  sol[2] = atan2(dy,dx) * 180.0 / PI;
}


#define MAXCOL 4

void fit_patch(int npt, real *x, real *y, real *v, real *sol, real scale)
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
  if (scale != 1.0) {
    sol[1] *= scale;
    sol[2] *= scale;
  }
}

