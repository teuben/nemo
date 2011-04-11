/*
 * VELFIT: (see also the MIRIAD program with the same name)
 *
 *	VELFIT fits a theoretical velocity pattern to an
 *	isovelocity image, weighted by the intensity image and a
 *	geometric factor specified by model input paramaters.
 *	The default option is to fit a rotation curve to an isovelocity
 *	image of a rotating disk. The rotation curve and rms for the
 *	fit are printed out and can be used to find the best fit to the
 *	other parameters. For more details see paper by Warner, Wright
 *	and Baldwin, 1973, MNRAS, 163,163.
 *
 *   30sep92  mchw  New task for Miriad.
 *   29aug02  pjt   Added frang= to prevent large divisions for models
 *   31aug02  pjt   Ieck, rms calculation wrong, arrays not reset to 0
 *    5-oct-2003  Created, embedding in rotcur was too painful    PJT
 *   22-oct-2005  Optional output map with rotation velocities  SNV/PJT
 *    3-feb-2011  more options to output map, added tab=           PJT
 *   10-feb-2011  cleanup and order of keywords aligned with rotcur PJT
 *
 *
 */

#include <nemo.h>
#include <image.h>
#include <spline.h>

string defv[] = {
  "in=???\n       input velocity field",
  "radii=\n       radii of the ring boundaries (Nring+1)",
  "pa=0\n         position angle of disk",
  "inc=45\n       inclination angle of disk",
  "vsys=0\n       systemic velocity",
  "center=\n      rotation center (mapcenter if left blank, 0,0=lower left)",
  "den=\n         input density image, unity if left blank",
  "frang=0\n      free angle around minor axis (2*frang is the total)",
  "blank=0.0\n    Value of the blank pixel to be ignored",
  "coswt=1\n      power of cos(theta) weighting",
  "wwb73=t\n      use the classic WWB73 method (fixed)",
  "mode=vtan\n    Output mode {vtan,vmod,vres,vtan/r,ome,vrad,dv/dr}",
  "out=\n         Optional output map of converted rotation speeds",
  "tab=\n         Optional output table of radii, velocities etc.",
  "VERSION=1.4\n  10-apr-2011 PJT",
  NULL,
};

string usage="fit rotation curve to coplanar disk (WWB73 method)";

string cvsid="$Id$";


#define MAXRING    1000


bool Qwwb73;
bool Qden = FALSE;
bool Qout = FALSE;
bool Qtab = FALSE;
imageptr denptr = NULL, velptr = NULL, outptr = NULL;

real rad[MAXRING];
int nrad;

real pixe[MAXRING], flux[MAXRING], vsum[MAXRING], vsqu[MAXRING], wsum[MAXRING];
real vrot[MAXRING], radius[MAXRING], coeff[3*MAXRING];

real pa, inc, vsys, xpos, ypos;
real undf;

    /* outmodes: the order of these is important, see mode= */
    /*    mode:    0    1    2    3      4   5    6     */
string outmodes = "vtan,vmod,vres,vtan/r,ome,vrad,dv/dr";

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

/*
 *  -1:     internal
 *   0:     first ring   (rad[0]..rad[1])
 *   nring: last ring
 *   -2:    outside
 */

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
  stream denstr, velstr, outstr, tabstr;
  real center[2], cospa, sinpa, cosi, sini, sint, cost, costmin, 
    x, y, xt, yt, r, den, vrot_s, dvdr_s;
  real vr, wt, frang, dx, dy, xmin, ymin, rmin, rmax, fsum, ave, tmp, rms;
  real sincosi, cos2i, tga, dmin, dmax, dval, vmod;
  int i, j, k, nx, ny, ir, nring, nundf, nout, nang, nsum, coswt;
  string outmode;
  int mode = -1;

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
    outmode = getparam("mode");
    mode = string_index(outmodes, outmode);
    if (mode < 0) error("Illegal mode=%s [%d], valid:",outmode,mode,outmodes);
    warning("New out= mode mode=%s [%d]",outmode,mode);
    Qout = TRUE;
    outstr = stropen(getparam("out"),"w");
    copy_image(velptr,&outptr);
  } 

  if (hasvalue("tab")) {
    Qtab = TRUE;
    tabstr = stropen(getparam("tab"),"w");
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

  cospa   = cos(pa*PI/180.0);
  sinpa   = sin(pa*PI/180.0);
  sini    = sin(inc*PI/180.0);
  cosi    = cos(inc*PI/180.0);
  costmin = sin(frang*PI/180.0);
  sincosi = sini*cosi;
  cos2i   = cosi*cosi;
    
  for (i=0; i<nring; i++)
    pixe[i] = flux[i] = vsum[i] = vsqu[i] = wsum[i] = 0.0;
  nundf = nout = nang = 0;

  ymin = Ymin(velptr);
  xmin = Xmin(velptr);
  dx = Dx(velptr);       dx = ABS(dx);    dx = -dx;
  dy = Dy(velptr);       dy = ABS(dy);
  rmin = -nx*dx*10.0;
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
      if (Qout) {      /* record the circular speed in diskplane */
	dval = MapValue(velptr,i,j)-vsys;
	if (dval < 0) dval = -dval;
	if (yt==0.0) {
	  dval = undf;
	} else {
	  tga = xt/yt;        /* X and Y are reversed  */
	  dval *= sqrt(cos2i+tga*tga)/sincosi;
	}
	if (mode==0)                      /* mode=vtan */
	  MapValue(outptr,i,j) = dval;
	else if (mode==3 || mode==4) {    /* mode=vtan/r or mode=ome */
	  MapValue(outptr,i,j) = dval/sqrt(sqr(xt/cosi)+sqr(yt));
	} 
      }
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
      vr = (MapValue(velptr,i,j)-vsys)/cost/sini;  /* rotation speed */
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

  /* set up rotation curve for differentation and get a spline */
  for (i=0; i<nring; i++) {
    radius[i] = 0.5*(rad[i]+rad[i+1]);
  }
  spline(coeff, radius, vrot, nring);


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
      sint = xt/r;
      if (ABS(cost) < costmin) {
	if (Qout) MapValue(outptr,i,j) = undf;
	continue;
      }
      den = Qden ? MapValue(denptr,i,j) : 1.0;
      wt = den;
      for (k=0; k<coswt; k++)  wt *= cost;
      wt = ABS(wt);
      if (ABS(cost) > costmin) {
	vrot_s = seval(r, radius, vrot, coeff, nring);  
	dvdr_s = spldif(r, radius, vrot, coeff, nring);
	vmod = vsys + vrot[ir]*cost*sini;     /* model */
	vr = MapValue(velptr,i,j)-vmod;       /* residual:   vobs-vmod */
	if (Qout) {
	  if (mode==1)                        /* mode=vmod */
	    MapValue(outptr,i,j) = vmod;
	  else if (mode==2)                   /* mode=vres */
	    MapValue(outptr,i,j) = vr;
	  else if (mode==5)                   /* mode=vrad */
	    MapValue(outptr,i,j) = vr/(sint*sini);
	  else if (mode==6)                   /* mode=dv/dr */
	    MapValue(outptr,i,j) = dvdr_s;
	}
	wsum[ir] += wt;
	vsum[ir] += wt*vr;
	vsqu[ir] += wt*vr*vr;
	dprintf(1,"%d %g %g %g %g\n",ir,xt,yt,vr,wt);
      } 
    }
  }

  /* report */

  if (Qtab) fprintf(tabstr,"# r v rms N tmp sqrt(F/N),  ave, vsum, wsum\n");
  fsum = 0.0;
  nsum = 0;
  for (i=0; i<nring; i++) {
    if (wsum[i] == 0.0) continue;
    nsum++;
    r = 0.5*(rad[i] + rad[i+1]);
    tmp = flux[i]/pixe[i];
    ave = vsum[i]/wsum[i];
    rms = vsqu[i]/wsum[i]-ave*ave;
    if (pixe[i]==0 || rms <  0) rms=0.0;
    rms = sqrt(rms);
    fsum += rms*rms;
    if (Qtab) fprintf(tabstr,"%g %g %g %g %g %g  %g %g %g\n",
	   r,vrot[i],rms,pixe[i],tmp,sqrt(fsum/nsum),    ave,vsum[i],wsum[i]);
  }

  dprintf(0,"Nundf=%d/%d Nout=%d Nang=%d (sum=%d)\n",
	  nundf,nx*ny,nout,nang,nout+nundf+nang);
  dprintf(0,"Rmin/max = %g %g\n",rmin,rmax);


  /* write output map(s), if needed */
  if (Qout) {
    dmin = dmax = MapValue(outptr,0,0); 
    for (j=0; j<ny; j++) {
      for (i=0; i<nx; i++) {
	dval = MapValue(outptr,i,j);
	if (dval > dmax) dmax = dval;
	if (dval < dmin) dmin = dval;
      }
    }
    dprintf(0,"Data min/max = %g %g\n",dmin,dmax);    
    MapMin(outptr) = dmin;
    MapMax(outptr) = dmax;
    write_image(outstr,outptr);
  }


}

