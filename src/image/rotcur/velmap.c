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
  "radii=\n       radii of the ring boundaries (Nring+1)",
  "pa=0\n         position angle of disk",
  "inc=45\n       inclination angle of disk",
  "center=\n      rotation center (mapcenter if left blank, 0,0=lower left)",
  "vsys=0\n       systemic velocity",
  "den=\n         input density image, unity if left blank",
  "frang=0\n      free angle around minor axis (2*frang is the total)",
  "blank=0.0\n    Value of the blank pixel to be ignored",
  "mode=v/r\n     Output mode (v, v/r, d*v/r, ...)",
  "out=\n         Optional output map",
  "tab=\n         Optional output table",
  "VERSION=0.3\n  8-feb-2011 PJT",
  NULL,
};


string usage="various mapping functions on velocity fields";

string cvsid="$Id$";


#ifndef MAXRING
#define MAXRING    2048
#endif

bool Qden = FALSE;
bool Qout = FALSE;
bool Qtab = FALSE;
imageptr denptr = NULL, velptr = NULL, outptr = NULL;

real rad[MAXRING];
int nrad;

real pixe[MAXRING], vsum[MAXRING], vsqu[MAXRING], wsum[MAXRING];
real vrot[MAXRING];

real pa, inc, vsys, xpos, ypos;
real undf;

               /*  0 1   2    */
string outmodes = "v,v/r,d*v/r";

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
  error("ring_index: should never gotten here %g in [%g : %g]",
	rad,r[0],r[n-1]);
}

nemo_main()
{
  stream denstr, velstr, outstr, tabstr;
  real center[2], cospa, sinpa, cosi, sini, cost, costmin, x, y, xt, yt, r, den;
  real vr, wt, frang, dx, dy, xmin, ymin, rmin, rmax, fsum, ave, tmp, rms;
  real sincosi, cos2i, tga, dmin, dmax, dval, dr, area;
  int i, j, k, nx, ny, ir, nring, nundf, nout, nang, nsum;
  string outmode;
  int mode = -1;

  velstr = stropen(getparam("in"),"r");

  read_image(velstr,&velptr);
  nx = Nx(velptr);
  ny = Ny(velptr);

  if (hasvalue("out")) {
    outmode = getparam("mode");
    mode = string_index(outmodes, outmode);
    if (mode < 0) error("Illegal mode=%s [%d], valid:",outmode,mode,outmodes);
    warning("New out= mode mode=%s [%d]",outmode,mode);
    Qout = TRUE;
    outstr = stropen(getparam("out"),"w");
    copy_image(velptr,&outptr);
  } 

  if (hasvalue("den")) {
    Qden = TRUE;
    denstr = stropen(getparam("den"),"r");
    read_image(denstr,&denptr);
  } else if (mode==2)
    error("Need den=");


  if (hasvalue("tab")) {
    Qtab = TRUE;
    tabstr = stropen(getparam("tab"),"w");
  } 

  nrad = nemoinpd(getparam("radii"),rad,MAXRING);
  if (nrad < 2) error("got no rings (%d), use radii=",nrad);
  nring = nrad-1;
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
  frang = getdparam("frang");

  cospa   = cos(pa*PI/180.0);
  sinpa   = sin(pa*PI/180.0);
  sini    = sin(inc*PI/180.0);
  cosi    = cos(inc*PI/180.0);
  costmin = sin(frang*PI/180.0);
  sincosi = sini*cosi;
  cos2i   = cosi*cosi;
    
  for (i=0; i<nring; i++)
    pixe[i] = vsum[i] = vsqu[i] = wsum[i] = 0.0;
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
      cost = yt/r;
      dval = MapValue(velptr,i,j);
      
      if (mode==1)
	dval /= r;
      else if (mode==2)
	dval *= MapValue(denptr,i,j) / r;

      if (outptr) {
	if  (ABS(cost) > costmin) {
	  MapValue(outptr,i,j) = dval;
	  if (dval > dmax) dmax = dval;
	  if (dval < dmin) dmin = dval;
	} else {
	  MapValue(outptr,i,j) = undf;
	}
      }

      /* now some ring accumulation, remnant of the velfit fitting */

      den = Qden ? MapValue(denptr,i,j) : 1.0;
      wt = den;
      wt = ABS(wt);
      if (ABS(cost) > costmin) {
	pixe[ir] += 1.0;
	wsum[ir] += wt;
	vsum[ir] += wt*dval;
	vsqu[ir] += wt*dval*dval;
      } else
	nang++;
    } /* i */
  } /* j */

  /* write output map(s), if needed */
  if (Qout) {
    dprintf(0,"Data min/max = %g %g\n",dmin,dmax);    
    MapMin(outptr) = dmin;
    MapMax(outptr) = dmax;
    write_image(outstr,outptr);
  }

  /* report on the rings */

  fsum = 0.0;
  nsum = 0;
  for (i=0; i<nring; i++) {
    if (wsum[i] == 0.0) continue;
    nsum++;
    r = 0.5*(rad[i] + rad[i+1]);
    dr = rad[i+1] - rad[i];
    area = PI*(sqr(rad[i+1]) - sqr(rad[i]));
    tmp = wsum[i]/pixe[i];
    ave = vsum[i]/wsum[i];
    rms = vsqu[i]/wsum[i]-ave*ave;
    if (rms <  0) rms=0.0;
    rms = sqrt(rms);
    fsum += rms*rms;

    if (Qtab) fprintf(tabstr,"%g %g %g %g %g %g ;; %g %g\n",
	   r,ave,rms,pixe[i],tmp,sqrt(fsum/nsum),    vsum[i],wsum[i]);
  }

  dprintf(0,"Nundf=%d/%d Nout=%d Nang=%d (sum=%d)\n",
	  nundf,nx*ny,nout,nang,nout+nundf+nang);
  dprintf(0,"Rmin/max = %g %g\n",rmin,rmax);

}

