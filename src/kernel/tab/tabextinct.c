/* TABEXTINCT:   extinct a spectrum
 *
 *
 *      13-may-05    first version, cloned off tabfilter
 *      16-may-05    added xmin,xmax
 *
 *  A note on intervals:  we need 'model' to encompass 'in' completely, 
 *                        but with xmin= (xQmin) and xmax= (xQmax) can now cheat
 *
 *     xmin                                                    xmax
 *       [-------------------- model ----------------------------]
 *
 *                 [----------- in ---------------]
 *               umin                           umax
 */

#include <stdinc.h> 
#include <getparam.h>
#include <spline.h>
#include <extstring.h>
#include <table.h>

string defv[] = {
  "model=???\n      Extinction model (table with C(wave) in col 1/2)",
  "in=???\n         Input spectrum table file",
  "xcol=1\n         Wavelength column for spectrum",
  "ycol=2\n         Flux column for spectrum",
  "xscale=1\n       Scale input wavelength to match the model",
  "yscale=1\n       Scale spectrum for your preferred units",
  "xmin=\n          Ignore points below this value",
  "xmax=\n          Ignore points above this value",
  "Av=1\n           Av to apply extinction curve with",
  "extinct=t\n      Extinction law, or some other linear law",
  "VERSION=0.4\n    16-may-05 PJT",
  NULL,

};

string usage="extinct a spectrum";

string cvsid="$Id$";


#define MAXDATA  16384

extern string *burststring(string, string);
extern void freestrings(string *);

/* we might want to re-use the filter (model) lookup technique used in tabfilter ?? */

string filtername(string shortname)
{
  int nsp;
  static char fullname[256];
  char line[MAX_LINELEN];
  string *sp, fpath;
  stream fstr;

  if (nemo_file_size(shortname) > 0) return shortname;

  fpath = getenv("NEMODAT");
  if (fpath == 0) error("NEMODAT does not exist");

  sprintf(fullname,"%s/filter/Fnames",fpath);
  dprintf(1,"Alias table %s\n",fullname);
  fstr = stropen(fullname,"r");
  while (fgets(line,MAX_LINELEN,fstr)) {
    if (line[0] == '#') continue;
    sp = burststring(line," \n");
    nsp = xstrlen(sp,sizeof(sp))-1;
    if (nsp > 1) {
      if (streq(shortname,sp[0])) {
	sprintf(fullname,"%s/filter/%s",fpath,sp[1]);
	dprintf(1,"Matching %s\n",fullname);
	freestrings(sp);
	return fullname;
      }
    }
    freestrings(sp);
  }
  return shortname;
}


void nemo_main()
{
  int colnr[2];
  real *coldat[2], *xdat, *ydat, xmin, xmax, ymin, ymax;
  real *udat, *vdat, umin, umax, vmin, vmax;
  real Av, C, xscale, yscale;
  stream instr;
  int i, n, ns, nmax;
  real *sdat, xQmin, xQmax;
  string spectrum = getparam("in");
  string model = getparam("model");
  bool Qextinct = getbparam("extinct");
  bool Qmin = hasvalue("xmin");
  bool Qmax = hasvalue("xmax");
  Av = getdparam("Av");
  xscale = getdparam("xscale");
  yscale = getdparam("yscale");
  
  nmax = nemo_file_lines(model,MAXLINES);
  xdat = coldat[0] = (real *) allocate(nmax*sizeof(real));
  ydat = coldat[1] = (real *) allocate(nmax*sizeof(real));
  colnr[0] = 1;  /* wavelenght in angstrom */
  colnr[1] = 2;  /* C */
  instr = stropen(model,"r");
  n = get_atable(instr,2,colnr,coldat,nmax);
  strclose(instr);
  
  for(i=0; i<n; i++) {
    dprintf(2,"%g %g\n",xdat[i],ydat[i]);
    if (i==0) {
      xmin = xmax = xdat[0];
      ymin = ymax = ydat[0];
    } else {
      if (xdat[i] <= xdat[i-1]) 
	error("Model %s must be sorted in wavelength",model);
      xmax = MAX(xmax,xdat[i]);
      ymax = MAX(ymax,ydat[i]);
      xmin = MIN(xmin,xdat[i]);
      ymin = MIN(ymin,ydat[i]);
    }
  }
  dprintf(1,"Model wavelength range: %g : %g\n",xmin,xmax);
  dprintf(1,"Model extinction range: %g : %g\n",ymin,ymax);
  
  /* setup a spline interpolation table into the extinction curve */
  sdat = (real *) allocate(sizeof(real)*n*3);
  spline(sdat,xdat,ydat,n);

  if (Qmin) xQmin = getdparam("xmin");
  if (Qmax) xQmax = getdparam("xmax");

  spectrum = getparam("in");
  nmax = nemo_file_lines(spectrum,MAXLINES);
  udat = coldat[0] = (real *) allocate(nmax*sizeof(real));
  vdat = coldat[1] = (real *) allocate(nmax*sizeof(real));
  colnr[0] = getiparam("xcol");
  colnr[1] = getiparam("ycol");
  instr = stropen(spectrum,"r");
  ns = get_atable(instr,2,colnr,coldat,nmax);
  strclose(instr);

  for(i=0; i<ns; i++) {
    dprintf(2,"%g %g\n",udat[i],vdat[i]);
    udat[i] *= xscale;
    vdat[i] *= yscale;
    if (i==0) {
      umin = umax = udat[0];
      vmin = vmax = vdat[0];
    } else {
      if (udat[i] <= udat[i-1])
	error("Spectrum %s must be sorted in wavelength",spectrum);
      umax = MAX(umax,udat[i]);
      vmax = MAX(vmax,vdat[i]);
      umin = MIN(umin,udat[i]);
      vmin = MIN(vmin,vdat[i]);
    }
  }
  dprintf(1,"Spectrum wavelength range: %g : %g\n",umin,umax);
  dprintf(1,"Spectrum response range: %g : %g\n",vmin,vmax);

  if (xmin > umin || xmax < umax)
    error("Spectrum in not embedded inside Extinction curve");
  if (Qmin && xQmin < umin) warning("xmin=%g less than minimum (%g)",xQmin,umin);
  if (Qmax && xQmax < umax) warning("xmax=%g greater than maximum (%g)",xQmax,umax);
  if (Qmin && xQmin > umax) error("xmin=%g greater than maximum (%g)",xQmin,umax);
  if (Qmax && xQmax > umin) error("xmax=%g less than minimum (%g)",xQmax,umin);

  for (i=0; i<ns; i++) {           /* loop over the spectrum */
    if (Qmin && udat[i] < xQmin) continue;
    if (Qmax && udat[i] > xQmax) continue;
    C = seval(udat[i],xdat,ydat,sdat,n);    /* model */
    dprintf(1,"%g %g %g %g\n",udat[i],vdat[i],C,vdat[i]*C);
    if (Qextinct)
      vdat[i] *= pow(10.0,-0.4*Av*C);       /* extinction law */
    else
      vdat[i] *= Av*C;                      /* simple linear law */
    printf("%g %g\n",udat[i],vdat[i]);
  }
}

