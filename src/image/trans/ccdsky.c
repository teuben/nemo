/* 
 * CCDSKY: lazy sky scaling calculator to help converting theory data
 *         into degrees and m/s for FITS
 *         Can optionally convert an image. 
 *
 *      17-aug-2012     created        Peter Teuben
 *      22-aug-2012     added sdv=     PJT
 *      28-aug-2012     implemented iscale=
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

#include <mks.h>

string defv[] = {
  "in=\n          Optional input image file",
  "out=\n         Output image file",
  "d=1,pc\n       Distance to object, optional unit [pc]",
  "r=1,AU\n       Length scale of object, optional unit [AU]",
  "v=1,km/s\n     Velocity scale of object, optional unit [km/s]",
  "sdv=1\n        Integrated Flux (must be in Jy.km/s)",
  "scale=1\n      Scale image values [not implemented]",
  "VERSION=1.2\n  28-aug-2012 PJT",
  NULL,
};

string usage = "lazy sky scaling calculator";

string cvsid = "$Id$";


#define CVI(x,y,z)  CubeValue(iptr,x,y,z)


#define MAXU 64


/*   
 * See also:  http://www.cv.nrao.edu/course/astr534/HILine.html
 */

static real HI_factor = 2.35e5;

/* 
 * Xco=2e20 cm-2/K km/s, and alpha_co=4.3 Msun/K km/s 
 * also included 1.36  factor due to Helium contribution to the mass. 
 * Use Xco=14e20 for 13CO, about 7 times that of 12CO
 *
 * Older material is also: http://ned.ipac.caltech.edu/level5/March09/Solomon/Solomon2.html
 */

static real CO_factor = 1.05e4;





/*
 *     get_nu :   parse a   NUMBER,UNIT    string
 *     input:   kv         e.g.   "1,pc"
 *              defunit    e.f.   "pc"
 *     output:  value 
 *              unit
 *
 */

void get_nu(string kv,real *value,string unit, string defunit)
{
  string *sp = burststring(kv,",");
  int nsp = xstrlen(sp,sizeof(string)) - 1;
  int nv;

  *value = 1.0;
  strcpy(unit,defunit);
  if (nsp == 0) return;

  nv = nemoinpr(sp[0],value,1);
  if (nv != 1) error("error parsing %s",sp[0]);
  if (nsp == 1) return;

  strcpy(unit,sp[1]);
}

/*
 * efactor:  scaling factor between two common units
 *           u1:    first unit
 *           u2:    second unit
 *           returned value:   how many u1's in a u2  (i.e.   u2/u1)      
 *
 *           @todo need a more general unit + prefix conversion, this
 *                 was just a quick hack
 */

real efactor(string u1, string u2)
{
  real s = -1.0;
  if (streq(u1,"AU")) {
    if (streq(u2,"AU")) 
      s = 1.0;
    else if (streq(u2,"pc")) 
      s = PC/AU;
    else if (streq(u2,"km")) 
      s = 1000.0/AU;
    else if (streq(u2,"Kpc") || streq(u2,"kpc"))
      s = 1e3*PC/AU;
    else if (streq(u2,"Mpc") || streq(u2,"mpc"))
      s = 1e6*PC/AU;
    else
      warning("Comparison_1 unit %s not understood for %s",u2,u1);
  } else if (streq(u1,"pc")) {
    if (streq(u2,"pc")) 
      s = 1.0;
    else if (streq(u2,"Kpc") || streq(u2,"kpc"))
      s = 1e3;
    else if (streq(u2,"Mpc") || streq(u2,"mpc"))
      s = 1e6;
    else if (streq(u2,"Gpc") || streq(u2,"Gpc"))
      s = 1e9;
    else if (streq(u2,"AU"))
      s = AU/PC;
    else
      warning("Comparison_2 unit %s not understood for %s",u2,u1);
  } else if (streq(u1,"Mpc")) {
    if (streq(u2,"pc")) 
      s = 1e-6;
    else if (streq(u2,"Kpc") || streq(u2,"kpc"))
      s = 1e-3;
    else if (streq(u2,"Mpc") || streq(u2,"mpc"))
      s = 1;
    else if (streq(u2,"Gpc") || streq(u2,"gpc"))
      s = 1e3;
    else
      warning("Comparison_3 unit %s not understood for %s",u2,u1);
  } else if (streq(u1,"m/s")) {
    if (streq(u2,"m/s")) 
      s = 1;
    else if  (streq(u2,"km/s")) 
      s = 1e3;
    else
      warning("Comparison_4 unit %s not understood for %s",u2,u1);      
  } else
    warning ("Comparison_5 unit %s not understood",u1);
  dprintf(1,"u1=%s u2=%s  s=%g\n",u1,u2,s);
  return s;
}


void nemo_main()
{
    stream  instr, outstr;
    int     ix, iy, iz, nx, ny, nz;
    imageptr iptr=NULL;
    real d,r,v,  rscale, vscale, iscale = getrparam("scale");
    real mass, sdv_scale,  sdv = getrparam("sdv");
    char ds[MAXU], rs[MAXU], vs[MAXU];

    /*  get the distance and scales      */

    get_nu(getparam("d"),&d,ds,"pc");
    get_nu(getparam("r"),&r,rs,"AU");
    get_nu(getparam("v"),&v,vs,"km/s");

    printf("d=%g %s\n",d,ds);
    printf("r=%g %s\n",r,rs);
    printf("v=%g %s\n",v,vs);
    printf("SdV=%g Jy.km/s\n",sdv);


    /* angles: convert to deg for FITS */
    rscale = (r * AU) / ( d * PC) * 180.0/PI;   
    rscale *= efactor("AU",rs);
    rscale /= efactor("pc",ds);
    printf("rscale=%g  (%g arcsec)\n",rscale,rscale*3600.0);

    /* velocities:  convert to m/s for FITS */
    vscale = v;
    vscale *= efactor("m/s",vs);
    printf("vscale=%g\n",vscale);

    /* mapvalues */
    printf("iscale=%g\n",iscale);

    /* sdv calculation */
    sdv_scale = efactor("Mpc",ds);
    mass = HI_factor * sqr( d * sdv_scale ) * sdv;
    printf("Mass(HI) = %g  \n",mass);
    mass = CO_factor * sqr( d * sdv_scale ) * sdv;
    printf("Mass(H2) = %g  (alpha=4.3; includes 1.36 He contribution)\n",mass);
    

    if (hasvalue("in") && hasvalue("out")) {      /* patch image if needed */
      instr = stropen(getparam("in"), "r");
      read_image( instr, &iptr);
      Xmin(iptr) *= rscale;
      Ymin(iptr) *= rscale;
      Zmin(iptr) *= vscale;
      Dx(iptr) *= rscale;
      Dy(iptr) *= rscale;
      Dz(iptr) *= vscale;
      if (iscale != 1.0) {
	nx = Nx(iptr);
	ny = Ny(iptr);
	nz = Nz(iptr);
	for (iz=0; iz<nz; iz++)
	  for (iy=0; iy<ny; iy++)
	    for (ix=0; ix<nx; ix++)
	      CubeValue(iptr,ix,iy,iz) *= iscale;
      }

      outstr = stropen(getparam("out"), "w");
      write_image(outstr, iptr);
    }
}

