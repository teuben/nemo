/* 
 * CCDSKY: lazy sky scaling calculator to help converting theory data
 *         into degrees and m/s for FITS
 *         Can optionally convert an image. 
 *
 *      17-aug-2012     PJT     created
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
  "d=1,pc\n       Distance to object",
  "r=1,AU\n       Length scale of object",
  "v=1,km/s\n     Velocity scale of object",
  "scale=1        Scale image values [not implemented]",
  "VERSION=1.0\n  17-aug-2012 PJT",
  NULL,
};

string usage = "lazy sky scaling calculator";

string cvsid = "$Id$";


#define CVI(x,y,z)  CubeValue(iptr,x,y,z)


#define MAXU 64



/*
 *     get_nu :   parse a   NUMBER,UNIT    string
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
      warning("Comparison unit %s not understood for %s",u2,u1);
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
      warning("Comparison unit %s not understood for %s",u2,u1);
  } else if (streq(u1,"m/s")) {
    if (streq(u2,"m/s")) 
      s = 1;
    else if  (streq(u2,"km/s")) 
      s = 1e3;
    else
      warning("Comparison unit %s not understood for %s",u2,u1);      
  } else
    warning ("Comparison unit %s not understood",u1);
  dprintf(1,"u1=%s u2=%s  s=%g\n",u1,u2,s);
  return s;
}

void nemo_main()
{
    stream  instr, outstr;
    int     nx, ny, nz;
    imageptr iptr=NULL;
    real d,r,v,  rscale, vscale, iscale = getrparam("scale");
    char ds[MAXU], rs[MAXU], vs[MAXU];

    /*  get the distance and scales      */

    get_nu(getparam("d"),&d,ds,"pc");
    get_nu(getparam("r"),&r,rs,"AU");
    get_nu(getparam("v"),&v,vs,"km/s");

    printf("d=%g %s\n",d,ds);
    printf("r=%g %s\n",r,rs);
    printf("v=%g %s\n",v,vs);


    /* angles: convert to deg for FITS */
    rscale = (r * AU) / ( d * PC) * 180.0/PI;   
    rscale *= efactor("AU",rs);
    rscale /= efactor("pc",ds);
    printf("rscale=%g  (%g arcsec)\n",rscale,rscale*3600.0);

    /* velocities:  convert to m/s for FITS */
    vscale = v;
    vscale *= efactor("m/s",vs);
    printf("vscale=%g\n",vscale);
    

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
	warning("scale=%g not implemented yet",iscale);
      }

      outstr = stropen(getparam("out"), "w");
      write_image(outstr, iptr);
    }
}

