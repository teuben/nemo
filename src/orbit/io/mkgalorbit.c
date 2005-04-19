/*
 * MKORBIT:	generate an orbit from observational input variables
 *  
 *  http://idlastro.gsfc.nasa.gov/ftp/pro/astro/gal_uvw.pro
 *
 *      18-apr-05   created        --   peter teuben
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <orbit.h>

string defv[] = {
    "out=???\n		  output filename (an orbit)",
    "lon=\n               Longitude (in degrees)",
    "lat=\n               Latitude (in degrees)",
    "pmlon=\n             Proper Motion in Longitude (mas/yr)",
    "pmlon=\n             Proper Motion in Latitude (mas/yr)",
    "vrad=\n              Radial velocity (km/s)",
    "dist=\n              Distance (in kpc)",
    "R0=8.0\n             LSR distance (in kpc)",
    "V0=220.0\n           LSR velocities (in km/s)",
    "lsr=f\n              Use LSR correction of sun going (-9,12,7)",
    "coordsys=gal\n       lon/lat coordinate system: equ=RA/DEC gal=GLON/GLAT",

    "time=0.0\n           time",
    "potname=zero\n	  optional potential(5NEMO)",
    "potpars=\n	          .. with optional parameters",
    "potfile=\n		  .. and optional datafile name",
    "headline=\n          random verbiage",
    "VERSION=0.1\n        18-apr-05 PJT",
    NULL,
};

string usage = "Make a galactic orbit with from given initial conditions";

string cvsid="$Id$";


string	infile,outfile;			/* file names */
stream  instr,outstr;			/* file streams */

orbitptr optr;
a_potential p;

double x,y,z,u,v,w;
double etot, lz;
double tnow, omega;
double R0, V0;
double lon, lat, pmlon, pmlat, dist, vrad;

int  Dpos, Dvel;
bool Qgal;

string coordsys;

double k = 4.74047;     /*  conversion factor for 1 AU/yr in km/s */

#define D2R   PI/180.0


void setparams()
{
  potproc_double pot;
  double pos[3],acc[3],epot;
  int ndim=3;
  int signcount;
  int lzsign;
  
  p.name = getparam("potname");
  p.pars = getparam("potpars");
  p.file = getparam("potfile");

  lon = getdparam("lon") * D2R;
  lat = getdparam("lat") * D2R;

  pmlon = getdparam("lon") * k;
  pmlat = getdparam("lat") * k;

  R0 = getdparam("R0");
  V0 = getdparam("V0");

  dist = getdparam("dist");
  vrad = getdparam("vrad");

  coordsys = getparam("coordsys");
  if (streq(coordsys,"gal"))
    Qgal = TRUE;
  else
    error("Only galactic coordinate system supported now");
    
  tnow = getdparam("time");
  outfile = getparam("out");
  
  pot = get_potential_double(p.name, p.pars, p.file);
  if (pot==NULL) 
    error("potential %s cannot be loaded",p.name);
  omega = get_pattern();
  dprintf(0,"Using pattern speed = %g\n",omega);

  /* sun is at (-R0,0,0) with LSR moving at (0,V0,0) */

  x = -R0 + dist * cos(lat) * cos(lon);
  y =       dist * cos(lat) * sin(lon);
  z =       dist * sin(lat);

  /* U,V,W are -vx,vy,vz */

  u =      pmlon * cos(lat)  + pmlat * sin(lat)    - vrad * cos(lat) * cos(lon);
  v = V0 + pmlon * cos(lat)  + pmlat * cos(lat)    + vrad * cos(lat) * sin(lon);
  w =                                              + vrad * sin(lat);

  pos[0] = x;
  pos[1] = y;
  pos[2] = z;
  (*pot)(&ndim,pos,acc,&epot,&tnow);
  epot -= 0.5*omega*omega*(x*x+y*y);
  etot = epot + 0.5*(u*u+v*v+w*w) - 0.5*omega*omega*(x*x+y*y);
  
  if (hasvalue("headline")) set_headline(getparam("headline"));
}


void nemo_main ()
{
  warning("new program, hasn't been tested out well");
  setparams();

  optr = NULL;				/* make an orbit */
  allocate_orbit (&optr, 3, 1);
  Masso(optr)  = 1.0;                     /* and set Mass */
  Torb(optr,0) = tnow;
  Xorb(optr,0) = x;			/* .. positions */
  Yorb(optr,0) = y;
  Zorb(optr,0) = z;
  Uorb(optr,0) = -u;			/* .. velocities : !!! notice vx = -u !!!!! */
  Vorb(optr,0) =  v;
  Worb(optr,0) =  w;
  I1(optr) = etot;			/*  energy (zero if not used) */
  I2(optr) = lz;                          /* angular momentum */
  
  dprintf(0,"pos: %f %f %f  \nvel: %f %f %f  \netot: %f\nlz=%f\n",
	  x,y,z,u,v,w,etot,lz);
  
  outstr = stropen (outfile,"w");		/* write to file */
  put_history(outstr);
  PotName(optr) = p.name;
  PotPars(optr) = p.pars;
  PotFile(optr) = p.file;
  write_orbit(outstr,optr);
  strclose(outstr);
}

