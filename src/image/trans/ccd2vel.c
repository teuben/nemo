/* 
 * CCD2VEL: convert a map of frequencies or wavelengths to doppler velocity,
 *          in any frame (radio, optical, relativistic)
 *
 *      27-apr-2017     created        Peter Teuben
 *
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

#include <mks.h>

string defv[] = {
  "in=\n            Optional input image file",
  "out=\n           Output image file",
  "from=rel\n       Input reference frame of reference, if converting from velocities (relativistic, radio, optical)",
  "to=rel\n         Output reference frame of reference (relativistic, radio, optical)",
  "restwave=\n      Input rest wavelength, in the units of the map, only if needed",
  "restfreq=\n      Input rest frequency, in the units of the map, only if needed",
  "clip=\n          Clipping of output (none, or two values needed)",
  "units=km/s\n     Desired output units (km/s, m/s)",
  "VERSION=0.3\n    28-apr-2017 PJT",
  NULL,
};

string usage = "convert map of wave/freq to doppler velocities";

string cvsid = "$Id$";

static real c = c_MKS / 1000.0;    /* use km/s for now */

#define CVI(x,y,z)  CubeValue(iptr,x,y,z)


#define ONE2ONE  0x00

#define WAV2REL  0x01
real wav2rel(real w,real w0)
{
  return c * (w*w-w0*w0) / (w*w+w0*w0);
}

#define FRQ2REL  0x02
real frq2rel(real f,real f0)
{
  return c * (f0*f0-f*f) / (f0*f0+f*f);
}

#define OPT2RAD  0x04
real opt2rad(real v)
{
  return v / (1 + v/c);
}

#define OPT2REL  0x08
real opt2rel(real v)
{
  real q = 1 + v/c;
  return c * (q*q-1)/(q*q+1);
}

#define RAD2OPT  0x10
real rad2opt(real v)
{
  return v / (1 - v/c);
}

#define RAD2REL  0x20
real rad2rel(real v)
{
  real q = 1 - v/c;
  return c * (1-q*q)/(1+q*q);
}

#define REL2OPT  0x40
real rel2opt(real v)
{
  real z = v/c;
  return c*(sqrt((1+z)/(1-z)) - 1);
}

#define REL2RAD 0x80
real rel2rad(real v)
{
  real z = v/c;
  return c*(1 - sqrt((1-z)/(1+z)));
}

#define WAV2OPT  0x0100
real wav2opt(real w,real w0)
{
  return c * (w/w0 - 1.0);
}

#define FRQ2RAD  0x0200
real frq2rad(real f,real f0)
{
  return c * (1.0 - f/f0);
}




void nemo_main()
{
    stream  instr, outstr;
    int     ix, iy, iz, nx, ny;
    int     mode1, mode2;
    imageptr iptr=NULL;
    string to    = getparam("to");
    string from  = getparam("from");
    bool Qclip, Qfirst, Qlinear = FALSE;
    real isca = 1.0, ioff = 0.0;
    real osca = 1.0, ooff = 0.0;
    real zsca = 1.0, zoff = 0.0;
    real vel, vmin, vmax, clip[2];
    real restwave = -1;
    real restfreq = -1;

    if (hasvalue("restwave")) restwave = getrparam("restwave");
    if (hasvalue("restfreq")) restfreq = getrparam("restfreq");
    Qclip = hasvalue("clip");
    if (Qclip) {
      if (nemoinpr(getparam("clip"),clip,2) != 2) error("Two values needed for clip=");
      dprintf(0,"Clipping on %g %g\n",clip[0], clip[1]);
    }

    if (restwave > 0)
      mode1 = WAV2REL;
    else if (restfreq > 0)
      mode1 = FRQ2REL;
    else if (streq(from,"opt"))
      mode1 = OPT2REL;
    else if (streq(from,"rad"))
      mode1 = RAD2REL;
    else if (streq(from,"rel"))
      mode1 = ONE2ONE;
    else
      error("Unknown from= frame");

    if (streq(to,"rad"))
      mode2 = REL2RAD;
    else if (streq(to,"opt"))
      mode2 = REL2OPT;
    else if (streq(to,"rel"))
      mode2 = ONE2ONE;
    else
      error("Unknown to= frame");

    dprintf(0,"MODES: 0x%x  and 0x%x\n",mode1,mode2);
    if (mode1 == OPT2REL && mode2 == REL2OPT) { mode1 = mode2 = ONE2ONE;}
    if (mode1 == RAD2REL && mode2 == REL2RAD) { mode1 = mode2 = ONE2ONE;}
    if (mode1 == RAD2REL && mode2 == REL2OPT) { mode1 = RAD2OPT; mode2 = ONE2ONE;}
    if (mode1 == OPT2REL && mode2 == REL2RAD) { mode1 = OPT2RAD; mode2 = ONE2ONE;}
    if (mode1 == WAV2REL && mode2 == REL2OPT) { mode1 = WAV2OPT; mode2 = ONE2ONE;}
    if (mode1 == FRQ2REL && mode2 == REL2RAD) { mode1 = FRQ2RAD; mode2 = ONE2ONE;}        
    dprintf(0,"MODES: 0x%x  and 0x%x\n",mode1,mode2);
    
    if (hasvalue("in") && hasvalue("out")) {    
      instr = stropen(getparam("in"), "r");
      read_image( instr, &iptr);
      nx = Nx(iptr);
      ny = Ny(iptr);
      Qfirst = TRUE;
      for (iy=0; iy<ny; iy++)
	for (ix=0; ix<nx; ix++) {
	  vel = MapValue(iptr,ix,iy);
	  if (Qlinear)
	    vel = zsca * vel + zoff;	  /* linear velocity conversion (radio,optical) */
	  else {
	    if (mode1 == WAV2REL)
	      vel = wav2rel(vel, restwave);
	    else if (mode1 == WAV2OPT)
	      vel = wav2opt(vel, restwave);
	    else if (mode1 == FRQ2REL)
	      vel = frq2rel(vel, restfreq);
	    else if (mode1 == FRQ2RAD)
	      vel = frq2rad(vel, restfreq);
	    else if (mode1 == OPT2REL)
	      vel = opt2rel(vel);
	    else if (mode1 == RAD2REL)
	      vel = rad2rel(vel);
	    else if (mode1 == OPT2RAD)
	      vel = opt2rad(vel);
	    else if (mode1 == RAD2OPT)
	      vel = rad2opt(vel);
 
	    if (mode2 == REL2OPT)
	      vel = rel2opt(vel);
	    else if (mode2 == REL2RAD)
	      vel = rel2rad(vel);
	  }

	  if (Qclip)
	    if (vel < clip[0] || vel > clip[1]) vel = 0.0;

	  if (Qfirst) {
	    if (vel != 0.0) {
	      vmin = vmax = vel;
	      Qfirst = FALSE;
	    }
	  }
	  if (vel != 0.0) {
	    vmin = MIN(vel,vmin);
	    vmax = MAX(vel,vmax);
	  }
	  MapValue(iptr,ix,iy) = vel;
	}
      MapMin(iptr) = vmin;
      MapMax(iptr) = vmax;
      dprintf(0,"MinMax: %g %g\n",vmin,vmax);
      outstr = stropen(getparam("out"), "w");
      write_image(outstr, iptr);
    }
}
