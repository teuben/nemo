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
  "frame=rel\n      Output frame of reference (relativistic, radio, optical)",
  "restwave=\n      Input rest wavelength, in the units of the map, only if needed",
  "restfreq=\n      Input rest frequency, in the units of the map, only if needed",
  "from=rel\n       Input frame of reference, if converting from velocities (relativistic, radio, optical)",
  "offset=0.0\n     Offset to subtract from map to be applied before conversion [not impl]",
  "scale=1.0\n      Scale applied to map after conversion [not impl]",
  "clip=\n          Clipping of output (none, or two values needed)",
  "units=km/s\n     Desired output units (km/s, m/s)",
  "VERSION=0.1\n    27-apr-2017 PJT",
  NULL,
};

string usage = "convert map of wave/freq to doppler velocities";

string cvsid = "$Id$";

static real c = c_MKS / 1000.0;    /* use km/s for now */

#define CVI(x,y,z)  CubeValue(iptr,x,y,z)


real rel_wave(real w,real w0)
{
  return c * (w*w-w0*w0) / (w*w+w0*w0);
}

real rel_freq(real f,real f0)
{
  return c * (f0*f0-f*f) / (f0*f0+f*f);
}


real opt2rad(real v)
{
  return v / (1 + v/c);
}

real opt2rel(real v)
{
  real q = 1 + v/c
  return (q*q-1)/(q*q+1);
}

real rad2opt(real v)
{
  return v / (1 - v/c);
}

real rad2rel(real v)
{
  real q = 1 - v/c
  return (1-q*q)/(1+q*q);
}

real rel2opt(real v)
{
  return -1.0;
}

real rel2rad(real v)
{
  return -1.0;
}





void nemo_main()
{
    stream  instr, outstr;
    int     ix, iy, iz, nx, ny;
    imageptr iptr=NULL;
    string frame = getparam("frame");
    bool Qclip, Qlinear, Qwave, Qfreq;
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
	
    Qlinear = Qwave = Qfreq = FALSE;
    if (streq(frame,"rel")) {
      dprintf(0,"Relativistic\n");
      if (restwave > 0) Qwave = TRUE;
      if (restfreq > 0) Qfreq = TRUE;
    } else
      Qlinear = TRUE;
    
    if (restwave > 0) {
      zsca = c / restwave;
      zoff = -c;
      dprintf(0,"Wave: %g : %g %g\n",restwave, zoff,zsca);
    } else if (restfreq > 0) {
      zsca = -c / restfreq;
      zoff = c;
      dprintf(0,"Freq: %g : %g %g\n",restfreq, zoff,zsca);
    } else
      warning("Did not supply restwave or restfreq");
      

    if (hasvalue("in") && hasvalue("out")) {      /* patch image if needed */
      instr = stropen(getparam("in"), "r");
      read_image( instr, &iptr);
      nx = Nx(iptr);
      ny = Ny(iptr);
      for (iy=0; iy<ny; iy++)
	for (ix=0; ix<nx; ix++) {
	  vel = MapValue(iptr,ix,iy);
	  vel = (vel + ioff) * isca;      /* before */
	  if (Qlinear)
	    vel = zsca * vel + zoff;	  /* linear velocity conversion (radio,optical) */
	  else if (Qwave)
	    vel = rel_wave(vel, restwave);
	  else if (Qfreq)
	    vel = rel_freq(vel, restfreq);
	  vel = (vel + ooff) * osca;      /* after */
	  if (Qclip)
	    if (vel < clip[0] || vel > clip[1]) vel = 0.0;
	  if (iy==0 && ix==0) {
	    vmin = vmax = vel;
	  } else {
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
