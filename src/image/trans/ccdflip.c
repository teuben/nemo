/* 
 * CCDFLIP: flip an image
 *	quick and dirty: 13-jul-89
 *
 *	7-mar-92    gcc happy
 *	21-may-92   fixed wrong getparam parameter "x"
 *      22-feb-97   allow flipping x and y if the image is square,
 *		    also fixed an initialization problem!
 *      17-jun-2019 allow flipping in z
 *      15-sep-2023 optionally fix WCS also instead of flipping data only
 *                      
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
        "in=???\n       Input image file",
	"out=???\n      Output file",
	"flip=x\n       Flip in x,y,z  or allow xy for square images",
	"wcs=t\n        Also fix WCS?",
	"VERSION=2.0\n  15-sep-2023 PJT",
	NULL,
};

string usage = "flip an image along certain axes";

#define X 1
#define Y 2
#define Z 4
#define XY 3

#define SWAPR(a,b)	          \
{                                 \
    register real _tmpr = (a);    \
    a = b; b = _tmpr;             \
}

#define SWAPS(a,b)	          \
{                                 \
    register char * _tmps = (a);  \
    a = b; b = _tmps;             \
}


void nemo_main()
{
    stream  instr, outstr;
    int     nx, ny, nz;        /* size of scratch map */
    int     ix, iy, iz, flip;
    imageptr iptr=NULL;        /* pointer to image */
    real    tmp, zzz;
    string  flipmode;
    bool    Qwcs = getbparam("wcs");

    flipmode = getparam("flip");
    if (streq(flipmode,"x"))
        flip = X;
    else if (streq(flipmode,"y"))
        flip = Y;
    else if (streq(flipmode,"xy"))
        flip = XY;
    else if (streq(flipmode,"z"))
        flip = Z;      
    else
        warning("No flipping done");

    instr = stropen(getparam("in"), "r");
    outstr = stropen(getparam("out"), "w");

    read_image( instr, &iptr);
    if (Axis(iptr) > 0 && !Qwcs)
      warning("WCS is not changed, only values flipped");

    nx = Nx(iptr);	
    ny = Ny(iptr);
    nz = Nz(iptr);

    if(flip==X) {
      // @todo do for all Z
      for (iy=0; iy<ny; iy++) {		    /* flip in x */
        for (ix=0; ix<nx/2; ix++) {
            tmp = MapValue(iptr,ix,iy);
            zzz = MapValue(iptr,nx-ix-1,iy);
            dprintf(1,"%d %d: %f %f\n",ix,iy,tmp,zzz);
            MapValue(iptr,ix,iy) = MapValue(iptr,nx-ix-1,iy);
            MapValue(iptr,nx-ix-1,iy) = tmp;
        }
      }
      if (Qwcs) {
	dprintf(0,"Fixing WCS in X\n");
	Dx(iptr) = -Dx(iptr);
	Xref(iptr) = Nx(iptr)-1 - Xref(iptr);
      }
    } else if (flip==Y) {
      // @todo do for all Z      
      for (iy=0; iy<ny; iy++) {		    /* flip in y */
        for (ix=0; ix<nx/2; ix++) {
            tmp = MapValue(iptr,iy,ix);
            zzz = MapValue(iptr,iy,nx-ix-1);
            dprintf(1,"%d %d: %f %f\n",ix,iy,tmp,zzz);
            MapValue(iptr,iy,ix) = MapValue(iptr,iy,nx-ix-1);
            MapValue(iptr,iy,nx-ix-1) = tmp;
        }
      }
      if (Qwcs) {
	dprintf(0,"Fixing WCS in Y\n");
	Dy(iptr) = -Dy(iptr);
	Yref(iptr) = Ny(iptr)-1 - Yref(iptr);
      }
    } else if (flip==XY) {
      if (nx != ny) error("Cannot flip non-square images yet in XY");
      // @todo do for all Z      
      for (iy=0; iy<ny; iy++) {		    /* swap the x and y axes */
        for (ix=iy+1; ix<nx; ix++) {
            tmp = MapValue(iptr,ix,iy);
            dprintf(1,"%d %d: %f \n",ix,iy,tmp);
            MapValue(iptr,ix,iy) = MapValue(iptr,iy,ix);
            MapValue(iptr,iy,ix) = tmp;
        }
      }
      SWAPR(Xmin(iptr),  Ymin(iptr));
      SWAPR(Dx(iptr),    Dy(iptr));
      SWAPS(Namex(iptr), Namey(iptr));
    } else if (flip==Z) {
      for (ix=0; ix<nx; ix++) {
	for (iy=0; iy<ny; iy++) {	
	  for (iz=0; iz<nz/2; iz++) {	    /* flip in z */
	    tmp = CubeValue(iptr,ix,iy,iz);
            zzz = CubeValue(iptr,ix,iy,nz-iz-1);
            dprintf(1,"%d %d: %f %f\n",ix,iy,iz,tmp,zzz);
            CubeValue(iptr,ix,iy,iz) = CubeValue(iptr,ix,iy,nz-iz-1);
            CubeValue(iptr,ix,iy,nz-iz-1) = tmp;
	  }
        }
      }
      if (Qwcs) {
	dprintf(0,"Fixing WCS in Z\n");
	Dz(iptr) = -Dz(iptr);
	Zref(iptr) = Nz(iptr)-1 - Zref(iptr);
      }
    }
    write_image(outstr, iptr);
}
