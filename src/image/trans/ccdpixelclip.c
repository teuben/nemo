/* 
 * CCDPIXELCLIP: clip an image, based on Joe DePasquale's PixelClip.js v1.1
 *
 *     2022 Joseph DePasquale - STScI
 *     Adapted from Copyright (C) 2019 Gerrit Barrere
 *
 *     This script scans through an image and replaces pixels with any 
 *     component above a specified threshold with the mean of its under-threshold 
 *     nearest neighbors. This can be used to replace blown-out 
 *     star cores in monochromatic telescope images. 
 *     The blown-out core is replaced with pixel values just outside the 
 *     blown-out region. This can be done early in the linear (pre-stretched) 
 *     state.  It restores dynamic range to that of the camera and recovers 
 *     the blown-out star cores. 
 *
 *	24-aug-2022    drafted for NEMO     - Peter Teuben
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
  "clip=0\n       clip threshold below which data are patched",
  "size=1\n       size around clip pixel to average around (-size..+size)",
  "matchclip=t\n  Also fix data where it matches the clip value?",
  "jwst=t\n       JWST mode? (only true is accepted now)",
  "VERSION=0.3\n  25-aug-2022 PJT",
  NULL,
};

string usage = "pixelclip an image - meant for filling in JWST clipped stars";


void nemo_main()
{
    stream   instr, outstr;
    int      nx, ny, nz; 
    int      ix, iy, iz, ix1, iy1;
    imageptr iptr=NULL; 
    real     clip = getrparam("clip");
    int      size = getiparam("size");
    int      countclip;
    real     sum0, sum1;
    bool     Qmatch = getbparam("matchclip");
    bool     Qjwst =  getbparam("jwst");

    if (!Qjwst) error("jwst=f mode not implemented yet");

    instr = stropen(getparam("in"), "r");
    outstr = stropen(getparam("out"), "w");

    read_image( instr, &iptr);
    nx = Nx(iptr);	
    ny = Ny(iptr);
    nz = Nz(iptr);

    countclip = 0;
    for (iz=0; iz<nz; iz++) {
      for (iy=0; iy<ny; iy++) {
	if (iy<size || iy>=ny-size) continue;
        for (ix=0; ix<nx; ix++) {
	  if (ix<size || ix>=nx-size) continue;
	  // @todo  implemented JWST mode
	  if ((Qmatch && CubeValue(iptr,ix,iy,iz) <= clip) || (CubeValue(iptr,ix,iy,iz) < clip)) {
	    sum0 = sum1 = 0.0;
	    for (iy1=-size; iy1<=size; iy1++) {
	      for (ix1=-size; ix1<=size; ix1++) {
		if ((CubeValue(iptr,ix+ix1,iy+iy1,iz) > clip) || (CubeValue(iptr,ix+ix1,iy+iy1,iz) >= clip)) {
		  sum0++;
		  sum1 += CubeValue(iptr,ix+ix1,iy+iy1,iz);
		} 
	      }
	    }
	    if (sum0 > 0) {
	      countclip++;
	      CubeValue(iptr,ix,iy,iz) = sum1/sum0;
	    } else
	      CubeValue(iptr,ix,iy,iz) = clip;
	    dprintf(1,"clip %d   %d %d %d -> %g %g\n",
		    countclip,ix,iy,iz,CubeValue(iptr,ix,iy,iz),sum0);
	  }
	}
      }
    }
    dprintf(0,"Clipped %d/%d %g%% values %s %g\n",
	    countclip,nx*ny*nz, (100.0*countclip)/(nx*ny*nz),
	    Qjwst ? (Qmatch ? "equal or below" : "below") :
	            (Qmatch ? "equal or above" : "above"),
	    clip);
    write_image(outstr, iptr);
}

