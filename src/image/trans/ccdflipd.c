/* 
 * CCDFLIPD: patch up an image where neighbors are within +/- of each other, i.e. flip the data
 *           cloned off ccdfill :-)
 *
 *	14-mar-2022	written	for N253 outflow vs. disk - very quick and dirty      pjt
 *                      
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
	"in=???\n       Input image file",
	"out=???\n      Output image file",
	"delta=20\n     Allowed delta between neighoring pixels",
	"m=1\n          Minimum number of neighbors needed to flip",
	"iter=50\n      Number of iterations",
	"VERSION=0.3\n  14-mar-2022 PJT",
	NULL,
};

string usage = "patch up neighbors that are +/- of each other within a delta";




void nemo_main(void)
{
  stream   instr, outstr;
  int      nx, ny, nz;  
  int      nflip, mflip;
  int      i,j,k, di, dj;
  imageptr iptr=NULL;
  int      iter, niter = getiparam("iter");
  real     delta = getrparam("delta");
  real     m1, m2, dmin, dmax;
  int      n=1, m = getiparam("m");
  
  instr = stropen(getparam("in"), "r");
  read_image( instr, &iptr);
  nx = Nx(iptr);
  ny = Ny(iptr);
  nz = Nz(iptr);
  dmin = MapMax(iptr);
  dmax = MapMin(iptr);
  dprintf(0,"Old data min/max: %g %g\n",dmax,dmin);
  if (nz>1) warning("Patching first plane only");

  
  outstr = stropen(getparam("out"), "w");

  for (iter=0; iter < niter; iter++) {
    nflip = 0;
    for (j=0; j<ny; j++) {          
      for (i=0; i<nx; i++) {
	m1 = MapValue(iptr,i,j);
	dmin = MIN(dmin,m1);
	dmax = MAX(dmax,m1);
	if (m1 > delta) {

	  // @todo  first grab how many are changed in the neighbor list
	  //        require a minimum number 
	  if (m==1) {
	    for (dj=-n; dj<=n; dj++) {               // look at all neighbors
	      if (j+dj<0 || j+dj>=ny) continue;
	      for (di=-n; di<=n; di++) {
		if (i+di<0 || i+di>=nx) continue;
		m2 = MapValue(iptr,i+di,j+dj);
		if (ABS(m1+m2)  < delta) {
		  MapValue(iptr,i+di,j+dj) = -m2;    // m1 was well above 0 
		  nflip++;
		  m2 = -m2;
		  dmin = MIN(dmin,m2);
		  dmax = MAX(dmax,m2);
		}
	      } // di
	    } // dj
	  } else {
	    mflip=0;
	    for (dj=-n; dj<=n; dj++) {               // look at all neighbors
	      if (j+dj<0 || j+dj>=ny) continue;
	      for (di=-n; di<=n; di++) {
		if (i+di<0 || i+di>=nx) continue;
		m2 = MapValue(iptr,i+di,j+dj);
		if (ABS(m1+m2)  < delta)
		  mflip++;
	      } // di
	    } // dj
	    if (mflip >= m) {
	      for (dj=-n; dj<=n; dj++) {               // look at all neighbors
		if (j+dj<0 || j+dj>=ny) continue;
		for (di=-n; di<=n; di++) {
		  if (i+di<0 || i+di>=nx) continue;
		  m2 = MapValue(iptr,i+di,j+dj);
		  if (ABS(m1+m2)  < delta) {
		    MapValue(iptr,i+di,j+dj) = -m2;    // m1 was well above 0 
		    nflip++;
		    m2 = -m2;
		    dmin = MIN(dmin,m2);
		    dmax = MAX(dmax,m2);
		  }
		} // di
	      } // dj
	    }
	      
	  }
	}
      } // i
    } // j
    dprintf(0,"Iter %d: found %d values to patch\n",iter+1,nflip);
    if (nflip == 0) break;
  } // iter

  dprintf(0,"New data min/max: %g %g\n",dmin,dmax);  
  MapMin(iptr) = dmin;
  MapMax(iptr) = dmax;
  write_image(outstr, iptr);
}

// @todo   how is this possible?
// Old data min/max: -452.245 702.553
// New data min/max: -581.688 527.792
