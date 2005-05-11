/* 
 * CCDINTPOL:   interpolate an image (like ccdfill) to patch up missing
 *              data
 *
 *	 8-jul-03       written, cloned off ccdfill,  Q&D again  PJT
 *                      only useable for regular good data rows/cols
 *                      use for the DH barfitting project
 *
 * some refs:
 *  http://www.library.cornell.edu/nr/bookcpdf/c3-6.pdf
 *  http://encyclopedia.thefreedictionary.com/Bilinear%20interpolation
 *  http://ct.radiology.uiowa.edu/~jiangm/courses/dip/html/node67.html
 *  http://bilinear-interpolation.wikiverse.org/
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
  "in=???\n       Input image file",
  "out=???\n      Output image file",
  "bad=0.0\n	  Value of a bad pixel to be patched",
  "goodx=\n       array of coluns in X that have good values (0=first)",
  "goody=\n       array of rows in Y that have good values (0=first)",
  "VERSION=0.2\n  29-oct-03 PJT",
  NULL,
};

string usage = "image bi-linear interpolation";

local void ini_fit(void), accum_fit(int, int, real);
local bool good_fit(int), lin_fit(int);
local real best_fit(void);

#define NLSQ  3
#define MAXN  256



void nemo_main(void)
{
  stream   instr, outstr;
  int      m, n, nx, ny, nz;        /* size of scratch map */
  int      ngood, ntry, nlin;
  int      i,j, di, dj;
  imageptr iptr=NULL, iptr1=NULL;      /* pointer to images */
  real     crit = getdparam("bad");
  int      irrx[MAXN],  irry[MAXN];         /* good rows and columns */
  int      nok, nbad, nirrx, nirry;
  int      ix,iy, i0,i1,j0,j1;
  
  instr = stropen(getparam("in"), "r");
  read_image( instr, &iptr);
  nx = Nx(iptr);	
  ny = Ny(iptr);
  nz = Nz(iptr);
  if (nz > 1) error("nz=%d :: Can't handle cubes yet",nz);
  create_image(&iptr1,nx,ny);
  copy_image(iptr,&iptr1);
  outstr = stropen(getparam("out"), "w");

  nirrx = nemoinpi(getparam("goodx"),irrx,MAXN);
  nirry = nemoinpi(getparam("goody"),irry,MAXN);

  if (nirrx == 0) {      /* if no goodx given, try to find them  */
    for (j=0; j<ny; j++) { 
      for (i=0; i<nx; i++) {
	if (MapValue(iptr,i,j) != crit) {
	  irrx[nirrx++] = i;
	}
      }
      if (nirrx) {
	dprintf(0,"found %d goodx values in row %d [%d .. %d]\n",
		nirrx,j,irrx[0],irrx[nirrx-1]);
	break;
      }
    }
  }

  if (nirry == 0) {      /* if no goody given, try to find them */
    for (i=0; i<nx; i++) {
      for (j=0; j<ny; j++) { 
	if (MapValue(iptr,i,j) != crit) {
	  irry[nirry++] = j;
	}
      }
      if (nirry) {
	dprintf(0,"found %d goody values in column %d [%d .. %d]\n",
		nirry,i,irry[0],irry[nirry-1]);
	break;
      }
    }
  }

  nbad=ngood=0;
  for (j=0; j<nirry; j++) { 
    for (i=0; i<nirrx; i++) {
      if (MapValue(iptr,irrx[i],irry[j]) != crit)
	nbad++;
      else
	nok++;
    }
  }
  if (nbad) warning("Not all good data is good: %d\n",nbad);
  if (nok)  dprintf(0,"ok, indeed found %d good values\n",ngood);
  

  for (j=0; j<ny; j++) {             /* copy array */
    for (i=0; i<nx; i++) {
      MapValue(iptr1,i,j) = MapValue(iptr,i,j);
    }
  }

  /* patch up rows - use linear interpolation */

  for (iy=0; iy<nirry; iy++) {
    for (ix=0; ix<nirrx-1; ix++) {
      i0=irrx[ix];
      i1=irrx[ix+1];
      j=irry[iy];
      for (i=i0+1; i<i1; i++) {
	MapValue(iptr1,i,j) = ( (i1-i)*MapValue(iptr,i0,j) + (i-i0)*MapValue(iptr,i1,j) )/(i1-i0);
      }
    }
  }

  /* patch up columns - use linear interpolation */

  for (ix=0; ix<nirrx; ix++) {
    for (iy=0; iy<nirry-1; iy++) {
      j0=irry[iy];
      j1=irry[iy+1];
      i=irrx[ix];
      for (j=j0+1; j<j1; j++) {
	MapValue(iptr1,i,j) = ( (j1-j)*MapValue(iptr,i,j0) + (j-j0)*MapValue(iptr,i,j1) )/(j1-j0);
      }
    }
  }

  /* patch up the middle pieces - use bi-linear interpolation */

#if 1
  for (iy=0; iy<nirry-1; iy++) {
    j0=irry[iy];
    j1=irry[iy+1];
    for (ix=0; ix<nirrx-1; ix++) {
      i0=irrx[ix];
      i1=irrx[ix+1];
      for (j=j0+1; j<j1; j++) {
	for (i=i0+1; i<i1; i++) {
	  if (ix==0 && iy==0) { 
	    /* debug first patch */
	    dprintf(1,"%d %d %d %d : %g %g %g %g  %g %g %g %g\n", i, j, i1-i0, j1-j0,
		   (i1-i)*(j1-j)*MapValue(iptr,i0,j0),
		   (i1-i)*(j-j0)*MapValue(iptr,i0,j1),
		   (i-i0)*(j1-j)*MapValue(iptr,i1,j0),
		   (i-i0)*(j-j0)*MapValue(iptr,i1,j1),
		   MapValue(iptr,i0,j0),
		   MapValue(iptr,i0,j1),
		   MapValue(iptr,i1,j0),
		   MapValue(iptr,i1,j1));

	  }
	  MapValue(iptr1,i,j) = (
				 (i1-i)*(j1-j)*MapValue(iptr,i0,j0) + 
				 (i1-i)*(j-j0)*MapValue(iptr,i0,j1) +
				 (i-i0)*(j1-j)*MapValue(iptr,i1,j0) +
				 (i-i0)*(j-j0)*MapValue(iptr,i1,j1) ) /
	    ((j1-j0)*(i1-i0));
	}
      }
    }
  }
#endif
  write_image(outstr, iptr1);
}


local int nsum;
local real lmat[NLSQ*NLSQ], lvec[NLSQ], lsol[NLSQ], la[NLSQ+1];

local void ini_fit(void)
{
    nsum = 0;
    lsq_zero(NLSQ,lmat,lvec);
}

local void accum_fit(int di, int dj, real val)
{
    nsum++;
    la[0] = 1.0;
    la[1] = (real) di;
    la[2] = (real) dj;
    la[3] = val;
    lsq_accum(NLSQ,lmat,lvec,la,1.0);
}

local bool good_fit(int m)
{
  return nsum>=m;
}

local bool lin_fit(int n)
{
  return nsum==2;
}


local real best_fit(void)
{
    lsq_solve(NLSQ,lmat,lvec,lsol);
    return lsol[0];
}
