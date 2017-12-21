/* 
 * CCDFILL: patch in the holes (or spikes) in a CCD frame
 *
 *	10-dec-93 	written	- very quick and dirty      pjt
 *	14-mar-95	protos's
 *	22-feb-97	1.1 added bad=				pjt
 *      19-feb-02       1.2 added all=, m=			pjt
 *       2-may-03       1.3 added iter=                         pjt
 *       2-jun-03       1.4 trying out xmirror/ymirror          pjt
 *       1-aug-05       1.5 make it work on cubes too           pjt
 *       2-may-2017     2.0 also allow deviate pixels values    pjt
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
	"n=1\n		Number of neighbor cells on all sides to use",
	"bad=0.0\n	Value of a bad pixel to be patched",
	"all=f\n        Force all points to be refitted?",
	"m=3\n          Minimum number of neighbor pixels needed",
	"iter=1\n       Number of iterations",
	"xmirror=f\n    Use points mirrored in X to get near a border",
	"ymirror=f\n    Use points mirrored in Y to get near a border",
	"spike=\n       Spike value above the neighbors to be filled",
	"VERSION=2.0\n  2-may-2017 PJT",
	NULL,
};

string usage = "patch up holes (or spikes) in an image by linear interpolation";

string cvsid = "$Id$";


local void ini_fit(void), accum_fit(int, int, real);
local bool good_fit(int), lin_fit(int);
local real best_fit(void);

#define NLSQ  3

bool bad_mv(real mv1,real mv2,real mv3,real mv4,real mv,real crit,real spike) {
  int i,ngood = 0;    /* have to loop manually over mv1...mv4 */
  real check[4];
  if (mv1 != crit) check[ngood++] = mv1;
  if (mv2 != crit) check[ngood++] = mv2;
  if (mv3 != crit) check[ngood++] = mv3;
  if (mv4 != crit) check[ngood++] = mv4;
  if (ngood==0) return FALSE;
  for (i=0; i<ngood; i++)
    if (ABS(check[i]-mv) < spike) return FALSE;
  return TRUE;
}


void nemo_main(void)
{
  stream   instr, outstr;
  int      m, n, nx, ny, nz;        /* size of scratch map */
  int      ngood, ntry, nlin, nspike;
  int      i,j,k, di, dj;
  imageptr iptr=NULL, iptr1=NULL;      /* pointer to images */
  real     mv, crit = getdparam("bad");
  bool     Qall = getbparam("all");
  bool     Qxm  = getbparam("xmirror");
  bool     Qym  = getbparam("ymirror");
  bool     Qspike = hasvalue("spike");
  int      iter, niter = getiparam("iter");
  real     spike;
  
  instr = stropen(getparam("in"), "r");
  n = getiparam("n");
  m = getiparam("m");
  if (m < NLSQ) error("Cannot choose m=%d < %d",m,NLSQ);
  if (Qspike) spike = getrparam("spike");
  
  read_image( instr, &iptr);
  
  nx = Nx(iptr);	
  ny = Ny(iptr);
  nz = Nz(iptr);
  create_image(&iptr1,nx,ny);     /* single temp plane to work in */
  
  outstr = stropen(getparam("out"), "w");

  /* first we find spikes away from the mean, and replace them with 0.0 to force a patch... */

  if (Qspike) {
    nspike = 0;
    for (k=0; k<nz; k++) {
      /* first make a copy of this plane to work in - patch in the cube */
      for (j=0; j<ny; j++) {                  
    	for (i=0; i<nx; i++) {
	  MapValue(iptr1,i,j) = CubeValue(iptr,i,j,k);
	}
      }
      
      for (j=1; j<ny-1; j++) {                  
    	for (i=1; i<nx-1; i++) {
	  mv = MapValue(iptr1,i,j);
	  if (mv == crit) continue;
	  if (bad_mv(MapValue(iptr1,i,j-1),
		     MapValue(iptr1,i,j+1),
		     MapValue(iptr1,i-1,j),
		     MapValue(iptr1,i+1,j),
		     mv, crit, spike)) {
	    CubeValue(iptr,i,j,k) = crit;
	    nspike++;
	    continue;
	  }
	}
      }
    }
    dprintf(0,"Found %d spike=%g values to reset\n",nspike, spike);
  }


  /* now fill the 0's */
  
  for (k=0; k<nz; k++) {
    for (iter=0; iter < niter; iter++) {
      
      for (j=0; j<ny; j++) {                  /* make copy of image */
    	for (i=0; i<nx; i++) {
	  MapValue(iptr1,i,j) = CubeValue(iptr,i,j,k);
	}
      }
      
      ntry = ngood = nlin = 0;
      for (j=1; j<ny-1; j++) {                  /* loop over image and patch */
    	for (i=1; i<nx-1; i++) {
	  if (MapValue(iptr1,i,j) == crit || Qall) {   /* try and patch up */
	    ntry++;
	    ini_fit();
	    for (dj=-n; dj<=n; dj++) {
	      if (j+dj<0 || j+dj>=ny) continue;
	      for (di=-n; di<=n; di++) {
		if (i+di<0 || i+di>=nx) continue;
		if (MapValue(iptr1,i+di,j+dj) != crit)
		  accum_fit(di,dj,MapValue(iptr1,i+di,j+dj));
	      }
	    }
	    if (good_fit(m)) {
	      CubeValue(iptr,i,j,k) = best_fit();
	      ngood++;
	    } else if (lin_fit(n)) {
	      if (MapValue(iptr1,i+1,j) != crit && MapValue(iptr1,i-1,j) != crit) {
		nlin++;
		CubeValue(iptr,i,j,k) = 0.5 * (MapValue(iptr1,i+1,j) + MapValue(iptr1,i-1,j));
	      }
	      if (MapValue(iptr1,i,j+1) != crit && MapValue(iptr1,i,j-1) != crit) {
		nlin++;
		CubeValue(iptr,i,j,k) = 0.5 * (MapValue(iptr1,i,j+1) + MapValue(iptr1,i,j-1));
	      }
	    }
	  }
	} /* i */
      } /* j */
      dprintf(0,"Found %d values to patch, successfull with %d, %d linear?\n",ntry,ngood,nlin);
      if (ngood == 0 || ntry==ngood) break;
    } /* iter */
  } /* k */
  write_image(outstr, iptr);
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
