/* 
 * CCDFILL: patch in the holes in a CCD frame
 *
 *	10-dec-93 	written	- very quick and dirty      pjt
 *	14-mar-95	protos's
 *	22-feb-97	1.1 added bad=				pjt
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
	"VERSION=1.1\n  22-feb-97 PJT",
	NULL,
};

string usage = "patch up holes in an image";


local void ini_fit(void), accum_fit(int, int, real);
local bool good_fit(void);
local real best_fit(void);

void nemo_main(void)
{
    stream   instr, outstr;
    int      n, nx, ny, nz;        /* size of scratch map */
    int      ngood, ntry;
    int      i,j, di, dj;
    imageptr iptr=NULL, iptr1=NULL;      /* pointer to images */
    real     crit = getdparam("bad");

    instr = stropen(getparam("in"), "r");
    n = getiparam("n");

    read_image( instr, &iptr);

    nx = Nx(iptr);	
    ny = Ny(iptr);
    create_image(&iptr1,nx,ny);

    outstr = stropen(getparam("out"), "w");


    for (j=0; j<ny; j++) {                  /* make copy of image */
    	for (i=0; i<nx; i++) {
    	    MapValue(iptr1,i,j) = MapValue(iptr,i,j);
	}
    }

    ntry = ngood = 0;
    for (j=0; j<ny; j++) {                  /* loop over image and patch */
    	for (i=0; i<nx; i++) {
            if (MapValue(iptr1,i,j) == crit) {      /* try and patch up */
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
                if (good_fit()) {
                    MapValue(iptr,i,j) = best_fit();
                    ngood++;
                } 
            }
	}
    }
    write_image(outstr, iptr);
    dprintf(0,"Found %d values to patch, successfull with %d\n",ntry,ngood);
}


#define NLSQ  3

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

local bool good_fit(void)
{
    return nsum>=NLSQ;
}

local real best_fit(void)
{
    lsq_solve(NLSQ,lmat,lvec,lsol);
    return lsol[0];
}
