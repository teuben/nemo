/* 
 * CCDSUB: get a subimage from an image
 *
 *	quick and dirty: 10-aug-93
 *	something...	  6-sep-93
 *	n[xy]aver=	  5-nov-93 	puzzling.....
 *		something is wrong here                     
 * 1.3  added some z stuff for Martin Bureau       1-may-2002      PJT
 * 1.4  added reorder=

    TODO:  wcs is wrong on output
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
        "in=???\n       Input image file",
	"out=???\n      Output image file",
	"x=\n		X coordinates to select (1..nx)",
	"y=\n		Y coordinates to select (1..ny)",
	"z=\n		Z coordinates to select (1..nz)",
	"aver=\n	Coordinates to average",
	"nxaver=1\n	Number X to aver (size remains same)",
	"nyaver=1\n	Number Y to aver (size remains same)",
	"nzaver=1\n	Number Z to aver (size remains same)",
	"skip=xyz\n	Coordinates to skip (** ignored **)",
	"reorder=\n     New coordinate ordering",
	"VERSION=1.4a\n  12-aug-2002 PJT",
	NULL,
};

string usage = "sub/average of an image";


#define MAXDIM 8192

int  ix[MAXDIM], iy[MAXDIM], iz[MAXDIM];

local int ax_index(string , int , int , int *);

void nemo_main()
{
    stream  instr, outstr;
    int     nx, ny, nz;        /* size of scratch map */
    int     nx1,ny1,nz1;
    int     nxaver, nyaver,nzaver;
    int     i,j,k, i0,j0,k0, i1,j1,k1;
    imageptr iptr=NULL, iptr1=NULL;      /* pointer to images */
    real    sum, tmp, zzz;
    bool    Qreorder;
    string  reorder;

    instr = stropen(getparam("in"), "r");
    nxaver=getiparam("nxaver");
    nyaver=getiparam("nyaver");
    nzaver=getiparam("nzaver");

    nx1 = nemoinpi(getparam("x"),ix,MAXDIM);
    ny1 = nemoinpi(getparam("y"),iy,MAXDIM);
    nz1 = nemoinpi(getparam("z"),iz,MAXDIM);
    if (nx1<0 || ny1<0 || nz1<0) error("Error parsing x,y,z=");

    read_image( instr, &iptr);

    nx = Nx(iptr);	
    ny = Ny(iptr);
    nz = Nz(iptr);
    if (hasvalue("reorder")) {
#if 0    
      reorder = getparam("reorder");
      if (strlen(reorder) != 3) error("Reorder must have 3 letters");
      for (i=0; i<3; i++)
	switch (reorder[i]) {
	case 'x','X':
	case 'y','Y':
	case 'z','Z':
	default:
	  error("...");
	}
#else
	  error("...this option not implemented yet...");
#endif
    } else {
      nx1 = ax_index("x",nx,nx1,ix);
      ny1 = ax_index("y",ny,ny1,iy);
      nz1 = ax_index("z",nz,nz1,iz);
    }

    outstr = stropen(getparam("out"), "w");

    if (nxaver>1 || nyaver>1 || nzaver>1) {  /*  averaging, but retaining size */
        dprintf(0,"Averaging map %d * %d * %d pixels; mapsize %d * %d * %d\n",
                   nxaver,nyaver,nzaver,nx,ny,nz);
        nx1 = nx/nxaver;
        ny1 = ny/nyaver;
        nz1 = nz/nzaver;
        for (k1=0; k1<nz1; k1++) {
	  k = k1*nzaver;
	  for (j1=0; j1<ny1; j1++) {
            j = j1*nyaver;
            for (i1=0; i1<nx1; i1++) {
                i = i1*nxaver;
                sum = 0.0;
                for (j0=0; j0<nyaver; j0++)
                    for (i0=0; i0<nxaver; i0++)
                        sum += CubeValue(iptr,i+i0, j+j0, k+k0);
                sum /= (real) (nxaver*nyaver*nzaver);
                for (k0=0; k0<nzaver; k0++)
		  for (j0=0; j0<nyaver; j0++)
                    for (i0=0; i0<nxaver; i0++)
                        CubeValue(iptr,i+i0, j+j0, k+k0) = sum;
            }
	  }
	}
        write_image(outstr, iptr);
    } else if (Qreorder) {            	/* reordering */
        create_cube(&iptr1,nx1,ny1,nz1);
        for (k=0; k<nz1; k++)
	  for (j=0; j<ny1; j++)
            for (i=0; i<nx1; i++)
                CubeValue(iptr1,i,j,k) = CubeValue(iptr,ix[i],iy[j],iz[k]);
        write_image(outstr, iptr1);
    } else {            	/* straight sub-sampling */
        create_cube(&iptr1,nx1,ny1,nz1);
        for (k=0; k<nz1; k++)
	  for (j=0; j<ny1; j++)
            for (i=0; i<nx1; i++)
                CubeValue(iptr1,i,j,k) = CubeValue(iptr,ix[i],iy[j],iz[k]);
	/* need to do WCS here */
        write_image(outstr, iptr1);
    }                    
}

/*
 * either initialize idx array, if not done, or normalize to 0..n-1
 */
 
int ax_index(string name, int n, int n1, int *idx)
{
    int i;
    
    if (n1==0) {		/* copy array */
        n1=n;
        for (i=0; i<n; i++) idx[i] = i;
    } else {
        for (i=0; i<n1; i++) {
            if (idx[i] < 1 || idx[i] > n)
                error("Index %d illegal in %s axis; max %d",idx[i],name,n);
            idx[i] -= 1;
        }
    }
    return n1;
}
