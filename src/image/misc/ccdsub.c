/* 
 * CCDSUB: get a subimage from an image
 *
 *	quick and dirty: 10-aug-93
 *	something...	  6-sep-93
 *	n[xy]aver=	  5-nov-93 	puzzling.....
 * 
		something is wrong here                     
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
        "in=???\n       Input image file",
	"out=???\n      Output image file",
	"x=\n		X coordinates to select",
	"y=\n		Y coordinates to select",
	"z=\n		Z coordinates to select",
	"aver=\n	Coordinates to average",
	"nxaver=1\n	Number X to aver (size remains same)",
	"nyaver=1\n	Number Y to aver (size remains same)",
	"skip=xyz\n	Coordinates to skip",
	"VERSION=1.2\n  5-nov-93 PJT",
	NULL,
};

string usage = "sub/average of an image";


#define MAXDIM 4096

int  ix[MAXDIM], iy[MAXDIM], iz[MAXDIM];

local int ax_index(string , int , int , int *);

void nemo_main()
{
    stream  instr, outstr;
    int     nx, ny, nz;        /* size of scratch map */
    int     nx1,ny1,nz1;
    int     nxaver, nyaver;
    int     i,j, i0,j0, i1,j1;
    imageptr iptr=NULL, iptr1=NULL;      /* pointer to images */
    real    sum, tmp, zzz;

    instr = stropen(getparam("in"), "r");
    nxaver=getiparam("nxaver");
    nyaver=getiparam("nyaver");

    nx1 = nemoinpi(getparam("x"),ix,MAXDIM);
    ny1 = nemoinpi(getparam("y"),iy,MAXDIM);
    nz1 = nemoinpi(getparam("z"),iz,MAXDIM);
    if (nx1<0 || ny1<0 || nz1<0) error("Error parsing x,y,z=");

    read_image( instr, &iptr);

    nx = Nx(iptr);	
    ny = Ny(iptr);

    nx1 = ax_index("x",nx,nx1,ix);
    ny1 = ax_index("y",ny,ny1,iy);

    outstr = stropen(getparam("out"), "w");

    if (nxaver>1 || nyaver>1) {      /*  averaging, but retaining size */
        dprintf(0,"Averaging map %d *%d pixels; mapsize %d * %d\n",
                   nxaver,nyaver,nx,ny);
        nx1 = nx/nxaver;
        ny1 = ny/nyaver;
        for (j1=0; j1<ny1; j1++) {
            j = j1*nyaver;
            for (i1=0; i1<nx1; i1++) {
                i = i1*nxaver;
                sum = 0.0;
                for (j0=0; j0<nyaver; j0++)
                    for (i0=0; i0<nxaver; i0++)
                        sum += MapValue(iptr,i+i0, j+j0);
                sum /= (real) (nxaver*nyaver);
                for (j0=0; j0<nyaver; j0++)
                    for (i0=0; i0<nxaver; i0++)
                        MapValue(iptr,i+i0, j+j0) = sum;
            }
        }
        write_image(outstr, iptr);
    } else {            	/* straight sub-sampling */
        create_image(&iptr1,nx1,ny1);
        for (j=0; j<ny1; j++)
            for (i=0; i<nx1; i++)
                MapValue(iptr1,i,j) = MapValue(iptr,ix[i],iy[j]);
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
