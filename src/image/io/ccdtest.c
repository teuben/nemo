/* 
 * CCDPRINT: print values at gridpoints, with optional numeric labels
 *
 *	 23-jan-89  V1.0 created        Peter Teuben
 *	 15-jul-89  V1.0a   - nemoinpX instead of nemoinp	PJT
 *	  8-jul-93  V1.1 blank default for x=,y=,z=             pjt
 *		    and  fixep various old style coding
 *	 14-sep-93  V1.2 Added index mode keyword offset=	pjt
 *	 14-oct-99  V1.2b fixed bug printing WCS labels		pjt
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
	"in=???\n        Input filename",
	"VERSION=1.2b\n     14-oct-99 PJT",
	NULL,
};

string usage = "print values at gridpoints of an image";


nemo_main()
{
    int     *ix, *iy, *iz, i, j, k, nx, ny, nz, nxpos, nypos, nzpos, offset;
    bool    newline, xlabel, ylabel, zlabel;
    real    f;
    string  infile;			        /* file name */
    stream  instr;				/* file stream */
    imageptr iptr=NULL;			      /* allocated dynamically */
    string   fmt, label;
    real     scale_factor, x, y, z;
    real     **a;
    
    instr = stropen (getparam("in"), "r");
    if (read_image (instr,&iptr) == 0)
    	error("Problem reading image");
    
    nx = Nx(iptr);	                        /* cube dimensions */
    ny = Ny(iptr);
    a = map2_image(iptr);
    for (j=0; j<ny; j++)
    for (i=0; i<nx; i++)
        printf(" %d %d -> %g\n",i,j,a[i][j]);

    strclose(instr);
}

