/* 
 * CCDPRINT: print values at gridpoints, with optional numeric labels
 *
 *	 23-jan-89  V1.0 created        Peter Teuben
 *	 15-jul-89  V1.0a   - nemoinpX instead of nemoinp	PJT
 *	  8-jul-93  V1.1 blank default for x=,y=,z=             pjt
 *		    and  fixep various old style coding
 *	 14-sep-93  V1.2 Added index mode keyword offset=	pjt
 *	 14-oct-99  V1.2b fixed bug printing WCS labels		pjt
 *       26-sep-2013  V1.3  print out 2D and 3D  correct for CDEF   PJT
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
	"in=???\n        Input filename",
	"VERSION=1.3\n   26-sep-2013 PJT",
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
    real     **a, ***b;
    
    instr = stropen (getparam("in"), "r");
    if (read_image (instr,&iptr) == 0)
    	error("Problem reading image");
    
    nx = Nx(iptr);	                        /* cube dimensions */
    ny = Ny(iptr);
    nz = Nz(iptr);
    if (nz>1) {
      printf("# nx,ny,nz: %d %d %d\n",nx,ny,nz);
      b = map3_image(iptr);                       /* get the first slice */
      for (k=0; k<nz; k++)
	for (j=0; j<ny; j++)
	  for (i=0; i<nx; i++)
	    printf(" B[%d][%d][%d] -> %g\n",k,j,i,b[k][j][i]); // CDEF
    } else {
      printf("# nx,ny %d %d\n",nx,ny);
      a = map2_image(iptr);                       /* get the first slice */
      for (j=0; j<ny; j++)
	for (i=0; i<nx; i++)
	  //printf(" A[%d][%d] -> %g\n",j,i,a[i][j]); // CDEF
	  printf(" A[%d][%d] -> %g\n",j,i,a[j][i]); // FORDEF
    }
    //strclose(instr);
}

