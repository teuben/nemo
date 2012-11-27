/* 
 * CCDPRINT: print info from an image
 *
 *	 27-nov-2012  V1.0 created        Peter Teuben
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
	"in=???\n         Input filename",
	"VERSION=1.0\n    27-nov-2012 PJT",
	NULL,
};

string usage = "print info from an image";


nemo_main()
{
  int     i, j, k, nx, ny, nz;
  string  infile;			        /* file name */
  stream  instr;				/* file stream */
  imageptr iptr=NULL;			      /* allocated dynamically */
    
  instr = stropen (getparam("in"), "r");
  if (read_image (instr,&iptr) == 0)
    error("Problem reading image");
    
  nx = Nx(iptr);	                        /* cube dimensions */
  ny = Ny(iptr);
  nz = Nz(iptr);
  
  printf("%d %d %d\n",nx,ny,nz);
  
  strclose(instr);
}

