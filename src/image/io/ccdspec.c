/* 
 * CCDSPEC: special version of ccdprint to take a spectrum and print ascii
 *          with optional stats and difference stats
 *
 *       11-feb-2021    Q&D
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>
#include <moment.h>

string defv[] = {
  "in=???\n          Input filename",
  "x=\n              Pixel in X to print (0=1st pixel)",
  "y=\n              Pixel in Y to print",
  "VERSION=0.2\n     11-feb-2021 PJT",
  NULL,
};

string usage = "print spectrum at a grid point";

string cvsid="$Id$";


void nemo_main(void)
{
    int      ix, iy, iz, nx, ny, nz;
    string   infile;			        /* file name */
    stream   instr;				/* file stream */
    imageptr iptr=NULL;			      /* allocated dynamically */
    real     scale_factor, x, y, z, f, f1, *data;
    Moment   m1, m2;

    instr = stropen (getparam("in"), "r");
    if (read_image (instr,&iptr) == 0)
      error("Problem reading image from in=",getparam("in"));
    strclose(instr);

    ini_moment(&m1, 2, 0);
    ini_moment(&m2, 2, 0);

    nx = Nx(iptr);	
    ny = Ny(iptr);
    nz = Nz(iptr);
    if (!hasvalue("x"))
      ix = nx/2;
    else
      ix = getiparam("x");
    if (!hasvalue("y"))
      iy = ny/2;
    else
      iy = getiparam("y");

    for (iz=0, f1=0; iz<nz; iz++) {
        z = Zmin(iptr) + (iz-Zref(iptr)) * Dz(iptr);
	f = CubeValue(iptr,ix,iy,iz);
	printf("%g %g   %g\n", z, f, f-f1);
	// moment analysis on the value and the difference from previous
	accum_moment(&m1, f, 1.0);
	if (iz>0) accum_moment(&m2, f-f1, 1.0);
	f1 = f;
    }
    printf("# %d %d mean/sig  %g %g    %g %g   %g\n", 
	   ix, iy,
	   mean_moment(&m1), sigma_moment(&m1),
 	   mean_moment(&m2), sigma_moment(&m2),
	   sigma_moment(&m2)/sigma_moment(&m1)/sqrt(2));
}

