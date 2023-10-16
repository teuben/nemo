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
  "z=\n              Pixel range in Z select (2 values, or none for all pixels)",
  "scale=1,1\n       Scaling factors for the two columns",
  "VERSION=0.6\n     21-dec-2022 PJT",
  NULL,
};

string usage = "print spectrum at a grid point of a cube";


void nemo_main(void)
{
    int      ix, iy, iz, nx, ny, nz, ns, nzr, zr[2];
    string   infile;			        /* file name */
    stream   instr;				/* file stream */
    imageptr iptr=NULL;			      /* allocated dynamically */
    real     sf[2], x, y, z, f, f1, *data;
    Moment   m1, m2;

    instr = stropen (getparam("in"), "r");
    if (read_image (instr,&iptr) == 0)
      error("Problem reading image from in=",getparam("in"));
    strclose(instr);

    ns = nemoinpr(getparam("scale"),sf,2);
    if (ns != 2) error("Need two values for scale=%s",getparam("scale"));

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
    
    nzr = nemoinpi(getparam("z"),zr,2);
    if (nzr == 0) {
      zr[0] = 0;
      zr[1] = nz-1;
    } else if (nzr != 2)
      error("Two integer values required for z=%s", getparam("z"));
    if (zr[0] > zr[1]) {
      iz = zr[0];
      zr[0] = zr[1];
      zr[1] = iz;
    }
    if (zr[0] < 0)    zr[0] = 0;
    if (zr[1] > nz-1) zr[1] = nz-1;

    /* simple header */
    printf("# ccdspec %s  %d/%d %d/%d\n", getparam("in"), ix, nx, iy, ny);
    printf("# restfreq %f\n",Restfreq(iptr));
    printf("# %s %s\n", Namez(iptr), Unit(iptr));
    
    /* write spectrum */
    for (iz=zr[0], f1=0; iz<=zr[1]; iz++) {
        z = (Zmin(iptr) + (iz-Zref(iptr)) * Dz(iptr)) * sf[0];
	f = CubeValue(iptr,ix,iy,iz) * sf[1];
	printf("%g %g\n", z, f);
	// moment analysis on the value and the difference from previous
	accum_moment(&m1, f, 1.0);
	if (iz>0) accum_moment(&m2, f-f1, 1.0);
	f1 = f;
    }
    printf("# %d %d %d min/mean/sig/max  %g %g %g %g    %g %g   %g\n", 
	   ix, iy, zr[1]-zr[0]+1,
	   min_moment(&m1), mean_moment(&m1), sigma_moment(&m1), max_moment(&m1),
 	   mean_moment(&m2), sigma_moment(&m2),
	   sigma_moment(&m2)/sigma_moment(&m1)/sqrt(2));
}

