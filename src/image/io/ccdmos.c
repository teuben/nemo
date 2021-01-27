/* 
 *	CCDMOS: mosaic a regular pattern of NX by NY "identical" 
 *      ccd images into a big one
 *
 *      1-jul-2005   Created for CARMA
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <extstring.h>
#include <image.h>

string defv[] = {
  "in=???\n       Input filename(s)",
  "out=???\n      Output filename",
  "nx=\n          Number of rows",
  "ny=\n          Number of columns",
  "onefile=f\n    All images in one file?",
  "VERSION=1.0\n  1-jul-05 PJT",
  NULL,
};

string usage="mosaic a regular pattern of images";

string cvsid="$Id$";

extern string *burststring(string,string);

void nemo_main(void)
{
    imageptr iptr=NULL, optr=NULL;
    string *filename;
    stream instr, outstr;
    int i, j, i0, j0, nfiles;
    int nx, ny, nx1, ny1, ix, iy, n;

    nx = hasvalue("nx") ? getiparam("nx") : 0;
    ny = hasvalue("ny") ? getiparam("ny") : 0;
    
    filename = burststring(getparam("in")," ,\n");
    nfiles = xstrlen(filename,sizeof(string))-1;
    if (nx==0 && ny==0)
      error("Need at least one of nx= or ny=");
    else if (nx==0) {
      nx = nfiles/ny;
      dprintf(0,"nx=%d ny=%d nfiles=%d\n",nx,ny,nfiles);
      if (nfiles % ny) error("nfiles=%d/ny=%d does not divide evenly",nfiles,ny);
    } else if (ny==0) {
      ny = nfiles/nx;
      dprintf(0,"nx=%d ny=%d nfiles=%d\n",nx,ny,nfiles);
      if (nfiles % nx) error("nfiles=%d/nx=%d does not divide evenly",nfiles,nx);
    }

    for (iy=0, n=0; iy<ny; iy++) {
      for (ix=0; ix<nx; ix++, n++) {
	instr = stropen (filename[n],"r");
	read_image (instr,&iptr);               /* read image */
	strclose(instr);                        /* close image file */
	if (n==0) {
	  nx1 = Nx(iptr);
	  ny1 = Ny(iptr);
	  create_image(&optr,nx*nx1,ny*ny1);
	}
	i0 = ix*nx1;
	j0 = iy*ny1;
	for (j=0; j<ny1; j++)
	  for (i=0; i<nx1; i++)
	    MapValue(optr,i+i0,j+j0) = MapValue(iptr,i,j);
	i++;
      }
    }
    outstr = stropen(getparam("out"),"w");
    write_image(outstr,optr);
    strclose(outstr);
}

