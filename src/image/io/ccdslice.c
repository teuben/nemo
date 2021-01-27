/* 
 *	CCDSLICE: take slices from a cube
 *
 *	 6-may-95  V1.0 created 	Peter Teuben
 *      29-dec-01  V1.0a   also compute MapMin/Max
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>

#include <image.h>

string defv[] = {	/* keywords + help string for user interface */
    "in=???\n       Input filename (image)",
    "out=???\n      Output filename (image)",
    "zvar=z\n       Slice variable (x,y,z)",
    "zrange=\n      Slices to select (1..n<zvar>)",
    "VERSION=1.0a\n 29-dec-01 PJT",
    NULL,
};

string usage="takes slices from a cube";

#ifndef MAXPLANE
#define MAXPLANE 2048
#endif

#define X_SLICE 0
#define Y_SLICE 1
#define Z_SLICE 2

static int need_init_minmax = 1;
static real new_min, new_max;

void new_minmax(real ival);
void slice(imageptr i, imageptr o, int mode, int *planes);


void nemo_main(void)
{
    imageptr iptr=NULL, optr=NULL;
    stream instr, outstr;
    int planes[MAXPLANE], mode, i;
    string zvar;
    int nx, ny, nz;

    instr = stropen (getparam("in"),"r");	/* get stream */
    read_image (instr,&iptr);               /* read image */
    strclose(instr);                        /* close image file */

    zvar = getparam("zvar");
    if (streq(zvar,"x")) {
        mode = X_SLICE;
        nx = Ny(iptr);
        ny = Nz(iptr);
        nz = Nx(iptr);
    } else if (streq(zvar,"y")) {
        mode = Y_SLICE;
        nx = Nx(iptr);
        ny = Nz(iptr);
        nz = Ny(iptr);
    } else if (streq(zvar,"z")) {
        mode = Z_SLICE;
        nx = Nx(iptr);
        ny = Ny(iptr);
        nz = Nz(iptr);
    } else 
        error("Illegal slice axis %s, must be x,y,z",zvar);


    if (hasvalue("zrange")) {
        nz = nemoinpi(getparam("zrange"),planes,MAXPLANE);
        if (nz<1) 
            error("Bad syntax %d for zrange=%s",nz,getparam("zrange"));
    } else {
        for (i=0; i<nz; i++)
            planes[i] = i+1;
    }
    create_cube(&optr,nx,ny,nz);
    slice(iptr,optr,mode,planes);

    outstr = stropen(getparam("out"),"w");
    write_image(outstr,optr);
    strclose(outstr);
}

void new_minmax(real ival) 
{
  if (need_init_minmax) {
    need_init_minmax = 0;
    new_min = new_max = ival;
    return;
  } else {
    new_min = MIN(new_min, ival);
    new_max = MAX(new_max, ival);
  }
}

void slice(imageptr i, imageptr o, int mode, int *planes)
{
    int x, y, z, iz;
    real ival;

    warning("Code not converted to fix reference pixel value");

    if (mode==X_SLICE) {
        for(iz=0; iz<Nz(o); iz++) {
            z = planes[iz]-1;
            if (z<0 || z>Nx(i)) 
                error("%d: illegal plane in x, max is %d",z,Nx(i));
            for (y=0; y<Ny(o); y++)
	      for (x=0; x<Nx(o); x++) {
		ival = CubeValue(i,z,x,y);
                CubeValue(o,x,y,iz) = ival;
		new_minmax(ival);
	      }
        }
        Namex(o) = Namey(i);
        Namey(o) = Namez(i);
        Namez(o) = Namex(i);
        Xmin(o) = Ymin(i);
        Ymin(o) = Zmin(i);
        Zmin(o) = Xmin(i);
        Dx(o) = Dy(i);
        Dy(o) = Dz(i);
        Dz(o) = Dx(i);
    } else if (mode==Y_SLICE) {
        for(iz=0; iz<Nz(o); iz++) {
            z = planes[iz]-1;
            if (z<0 || z>Ny(i)) 
                error("%d: illegal plane in y, max is %d",z,Ny(i));
            for (y=0; y<Ny(o); y++)
	      for (x=0; x<Nx(o); x++) {
		ival = CubeValue(i,x,z,y);
                CubeValue(o,x,y,iz) = ival;
		new_minmax(ival);
	      }
        }
        Namex(o) = Namex(i);
        Namey(o) = Namez(i);
        Namez(o) = Namey(i);
        Xmin(o) = Xmin(i);
        Ymin(o) = Zmin(i);
        Zmin(o) = Ymin(i);
        Dx(o) = Dx(i);
        Dy(o) = Dz(i);
        Dz(o) = Dy(i);
    } else if (mode==Z_SLICE) {
        for(iz=0; iz<Nz(o); iz++) {
            z = planes[iz]-1;
            if (z<0 || z>Nz(i)) 
                error("%d: illegal plane in z, max is %d",z,Nz(i));
            for (y=0; y<Ny(o); y++)
	      for (x=0; x<Nx(o); x++) {
		ival = CubeValue(i,x,y,z);
                CubeValue(o,x,y,iz) = ival;
		new_minmax(ival);
	      }
        }
        Namex(o) = Namex(i);
        Namey(o) = Namey(i);
        Namez(o) = Namez(i);
        Xmin(o) = Xmin(i);
        Ymin(o) = Ymin(i);
        Zmin(o) = Zmin(i);
        Dx(o) = Dx(i);
        Dy(o) = Dy(i);
        Dz(o) = Dz(i);
    }
    MapMin(o) = new_min;
    MapMax(o) = new_max;
}
