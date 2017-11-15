/* 
 * CCDBORDER: put a border around an image 
 *
 *  0.1    PJT    quick and dirty again
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
  "in=???\n       Input image file",
  "out=???\n      Output image file",
  "x=0,0\n	  X pixels to add before and after",
  "y=0,0\n	  Y pixels to add before and after",
  "z=0,0\n	  Z pixels to add before and after",
  "value=0\n      value in the border region",
  "VERSION=0.2\n  19-apr-2017 PJT",
  NULL,
};

string usage = "add border around an image/cube";

string cvsid="$Id$";




#define LOOP(i,n)     for(i=0;i<n;i++)
#define CV(p,i,j,k)   CubeValue(p,i,j,k)

void nemo_main()
{
    stream  instr, outstr;
    int     nx, ny, nz;        /* size of scratch map */
    int     nx1,ny1,nz1;
    int     xb[2], yb[2], zb[2];
    int     i,j,k, i0,j0,k0, i1,j1,k1, l;
    imageptr iptr=NULL, optr=NULL;      /* pointer to images */
    real    value;
    real    *row;


    instr = stropen(getparam("in"), "r");
    value = getrparam("value");
    
    if (nemoinpi(getparam("x"),xb,2) != 2) error("Need 2 values for x");
    if (nemoinpi(getparam("y"),yb,2) != 2) error("Need 2 values for y");
    if (nemoinpi(getparam("z"),zb,2) != 2) error("Need 2 values for z");

    /* TODO: make sure no negative values for the borders */

    read_image( instr, &iptr);
    outstr = stropen(getparam("out"), "w");

    nx = Nx(iptr);	                   /* old cube size */
    ny = Ny(iptr);      
    nz = Nz(iptr);      
    nx1 = nx + xb[0] + xb[1];
    ny1 = ny + yb[0] + yb[1];
    nz1 = nz + zb[0] + zb[1];

    create_cube(&optr, nx1, ny1, nz1);
    LOOP(k,nz1) LOOP(j,ny1) LOOP(i,nx1) CV(optr,i,j,k) = 0.0;
    LOOP(k,nz)  LOOP(j,ny)  LOOP(i,nx)  CV(optr,i+xb[0],j+yb[0],k+zb[0]) = CV(iptr,i,j,k);


    Xmin(optr) = Xmin(iptr) - xb[0] * Dx(iptr);
    Dx(optr)   = Dx(iptr);

    Ymin(optr) = Ymin(iptr) - yb[0] * Dy(iptr);
    Dy(optr)   = Dy(iptr);

    Zmin(optr) = Zmin(iptr) - zb[0] * Dz(iptr);
    Dz(optr)   = Dz(iptr);
    dprintf(0,"WCS Corner: %g %g %g\n",Xmin(optr),Ymin(optr),Zmin(optr));

    Namex(optr) = Namex(iptr);
    Namey(optr) = Namey(iptr);
    Namez(optr) = Namez(iptr);
    
    write_image(outstr, optr);

}
