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
  "VERSION=0.1\n  19-apr-2017 PJT",
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


    Xmin(optr) = Xmin(iptr) + xb[0] * Dx(iptr);
    Dx(optr)   = Dx(iptr);

    Ymin(optr) = Ymin(iptr) + yb[0] * Dy(iptr);
    Dy(optr)   = Dy(iptr);

    Zmin(optr) = Zmin(iptr) + zb[0] * Dz(iptr);
    Dz(optr)   = Dz(iptr);
    dprintf(0,"WCS Corner: %g %g %g\n",Xmin(optr),Ymin(optr),Zmin(optr));

    Namex(optr) = Namex(iptr);
    Namey(optr) = Namey(iptr);
    Namez(optr) = Namez(iptr);
    
    write_image(outstr, optr);

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


void ax_copy(imageptr i0, imageptr i1)
{
  Dx(i1) = Dx(i0);
  Dy(i1) = Dy(i0);
  Dz(i1) = Dz(i0);
  Xmin(i1) = Xmin(i0);
  Ymin(i1) = Ymin(i0);
  Zmin(i1) = Zmin(i0);
  Xref(i1) = Xref(i0);
  Yref(i1) = Yref(i0);
  Zref(i1) = Zref(i0);
  Beamx(i1) = Beamx(i0);
  Beamy(i1) = Beamy(i0);
  Beamz(i1) = Beamz(i0);
  if (Namex(i0))  Namex(i1) = strdup(Namex(i0));
  if (Namey(i0))  Namey(i1) = strdup(Namey(i0));
  if (Namez(i0))  Namez(i1) = strdup(Namez(i0));
}
