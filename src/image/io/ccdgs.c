/* 
 *	CCDGS:   play with gather-scatter operations on large blocks of data
 *
 *      gather:   loop(i,N) x[i] = y[idx[i]];
 *      scatter:  loop(i,N) y[idx[i]] = x[i];
 *
 *      24-jul-2013   Tinkertoy'd
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>

/* for readv() and writev() */
#include <sys/uio.h>

string defv[] = {
  "in=???\n       Input filename (2d/3d block of data)",
  "out=???\n      Output filename (2d/3d block of data)",
  "size=10,10,1\n Size in nx,ny[,nz]",
  "offset=0,0\n   Offset in input and output file (**not used**)",
  "VERSION=0.1\n  31-jul-2013 PJT",
  NULL,
};

string usage="gather-scatter operations on large blocks of data";

void gs_xy(stream instr, stream outstr,int offset,int nx,int ny);
void gs_yx(stream instr, stream outstr,int offset,int nx,int ny);

void gs_xyz(stream instr, stream outstr,int offset,int nx,int ny,int nz);
void gs_yxz(stream instr, stream outstr,int offset,int nx,int ny,int nz);
void gs_zxy(stream instr, stream outstr,int offset,int nx,int ny,int nz);
void gs_zyx(stream instr, stream outstr,int offset,int nx,int ny,int nz);
void gs_yzx(stream instr, stream outstr,int offset,int nx,int ny,int nz);
void gs_yxz(stream instr, stream outstr,int offset,int nx,int ny,int nz);

void Print2DArray(real *A, int nr, int nc);

nemo_main()
{
  stream instr, outstr;
  int i, j, i0, j0;
  int size[3], nsize;
  int offset[2];      // for now, ignore (i.e. assume 0,0)
  int n, nx,ny,nz;
  
  nsize = nemoinpi(getparam("size"),size,3);
  if (nsize < 1) error("problem specifying size=%s",getparam("size"));
  offset = getiparam("offset");
  nx = size[0];
  ny = size[1];
  nz = size[2];
  if (nz>1) error("z not yet supported");
  
  instr  = stropen(getparam("in"),"r");	 
  outstr = stropen(getparam("out"),"w");
  
  if (nsize==2) {
    gs_yx(instr,outstr,offset,nx,ny);
  }

}

void gs_xy(stream instr, stream outstr,int offset,int nx,int ny)
{
  /* no work needed */
}

void gs_yx(stream instr, stream outstr,int offset,int nx,int ny)
{
  int ix,iy;
  int nxy = nx*ny;

  


  for(ix=0; ix<nx; ix++) {
    for(iy=0; iy<ny; iy++) {
    }
  }
}

void gs_zxy(stream instr, stream outstr,int offset,int nx,int ny,int nz)
{
  int ix,iy,iz;
  
}

// A utility function to print a 2D array of size nr x nc and base address A
void Print2DArray(real *A, int nr, int nc)
{
    for(int r = 0; r < nr; r++)
    {
        for(int c = 0; c < nc; c++)
            printf("%f", *(A + r*nc + c));
 
        printf("\n");
    }
 
    printf("\n\n");
}
