/* 
 *	CCDGS:   play with gather-scatter operations on large blocks of data
 *
 *      gather:   loop(i,N) x[i] = y[idx[i]];
 *      scatter:  loop(i,N) y[idx[i]] = x[i];
 *
 *      24-jul-2013   Tinkertoy'd
 *      11-nov-2013   APT -> ATP cube;  goal is A=10  P=1e6  T=1e4 (=100GP)
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
  "bs=4\n         Block size per item",
  "mode=xy\n      Mode of copy",
  "offset=0\n     Offset (not used)",
  "VERSION=0.2a\n  26-jan-2021 PJT",
  NULL,
};

string usage="benchmark gather-scatter operations on large blocks of data";

void gs_xy(stream instr, stream outstr,int offset,int bs, int nx,int ny);
void gs_yx(stream instr, stream outstr,int offset,int bs, int nx,int ny);

void gs_xyz(stream instr, stream outstr,int offset,int bs, int nx,int ny,int nz);
void gs_yxz(stream instr, stream outstr,int offset,int bs, int nx,int ny,int nz);
void gs_zxy(stream instr, stream outstr,int offset,int bs, int nx,int ny,int nz);
void gs_zyx(stream instr, stream outstr,int offset,int bs, int nx,int ny,int nz);
void gs_yzx(stream instr, stream outstr,int offset,int bs, int nx,int ny,int nz);
void gs_yxz(stream instr, stream outstr,int offset,int bs, int nx,int ny,int nz);

void Print2DArray(real *A, int nr, int nc);

void nemo_main(void)
{
  stream instr, outstr;
  int i, j, i0, j0;
  int size[3], nsize;
  int offset, bs;
  int n, nx,ny,nz;
  char *mode;
  
  nsize = nemoinpi(getparam("size"),size,3);
  if (nsize < 1) error("problem specifying size=%s",getparam("size"));
  offset = getiparam("offset");
  bs = getiparam("bs");
  nx = size[0];
  ny = size[1];
  nz = (nsize < 3) ? 1 : size[2];
  if (nz>1) error("z=%d not yet supported",nz);
  mode = getparam("mode");
  
  instr  = stropen(getparam("in"),"r");	 
  outstr = stropen(getparam("out"),"w");
  
  if (streq(mode,"xy"))
    gs_xy(instr,outstr,offset,bs,nx,ny);
  else if (streq(mode,"yx"))
    gs_yx(instr,outstr,offset,bs,nx,ny);
  else
    warning("Mode %s not implemented",mode);

}

void gs_xy(stream instr, stream outstr,int offset,int bs, int nx,int ny)
{
  int nxy = nx*ny;
  int i, n;
  char *buf = (char *) allocate(sizeof(char)*bs);
  dprintf(0,"XY: straight copy %d blocks\n",nxy);
  for (i=0; i<nxy; i++) {
    n=fread(buf,bs,1,instr);
    if (n!=1) error("Not read at i=%d/%d",i,nxy);
    fwrite(buf,bs,1,outstr);
  }
}

void gs_yx(stream instr, stream outstr,int offset,int bs, int nx,int ny)
{
  int ix,iy,off,n;
  char *buf = (char *) allocate(sizeof(char)*bs);
  
  /* read sequential , write crazy */
  dprintf(0,"YX: sequential read, seed'd write (%d,%d) %d-blocks\n",nx,ny,bs);
  for(iy=0; iy<ny; iy++) {
    for(ix=0; ix<nx; ix++) {
      n=fread(buf,bs,1,instr);
      if (n!=1) error("Not read at ix,iy=%d,%d",ix,iy);
      //off=iy*nx+ix
      off=(ix*ny+iy)*bs;
      dprintf(1,"%d %d -> %d\n",ix,iy,off);
      fseek(outstr,off,SEEK_SET);
      fwrite(buf,bs,1,outstr);
    }
  }
}

void gs_zxy(stream instr, stream outstr,int offset,int bs,int nx,int ny,int nz)
{
  int ix,iy,iz;
  
}

// A utility function to print a 2D array of size nr x nc and base address A
void Print2DArray(real *A, int nr, int nc)
{
    int r, c;
    for(r = 0; r < nr; r++) {
      for(c = 0; c < nc; c++)
	printf("%f", *(A + r*nc + c));
      printf("\n");
    }
     printf("\n\n");
}
