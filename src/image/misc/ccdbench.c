/*
 * CCDBENCH: take moments in a cube along different axes
 *
 *       18-mar-2024    simple gather
 *                      
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>
#include <moment.h>
#include <mdarray.h>

string defv[] = {
  "in=???\n       Input image file",
  "order=xyz\n    Order to loop over cube, where the last axis is innermost loop",
  "iter=1\n       How often to loop",
  "mode=0\n       0=CubeValue(i,j,k)  1=cube[i][k][k]",
  "test=f\n       Show memory addresses for a small 2x3 array",
  "VERSION=0.3\n  13-apr-2024 PJT",
  NULL,
};

string usage = "benchmark taking moments along an axis of an image";

#define LOOP(i,n)     for(i=0;i<n;i++)
#define WHILE(iter)   while(iter--)

local void test();



void nemo_main()
{
  stream   instr;
  string   order = getparam("order");
  imageptr iptr = NULL;
  int      i,j,k,nx,ny,nz;
  int      iter = getiparam("iter");
  real     sum = 0.0;
  int      mode = getiparam("mode");

  if (getbparam("test")) {
    test();
    return;
  }
  
  instr = stropen(getparam("in"), "r");
    
  read_image( instr, &iptr);
  nx = Nx(iptr);	
  ny = Ny(iptr);
  nz = Nz(iptr);

  if (mode==1) {
    real ***cube = map3_image(iptr);
    dprintf(0,"cube[0][0][0]=%g\n",cube[0][0][0]);
    if (streq(order,"zyx"))
      WHILE(iter) LOOP(k,nz) LOOP(j,ny) LOOP(i,nx) sum += cube[i][j][k];
    else if (streq(order,"zxy"))
      WHILE(iter) LOOP(k,nz) LOOP(i,nx) LOOP(j,ny) sum += cube[i][j][k];
    else if (streq(order,"yxz"))
      WHILE(iter) LOOP(j,ny) LOOP(i,nx) LOOP(k,nz) sum += cube[i][j][k];
    else if (streq(order,"yzx"))
      WHILE(iter) LOOP(j,ny) LOOP(k,nz) LOOP(i,nx) sum += cube[i][j][k];
    else if (streq(order,"xyz"))
      WHILE(iter) LOOP(i,nx) LOOP(j,ny) LOOP(k,nz) sum += cube[i][j][k];
    else if (streq(order,"xzy"))
      WHILE(iter) LOOP(i,nx) LOOP(k,nz) LOOP(j,ny) sum += cube[i][j][k];
    else
      warning("Cannot do order %s",order);
  } else if (mode==0) {
    if (streq(order,"zyx"))
      WHILE(iter) LOOP(k,nz) LOOP(j,ny) LOOP(i,nx) sum += CubeValue(iptr,i,j,k);
    else if (streq(order,"zxy"))
      WHILE(iter) LOOP(k,nz) LOOP(i,nx) LOOP(j,ny) sum += CubeValue(iptr,i,j,k);    
    else if (streq(order,"yxz"))
      WHILE(iter) LOOP(j,ny) LOOP(i,nx) LOOP(k,nz) sum += CubeValue(iptr,i,j,k);    
    else if (streq(order,"yzx"))
      WHILE(iter) LOOP(j,ny) LOOP(k,nz) LOOP(i,nx) sum += CubeValue(iptr,i,j,k);    
    else if (streq(order,"xyz"))
      WHILE(iter) LOOP(i,nx) LOOP(j,ny) LOOP(k,nz) sum += CubeValue(iptr,i,j,k);    
    else if (streq(order,"xzy"))
      WHILE(iter) LOOP(i,nx) LOOP(k,nz) LOOP(j,ny) sum += CubeValue(iptr,i,j,k);    
    else
      warning("Cannot do order %s",order);
  } else
    warning("Cannot do mode %d",mode);
  printf("%s %d %s %g\n",order,getiparam("iter"),getparam("in"),sum);
}

// USE_IARRAY vs.   CDEF (row-major) or FORDEF (column-major)
// CDEF is the default, i.e.  data[row][col] , row after row is stored
local void test()
{
  static double is4[2][3];  // [ny][nx] - nx runs fastest
  printf("test:\n");
  printf("DD:  0x%x 0x%x  0x%x (0x%x)\n", &is4[0][0], &is4[0][1], &is4[1], is4[1]);

  imageptr iptr = NULL;
  create_image(&iptr, 3, 2);   // (nx,ny) - nx runs fastest
  printf("CV:  0x%x 0x%x  0x%x\n", &MapValue(iptr,0,0), &MapValue(iptr,0,1),  &MapValue(iptr,1,0));

  mdarray2 md2 = allocate_mdarray2(2,3);
  printf("MD:  0x%x 0x%x  0x%x (0x%x)\n", &md2[0][0], &md2[0][1], &md2[1], md2[1]);
}
