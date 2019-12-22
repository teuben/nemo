/*
 * mdarray: a simple multi-dimensional array (allocator)
 *          cannot be used in a multi-threaded environment
 *
 * See also: Iliffe Vector
 *           Shortridge ADASS 2019 proceedings - https://github.com/KnaveAndVarlet/ADASS2019
 *
 * Note:  MDMAXDIM defines the highest dimension we've implemented (currently 7)
 *        see mdarray.h
 *
 *  5-may-2003   Created                                 Peter Teuben
 * 18-feb-2006   use sequential memory for the data      PJT
 *
 */

#include <nemo.h>
#include <mdarray.h>

static int top = 0;    /* highest mdarray will have this at 0, each call increases top */
static real *p = 0;    /* keep pointer to the big alloc incrementing by n1 */


mdarray1 allocate_mdarray1(int n1)
{
  real *p1;
  if (top) {
    dprintf(7,"  top @ 0x%x\n",p);
    p1 = p;
    p += n1;
    return p1;
  } else
    return (mdarray1 )allocate(n1*sizeof(real));
}

mdarray2 allocate_mdarray2(int n2,int n1)
{
  mdarray2 x = (mdarray2)allocate(sizeof(mdarray2)*n2);
  int i;
  if (top == 0) 
    p = (real *) allocate(sizeof(real)*n1*n2);
  top++;
  for (i=0;i<n2;i++)
    x[i] = allocate_mdarray1(n1);
  top--;
  return x;
}

mdarray3 allocate_mdarray3(int n3,int n2,int n1)
{
  mdarray3 x = (mdarray3)allocate(sizeof(mdarray3)*n3);
  int i;
  if (top == 0)
    p = (real *) allocate(sizeof(real)*n1*n2*n3);
  top++;
  for (i=0;i<n3;i++)
    x[i] = allocate_mdarray2(n2,n1);
  top--;
  return x;
}

mdarray4 allocate_mdarray4(int n4,int n3,int n2,int n1)
{
  mdarray4 x = (mdarray4)allocate(sizeof(mdarray4)*n4);
  int i;
  if (top == 0)
    p = (real *) allocate(sizeof(real)*n1*n2*n3*n4);
  top++;
  for (i=0;i<n4;i++)
    x[i] = allocate_mdarray3(n3,n2,n1);
  top--;
  return x;
}

mdarray5 allocate_mdarray5(int n5,int n4,int n3,int n2,int n1)
{
  mdarray5 x = (mdarray5)allocate(sizeof(mdarray5)*n5);
  int i;
  if (top == 0)
    p = (real *) allocate(sizeof(real)*n1*n2*n3*n4*n5);
  top++;
  for (i=0;i<n5;i++)
    x[i] = allocate_mdarray4(n4,n3,n2,n1);
  top--;
  return x;
}

mdarray6 allocate_mdarray6(int n6,int n5,int n4,int n3,int n2,int n1)
{
  mdarray6 x = (mdarray6)allocate(sizeof(mdarray6)*n6);
  int i;
  if (top == 0)
    p = (real *) allocate(sizeof(real)*n1*n2*n3*n4*n5*n6);
  top++;
  for (i=0;i<n6;i++)
    x[i] = allocate_mdarray5(n5,n4,n3,n2,n1);
  top--;
  return x;
}

mdarray7 allocate_mdarray7(int n7,int n6,int n5,int n4,int n3,int n2,int n1)
{
  mdarray7 x = (mdarray7)allocate(sizeof(mdarray7)*n7);
  int i;
  if (top == 0)
    p = (real *) allocate(sizeof(real)*n1*n2*n3*n4*n5*n6*n7);
  top++;
  for (i=0;i<n7;i++)
    x[i] = allocate_mdarray6(n6,n5,n4,n3,n2,n1);
  top--;
  return x;
}

mdarray8 allocate_mdarray8(int n8,int n7,int n6,int n5,int n4,int n3,int n2,int n1)
{
  mdarray8 x = (mdarray8)allocate(sizeof(mdarray8)*n8);
  int i;
  if (top == 0)
    p = (real *) allocate(sizeof(real)*n1*n2*n3*n4*n5*n6*n7*n8);
  top++;
  for (i=0;i<n8;i++)
    x[i] = allocate_mdarray7(n7,n6,n5,n4,n3,n2,n1);
  top--;
  return x;
}

void free_mdarray1(mdarray1 x, int n1)
{
  if (top==0) free(x);
}

void free_mdarray2(mdarray2 x, int n2, int n1)
{
  int i;
  if (top==0) free(x[0]);
  top++;
  for (i=0;i<n2;i++)
    free_mdarray1(x[i],n1);
  top--;
  free(x);
}

void free_mdarray3(mdarray3 x, int n3, int n2, int n1)
{
  int i;
  if (top==0) free(x[0][0]);
  top++;
  for (i=0;i<n3;i++)
    free_mdarray2(x[i],n2,n1); 
  top--;
  free(x);
}

void free_mdarray4(mdarray4 x, int n4, int n3, int n2, int n1)
{
  int i;
  if (top==0) free(x[0][0][0]);
  top++;
  for (i=0;i<n4;i++)
    free_mdarray3(x[i],n3,n2,n1);
  top--;
  free(x);
}

void free_mdarray5(mdarray5 x, int n5, int n4, int n3, int n2, int n1)
{
  int i;
  if (top==0) free(x[0][0][0][0]);
  top++;
  for (i=0;i<n5;i++)
    free_mdarray4(x[i],n4,n3,n2,n1);
  top--;
  free(x);
}

void free_mdarray6(mdarray6 x, int n6, int n5, int n4, int n3, int n2, int n1)
{
  int i;
  if (top==0) free(x[0][0][0][0][0]);
  top++;
  for (i=0;i<n6;i++)
    free_mdarray5(x[i],n5,n4,n3,n2,n1);
  top--;
  free(x);
}

void free_mdarray7(mdarray7 x, int n7, int n6, int n5, int n4, int n3, int n2, int n1)
{
  int i;
  if (top==0) free(x[0][0][0][0][0][0]);
  top++;
  for (i=0;i<n7;i++)
    free_mdarray6(x[i],n6,n5,n4,n3,n2,n1);
  top--;
  free(x);
}

void free_mdarray8(mdarray8 x, int n8, int n7, int n6, int n5, int n4, int n3, int n2, int n1)
{
  int i;
  if (top==0) free(x[0][0][0][0][0][0][0]);
  top++;
  for (i=0;i<n8;i++)
    free_mdarray7(x[i],n7,n6,n5,n4,n3,n2,n1);
  top--;
  free(x);
}



/* =============================================================================== */
#ifdef TESTBED

string defv[] = {
  "dim=10,20,5\n    Dimensions of array A[dim1][dim2][dim3]....",
  "flip=f\n         Reverse traversal through array, for benchmarking",
  "iter=1\n         Number of times to do the work, for benchmarking",
  "free=t\n         Free-Allocate the array each iter?",
  "layout=f\n       Show memory layout",
  "static=f\n       Benchmark on a static array",
  "transpose=f\n    Transpose a 2D array",
  "axb=f\n          Test an Ax=b",
  "VERSION=2.4\n    13-dec-2019 PJT",
  NULL,
};

string usage = "testing nested multidimensional arrays in C";

init1(int *dim, mdarray1 x, bool flip)
{
  int i;

  for (i=0; i<dim[0]; i++)
    x[i] = i;
  
}

work1(int *dim, mdarray1 x, bool flip)
{
  int i;
  real sum=0;

  for (i=0; i<dim[0]; i++)
    sum += x[i]; 
  dprintf(1,"sum1=%g\n",sum);
}


init2(int *dim, mdarray2 x,bool flip)
{
  int i,j;

  if (flip) {
    for (i=0; i<dim[0]; i++)
      for (j=0; j<dim[1]; j++)
	x[j][i] = (real)i+j;
  } else {
    for (j=0; j<dim[1]; j++)
      for (i=0; i<dim[0]; i++)
	x[j][i] = (real)i+j;
  }
}

work2(int *dim, mdarray2 x,bool flip)
{
  int i,j;
  real sum=0;

  if (flip) {
    for (i=0; i<dim[0]; i++)
      for (j=0; j<dim[1]; j++)
	sum += x[j][i];
  } else {
    for (j=0; j<dim[1]; j++)
      for (i=0; i<dim[0]; i++)
	sum += x[j][i];
  }
  dprintf(1,"sum2=%g\n",sum);
}

show2(int *dim, mdarray2 x)
{
  int i,j, nsize = 1;
  real *p;
  for (i=0; i<2; i++)  nsize *= dim[i];
  
  for (j=0; j<dim[1]; j++)
    for (i=0; i<dim[0]; i++)
      printf("%g ",x[j][i]);
  printf("\n");

  p = &x[0][0];
  for (i=0; i<nsize; i++)
    printf("%g ",p[i]);
  printf("\n");
}


init3(int *dim, mdarray3 x, bool flip)
{
  int i,j,k;


  if (flip) {
    for (i=0; i<dim[0]; i++)
      for (j=0; j<dim[1]; j++)
	for (k=0; k<dim[2]; k++)
	  x[k][j][i] = (real)i+j+k;
  } else {
    for (k=0; k<dim[2]; k++)
      for (j=0; j<dim[1]; j++)
	for (i=0; i<dim[0]; i++)
	  x[k][j][i] = (real)i+j+k;
  }
}

work3(int *dim, mdarray3 x, bool flip)
{
  int i,j,k;
  real sum=0;

  if (flip) {
    for (i=0; i<dim[0]; i++)
      for (j=0; j<dim[1]; j++)
	for (k=0; k<dim[2]; k++)
	  sum += x[k][j][i];	  
  } else {
    for (k=0; k<dim[2]; k++)
      for (j=0; j<dim[1]; j++)
	for (i=0; i<dim[0]; i++)
	  sum += x[k][j][i];	  
  }
  dprintf(1,"sum3=%g\n",sum);
}

show3(int *dim, mdarray3 x)
{
  int i,j,k, nsize = 1;
  real *p;
  for (i=0; i<3; i++)  nsize *= dim[i];
  
  for (k=0; k<dim[2]; k++)
    for (j=0; j<dim[1]; j++)
      for (i=0; i<dim[0]; i++)
	printf("%g ",x[k][j][i]);
  printf("\n");

  p = &x[0][0][0];
  for (i=0; i<nsize; i++)
    printf("%g ",p[i]);
  printf("\n");
}


sinit3(int *dim, real x[][4][4], bool flip)
{
  int i,j,k;

  if (flip) {
    for (i=0; i<dim[0]; i++)
      for (j=0; j<dim[1]; j++)
	for (k=0; k<dim[2]; k++)
	  x[k][j][i] = (real)i+j+k;
  } else {
    for (k=0; k<dim[2]; k++)
      for (j=0; j<dim[1]; j++)
	for (i=0; i<dim[0]; i++)
	  x[k][j][i] = (real)i+j+k;
  }
}

swork3(int *dim, real x[][4][4], bool flip)
{
  int i,j,k;
  real sum=0;

  if (flip) {
    for (i=0; i<dim[0]; i++)
      for (j=0; j<dim[1]; j++)
	for (k=0; k<dim[2]; k++)
	  sum += x[k][j][i];	  
  } else {
    for (k=0; k<dim[2]; k++)
      for (j=0; j<dim[1]; j++)
	for (i=0; i<dim[0]; i++)
	  sum += x[k][j][i];	  
  }
  dprintf(1,"sum3=%g\n",sum);
}

init4(int *dim, mdarray4 x, bool flip)
{
  int i,j,k,l;

  if (flip) {
    for (i=0; i<dim[0]; i++)
      for (j=0; j<dim[1]; j++)
	for (k=0; k<dim[2]; k++)
	  for (l=0; l<dim[3]; l++)
	    x[l][k][j][i] = (real)i+j+k+l;
  } else {
    for (l=0; l<dim[3]; l++)
      for (k=0; k<dim[2]; k++)
	for (j=0; j<dim[1]; j++)
	  for (i=0; i<dim[0]; i++)
	    x[l][k][j][i] = (real)i+j+k+l;
  }
}

work4(int *dim, mdarray4 x, bool flip)
{
  int i,j,k,l;
  real sum=0;

  if (flip) {
    for (i=0; i<dim[0]; i++)
      for (j=0; j<dim[1]; j++)
	for (k=0; k<dim[2]; k++)
	  for (l=0; l<dim[3]; l++)
	    sum += x[l][k][j][i];
  } else {
    for (l=0; l<dim[3]; l++)
      for (k=0; k<dim[2]; k++)
	for (j=0; j<dim[1]; j++)
	  for (i=0; i<dim[0]; i++)
	    sum += x[l][k][j][i];
  }
  dprintf(1,"sum4=%g\n",sum);
}

work5(int *dim, mdarray5 x, bool flip)
{
  int i,j,k,l,m;

  if (flip) {
    for (i=0; i<dim[0]; i++)
      for (j=0; j<dim[1]; j++)
	for (k=0; k<dim[2]; k++)
	  for (l=0; l<dim[3]; l++)
	    for (m=0; m<dim[4]; m++)
	      x[m][l][k][j][i] = i+j+k+l+m;
  } else {
    for (m=0; m<dim[4]; m++)
      for (l=0; l<dim[3]; l++)
	for (k=0; k<dim[2]; k++)
	  for (j=0; j<dim[1]; j++)
	    for (i=0; i<dim[0]; i++)
	      x[m][l][k][j][i] = i+j+k+l+m;
  }
}

work6(int *dim, mdarray6 x, bool flip)
{
  int i,j,k,l,m,n;

  if (flip) {
    for (i=0; i<dim[0]; i++)
      for (j=0; j<dim[1]; j++)
	for (k=0; k<dim[2]; k++)
	  for (l=0; l<dim[3]; l++)
	    for (m=0; m<dim[4]; m++)
	      for (n=0; n<dim[5]; n++)
		x[n][m][l][k][j][i] = i+j+k+l+m+n;
  } else {
    for (n=0; n<dim[5]; n++)
      for (m=0; m<dim[4]; m++)
	for (l=0; l<dim[3]; l++)
	  for (k=0; k<dim[2]; k++)
	    for (j=0; j<dim[1]; j++)
	      for (i=0; i<dim[0]; i++)
		x[n][m][l][k][j][i] = i+j+k+l+m+n;
  }
}

work7(int *dim, mdarray7 x, bool flip)
{
  int i,j,k,l,m,n,o;

  if (flip) {
    for (i=0; i<dim[0]; i++)
      for (j=0; j<dim[1]; j++)
	for (k=0; k<dim[2]; k++)
	  for (l=0; l<dim[3]; l++)
	    for (m=0; m<dim[4]; m++)
	      for (n=0; n<dim[5]; n++)
		for (o=0; o<dim[6]; o++)
		  x[o][n][m][l][k][j][i] = i+j+k+l+m+n+o;
  } else {
    for (o=0; o<dim[6]; o++)
      for (n=0; n<dim[5]; n++)
	for (m=0; m<dim[4]; m++)
	  for (l=0; l<dim[3]; l++)
	    for (k=0; k<dim[2]; k++)
	      for (j=0; j<dim[1]; j++)
		for (i=0; i<dim[0]; i++)
		  x[o][n][m][l][k][j][i] = i+j+k+l+m+n+o;
  }
}

static real p1 = 1;
static real p2 = 2;
static real p3[2][2][2] = {11,12,13,14,15,16,17,18};
static real p4 = 20;
static real p5[2][2] = {101,102,103,104};
static real p6 = 201;

void try()
{
  real *t;
  t = &p1;

  printf("static C arrays:\n");
  printf("p1:\n");
  printf("0x%x\n",&p1);
  printf("p2:\n");
  printf("0x%x\n",&p2);
  printf("p3:\n");
  printf("0x%x\n",p3);
  printf("0x%x\n",p3[0]);
  printf("0x%x\n",p3[0][0]);
  printf("0x%x\n",&p3[0][0][0]);
  printf("p4:\n");
  printf("0x%x\n",&p4);
  printf("p5:\n");
  printf("0x%x\n",p5);
  printf("0x%x\n",p5[0]);
  printf("0x%x\n",&p5[0][0]);
  printf("p6:\n");
  printf("0x%x\n",&p6);

  printf("t=%g\n",*t++);
  printf("t=%g\n",*t++);
  printf("t=%g\n",*t++);
  printf("t=%g\n",*t++);
  printf("t=%g\n",*t++);
  printf("t=%g\n",*t++);
}

#define MAXDIM 100

nemo_main()
{
  int i,j,dim[MAXDIM];
  int ndim = nemoinpi(getparam("dim"),dim,MAXDIM);
  bool flip = getbparam("flip"), free=getbparam("free"), trans=getbparam("transpose");
  bool axb = getbparam("axb");
  int iter=getiparam("iter");
  bool layout = getbparam("layout");
  bool statbench = getbparam("static");
  real a2[2][2];
  real a3[4][4][4];
  real a7[3][3][3][3][3][3][3];
  real *p;
  int i1,i2,i3,i4, nd;
  int nsize;

  mdarray1 x1;
  mdarray2 x2, x2t;
  mdarray3 x3;
  mdarray4 x4;
  mdarray5 x5;
  mdarray6 x6;
  mdarray7 x7;

#if 1
  mdarray2  A;
  mdarray1  x,b;

  if (axb) {
    A  = allocate_mdarray2(2,3);       /*  Ax=b ;   A[2][3] */
    x  = allocate_mdarray1(3);         /*           x[3]    */
    b  = allocate_mdarray1(2);         /*           b[2]    */
    A[0][0] = 1;
    A[0][1] = 2;
    A[0][2] = 3;
    A[1][0] = 4;
    A[1][1] = 5;
    A[1][2] = 6;
    x[0] = 3;
    x[1] = 2;
    x[2] = 1;
    for (i=0; i<2; i++)
      for (j=0, b[i]=0; j<3; j++) b[i] += A[i][j] * x[j];
    for (i=0; i<2; i++)
      printf("b[%d] = %g\n",i,b[i]);
	
    stop(0);
  }
#else
#endif

  if (trans) {
    if (ndim != 2) error("transposing needs 2 dimensional array");
    x2  = allocate_mdarray2(dim[1],dim[0]);
    x2t = allocate_mdarray2(dim[0],dim[1]);
    while (iter-- > 0) {
      init2(dim,x2,flip);
      for (i1=0; i1<nd; i1++)
	for (i2=0; i2<nd; i2++) {
	  x2t[i2][i1] = x2[i1][i2];
      }
    }
    if (free) {
      free_mdarray2(x2, dim[1],dim[0]);
      free_mdarray2(x2t,dim[0],dim[1]);
    }
    stop(0);
  }

  if (layout) {
    try();
    printf("a2 static\n");
    printf("0x%x\n",a2);
    printf("0x%x\n",a2[0]);
    printf("0x%x\n",&a2[0][0]);
    for (i1=0; i1<2; i1++)
      for (i2=0; i2<2; i2++)
	printf("  0x%x\n",&a2[i1][i2]);

    nd = 2;

    printf("x2 dynamic %d\n",nd);
    x2 = allocate_mdarray2(nd,nd);
    printf("0x%x\n",x2);
    printf("0x%x\n",x2[0]);
    printf("0x%x\n",&x2[0][0]);
    p = &x2[0][0];
    for (i1=0; i1<nd; i1++)
      for (i2=0; i2<nd; i2++) {
	printf("  0x%x  0x%d\n",&x2[i1][i2],&x2[i1][i2]-p);
	p = &x2[i1][i2];
      }
    
    printf("x3 dynamic %d\n",nd);
    x3 = allocate_mdarray3(nd,nd,nd);
    printf("0x%x\n",x3);
    printf("0x%x\n",x3[0]);
    printf("0x%x\n",x3[0][0]);
    printf("0x%x\n",&x3[0][0][0]);
    p = &x3[0][0][0];
    for (i1=0; i1<nd; i1++)
      for (i2=0; i2<nd; i2++)
	for (i3=0; i3<nd; i3++) {
	  printf("  0x%x 0x%x\n",&x3[i1][i2][i3],&x3[i1][i2][i3]-p);
	  p = &x3[i1][i2][i3];
	}

    printf("x4 dynamic %d\n",nd);
    x4 = allocate_mdarray4(nd,nd,nd,nd);
    printf("0x%x\n",x4);
    printf("0x%x\n",x4[0]);
    printf("0x%x\n",x4[0][0]);
    printf("0x%x\n",x4[0][0][0]);
    printf("0x%x\n",&x4[0][0][0][0]);
    p = &x4[0][0][0][0];
    for (i1=0; i1<nd; i1++)
      for (i2=0; i2<nd; i2++)
	for (i3=0; i3<nd; i3++)
	  for (i4=0; i4<nd; i4++) {
	    printf("  0x%x 0x%x\n",&x4[i1][i2][i3][i4],&x4[i1][i2][i3][i4]-p);
	    p = &x4[i1][i2][i3][i4];
	  }

    stop(0); 
  }

  printf("Working with ndim=%d MDArray",ndim);
  nsize = 1;
  for (i=ndim-1; i>=0; i--) {
    printf("[%d]",dim[i]);
    nsize *= dim[i];
  }
  printf("  Total size %d\n",nsize);

  if (ndim==1) {
    x1 = allocate_mdarray1(dim[0]);
    while (iter-- > 0) {
      init1(dim,x1,flip);
      work1(dim,x1,flip);
    }
    if (free) free_mdarray1(x1,dim[0]);
  } else if (ndim==2) {
    x2 = allocate_mdarray2(dim[1],dim[0]);
    while (iter-- > 0) {
      init2(dim,x2,flip);
      work2(dim,x2,flip);
    }
    show2(dim,x2);
    if (free) free_mdarray2(x2,dim[1],dim[0]);
  } else if (ndim==3) {
    if (statbench) {
      warning("dim=4,4,4 is needed for this static test");
      while(iter-- > 0) {
	sinit3(dim,a3,flip);
	swork3(dim,a3,flip);
      }
    } else {
      x3 = allocate_mdarray3(dim[2],dim[1],dim[0]);
      while (iter-- > 0) {
	init3(dim,x3,flip);
	work3(dim,x3,flip);
      }
      show3(dim,x3);      
      if (free) free_mdarray3(x3,dim[2],dim[1],dim[0]);
    }
  } else if (ndim==4) {
    x4 = allocate_mdarray4(dim[3],dim[2],dim[1],dim[0]);
    while (iter-- > 0) {
      init4(dim,x4,flip);
      work4(dim,x4,flip);
    }
    if (free) free_mdarray4(x4,dim[3],dim[2],dim[1],dim[0]);
  } else if (ndim==5) {
    x5 = allocate_mdarray5(dim[4],dim[3],dim[2],dim[1],dim[0]);
    work5(dim,x5,flip);
    if (free) free_mdarray5(x5,dim[4],dim[3],dim[2],dim[1],dim[0]);
  } else if (ndim==6) {
    x6 = allocate_mdarray6(dim[5],dim[4],dim[3],dim[2],dim[1],dim[0]);
    work6(dim,x6,flip);
    if (free) free_mdarray6(x6,dim[5],dim[4],dim[3],dim[2],dim[1],dim[0]);
  } else if (ndim==7) {
    x7 = allocate_mdarray7(dim[6],dim[5],dim[4],dim[3],dim[2],dim[1],dim[0]);
    work7(dim,x7,flip);
    if (free) free_mdarray7(x7,dim[6],dim[5],dim[4],dim[3],dim[2],dim[1],dim[0]);
  } else if (ndim==8) {
    warning("not written, but it's supported");
  } else
    error("dimension %d not supported yet (MDMAXDIM=%d)",ndim,MDMAXDIM);
}


#endif
