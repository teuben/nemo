/*
 * mdbench:  benchmarking 4D mdarray's
 *
 * 18-oct-2011   Created from mdarray testbed          Peter Teuben
 * 21-dec-2019   Tried out some OMP directives         PJT
 *
 */


#include <nemo.h>
#include <mdarray.h>

string defv[] = {
  "dim=10,20,30,40\n Dimensions of array A[dim1][dim2][dim3][dim4]....",
  "work=0\n          Work level (0=init  1=trans 2=reduc)",
  "ntest=1\n         How often to run the benchmark test",

  "flip=f\n          Reverse traversal through array, for benchmarking",
  "iter=1\n          Number of times to do the work, for benchmarking",
  "free=f\n          Free things we don't need anymore"
  "nprocs=-1\n       No OMP enabled",
  "VERSION=1.2\n     11-feb-2024 PJT",
  NULL,
};

string usage = "multidimensional array benchmarking (4D now)";

void init1(int *dim, mdarray1 x, bool flip)
{
  int i;

  for (i=0; i<dim[0]; i++)
    x[i] = i;
  
}

void work1(int *dim, mdarray1 x, bool flip)
{
  int i;
  real sum=0;

  for (i=0; i<dim[0]; i++)
    sum += x[i]; 
  dprintf(1,"sum1=%g\n",sum);
}


void init2(int *dim, mdarray2 x, bool flip)
{
  int i,j;

  if (flip) {
    for (i=0; i<dim[0]; i++)
      for (j=0; j<dim[1]; j++)
	x[j][i] = (real)i+j*10;
  } else {
    for (j=0; j<dim[1]; j++)
      for (i=0; i<dim[0]; i++)
	x[j][i] = (real)i+j*10;
  }
}

void work2(int *dim, mdarray2 x, bool flip)
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

void init3(int *dim, mdarray3 x, bool flip)
{
  int i,j,k;


  if (flip) {
    for (i=0; i<dim[0]; i++)
      for (j=0; j<dim[1]; j++)
	for (k=0; k<dim[2]; k++)
	  x[k][j][i] = (real)i+j*10+k*100;
  } else {
    for (k=0; k<dim[2]; k++)
      for (j=0; j<dim[1]; j++)
	for (i=0; i<dim[0]; i++)
	  x[k][j][i] = (real)i+j*10+k*100;
  }
}

void work3(int *dim, mdarray3 x, bool flip)
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

void sinit3(int *dim, real x[][4][4], bool flip)
{
  int i,j,k;

  if (flip) {
    for (i=0; i<dim[0]; i++)
      for (j=0; j<dim[1]; j++)
	for (k=0; k<dim[2]; k++)
	  x[k][j][i] = (real)i+j*10+k*100;
  } else {
    for (k=0; k<dim[2]; k++)
      for (j=0; j<dim[1]; j++)
	for (i=0; i<dim[0]; i++)
	  x[k][j][i] = (real)i+j*10+k*100;
  }
}

void swork3(int *dim, real x[][4][4], bool flip)
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

void init4(int *dim, mdarray4 x, bool flip)
{
  int i,j,k,l;

  if (flip) {
    for (i=0; i<dim[0]; i++)
      for (j=0; j<dim[1]; j++)
	for (k=0; k<dim[2]; k++)
	  for (l=0; l<dim[3]; l++)
	    x[l][k][j][i] = (real)i+j*10+k*100+l*1000;
  } else {
    for (l=0; l<dim[3]; l++)
      for (k=0; k<dim[2]; k++)
	for (j=0; j<dim[1]; j++)
	  for (i=0; i<dim[0]; i++)
	    x[l][k][j][i] = (real)i+j*10+k*100+l*1000;
  }
}

void work4(int *dim, mdarray4 x, bool flip)
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

void work5(int *dim, mdarray5 x, bool flip)
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

void work6(int *dim, mdarray6 x, bool flip)
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

void work7(int *dim, mdarray7 x, bool flip)
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

void nemo_main()
{
  int i,dim[MDMAXDIM];
  int ndim  = nemoinpi(getparam("dim"),dim,MDMAXDIM);
  bool flip = getbparam("flip");
  bool free = getbparam("free");
  int iter  = getiparam("iter");
  int ntest = getiparam("ntest");
  int iwork = getiparam("work");
  int nprocs= getiparam("nprocs");
#if 1
  int test1 = 20, test2=10, test3[test2][test1];   // in C99 this is now allowed   +x50
#else
  int test1 = 20, test2=10, test3[10][20];         // old K&R style
#endif
  real sum;
  int i1,i2,i3,i4;

  mdarray1 x1;
  mdarray2 x2;
  mdarray3 x3, y3;
  mdarray4 x4, y4;
  mdarray5 x5;
  mdarray6 x6;
  mdarray7 x7;

  dprintf(0,"Using single CPU, no OMP enables\n");
  if (nprocs>1) warning("No OMP was enabled");

  /* C99 now does it the way I wanted it to work */
  dprintf(1,"pointer test3: 0x%x 0x%x 0x%x    0x%x   0x%x 0x%x\n",test3,test3[0],&test3[0][0],&test3[0][1],test3[1],&test3[1][0]);
  // @todo test if address differences are the right amount of sizeof(int) and test1*sizeof(int)
  void *p1 = test3;
  void *p2 = test3[0];
  void *p3 = &test3[0][0];
  if (p1 != p2) warning("p1!=p2");
  if (p1 != p3) warning("p1!=p3");
  void *p4 = &test3[0][1];
  if (p4-p3 != sizeof(int)) warning("p4 != p3+%d",sizeof(int));
  void *p5 = test3[1];  
  void *p6 = &test3[1][0];
  if (p5 != p6) warning("p5!=p6");
  if (p4-p3 != sizeof(int)) warning("p4 != p3+%d",sizeof(int));
  if (p5-p1 != test1*sizeof(int)) warning("p5 != p1+%d",test1*sizeof(int));

  if (ndim != 4) error("ndim=4 for now");

  x4 = allocate_mdarray4(dim[3],dim[2],dim[1],dim[0]);
#if 0
  warning("4DIM test wrong");
  //  80^4 * 10 -> 20.3"
  //  90^4 * 1  ->  4.1
  // 100^4 * 1  ->  8.2
  for (i=0; i<ntest; i++) {
    for (i1=0; i1<dim[0]; i1++)
      for (i2=0; i2<dim[1]; i2++)
	for (i3=0; i3<dim[2]; i3++)
	  for (i4=0; i4<dim[3]; i4++) {
	    x4[i4][i3][i2][i1] = i4 + 100*(i3+100*(i2 + 100*i1));
	  }
  }
#else
  warning("4DIM test right");
  //  80^4 * 10 -> 3.832"
  //  90^4 * 1  ->  0.63
  // 100^4 * 1  ->  1.00
  for (i=0; i<ntest; i++) {
    for (i4=0; i4<dim[3]; i4++) 
      for (i3=0; i3<dim[2]; i3++)
	for (i2=0; i2<dim[1]; i2++)
	  for (i1=0; i1<dim[0]; i1++) {
	    x4[i4][i3][i2][i1] = i4 + 100*(i3+100*(i2 + 100*i1));
	  }
  }
#endif
  if (iwork>0) {
#if 1
    warning("4DIM  x4->y4->y3");
    // transpose a hypercube
    // 80^4:  0.38   1.08
    // 90^4:  0.62   1.60
    // 95^4:  0.74   2.00 
    y4 = allocate_mdarray4(dim[3],dim[1],dim[0],dim[2]);
    for (i=0; i<ntest; i++) {
      for (i4=0; i4<dim[3]; i4++) 
	for (i3=0; i3<dim[2]; i3++)
	  for (i2=0; i2<dim[1]; i2++)
	    for (i1=0; i1<dim[0]; i1++) {
	      y4[i4][i2][i1][i3] = x4[i4][i3][i2][i1];
	    }//i1
    }//i
    if (free) free_mdarray4(x4,dim[3],dim[2],dim[1],dim[0]);
    if (iwork>1) {
      y3 = allocate_mdarray3(dim[3],dim[1],dim[0]);
      for (i=0; i<ntest; i++) {
	for (i4=0; i4<dim[3]; i4++) 
	  for (i2=0; i2<dim[1]; i2++)
	    for (i1=0; i1<dim[0]; i1++) {
	      sum = 0.0;
	      for (i3=0; i3<dim[2]; i3++) {
		sum += y4[i4][i2][i1][i3];
	      }//i3
	      y3[i4][i2][i1] = sum;
	    }//i1
      }//i
      // report
      sum = 0.0;
      for (i4=0; i4<dim[3]; i4++) 
	for (i2=0; i2<dim[1]; i2++)
	  for (i1=0; i1<dim[0]; i1++)
	    sum += y3[i4][i2][i1];
      dprintf(0,"y3-sum=%g\n",sum);
    }//iwork>1
#else
    // instead of x4 -> y4 -> y3
    // this reduced x4->y3
    // 90^4 now goes in 1.5" instead of 2.2" in x4-y4-y3 chain
    warning("4DIM  x4->y3");
    y3 = allocate_mdarray3(dim[3],dim[1],dim[0]);
    for (i=0; i<ntest; i++) {
      for (i4=0; i4<dim[3]; i4++) 
	for (i2=0; i2<dim[1]; i2++)
	  for (i1=0; i1<dim[0]; i1++) {
	    sum = 0.0;
	    for (i3=0; i3<dim[2]; i3++) {
	      sum += x4[i4][i3][i2][i1];
	    }//i3
	    y3[i4][i2][i1] = sum;
	  }//i1
    }//i
#endif
  }//iwork>0
}

