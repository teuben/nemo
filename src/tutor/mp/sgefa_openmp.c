# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <omp.h>
# include <time.h>

int main ( void );
void test01 ( int n );
void test02 ( int n );
void test03 ( int n );

int isamax ( int n, float x[], int incx );
void matgen ( int lda, int n, float a[], float x[], float b[] );
void msaxpy ( int nr, int nc, float a[], int n, float x[], float y[] );
void msaxpy2 ( int nr, int nc, float a[], int n, float x[], float y[] );
int msgefa ( float a[], int lda, int n, int ipvt[] );
int msgefa2 ( float a[], int lda, int n, int ipvt[] );
void saxpy ( int n, float a, float x[], int incx, float y[], int incy );
float sdot ( int n, float x[], int incx, float y[], int incy );
int sgefa ( float a[], int lda, int n, int ipvt[] );
void sgesl ( float a[], int lda, int n, int ipvt[], float b[], int job );
void sscal ( int n, float a, float x[], int incx );
void sswap ( int n, float x[], int incx, float y[], int incy );
void timestamp ( );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for the SGEFA_OPENMP test program.

  Discussion:

    We want to compare methods of solving the linear system A*x=b.

    The first way uses the standard sequential algorithm "SGEFA".

    The second way uses a variant of SGEFA that has been modified to
    take advantage of OpenMP.

    The third way reruns the variant code, but with OpenMP turned off.

  Modified:

    17 April 2009

  Author:

    John Burkardt
*/
{
  int n;

  timestamp ( );

  printf ( "\n" );
  printf ( "SGEFA_OPENMP\n" );
  printf ( "  C + OpenMP version\n" );

  printf ( "\n" );
  printf ( "  Number of processors available = %d\n", omp_get_num_procs ( ) );
  printf ( "  Number of threads =              %d\n", omp_get_max_threads ( ) );

  printf ( "\n" );
  printf ( " Algorithm        Mode          N    Error       Time\n" );

  printf ( "\n" );
  n = 10;
  test01 ( n );
  test02 ( n );
  test03 ( n );

  printf ( "\n" );
  n = 100;
  test01 ( n );
  test02 ( n );
  test03 ( n );

  printf ( "\n" );
  n = 1000;
  test01 ( n );
  test02 ( n );
  test03 ( n );

  printf ( "\n" );
  printf ( "SGEFA_OPENMP\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( int n )

/******************************************************************************/
/*
  Purpose:

    TEST01 runs the sequential version of SGEFA.

  Modified:

    07 April 2008

  Author:

    John Burkardt
*/
{
  float *a;
  float *b;
  float err;
  int i;
  int info;
  int *ipvt;
  int job;
  int lda;
  double wtime;
  float *x;
/*
  Generate the linear system A * x = b.
*/
  lda = n;
  a = ( float * ) malloc ( lda * n * sizeof ( float ) );
  b = ( float * ) malloc ( n * sizeof ( float ) );
  x = ( float * ) malloc ( n * sizeof ( float ) );

  matgen ( lda, n, a, x, b );
/*
  Factor the linear system.
*/
  ipvt = ( int * ) malloc ( n * sizeof ( int ) );

  wtime = omp_get_wtime ( );
  info = sgefa ( a, lda, n, ipvt );
  wtime = omp_get_wtime ( ) - wtime;

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST01 - Fatal error!\n" );
    printf ( "  SGEFA reports the matrix is singular.\n" );
    exit ( 1 );
  }
/*
  Solve the linear system.
*/
  job = 0;
  sgesl ( a, lda, n, ipvt, b, job );

  err = 0.0;
  for ( i = 0; i < n; i++ )
  {
    err = err + fabs ( x[i] - b[i] );
  }
  printf ( "  Original  Sequential   %8d  %10.4e  %10.4e\n", n, err, wtime );

  free ( a );
  free ( b );
  free ( ipvt );
  free ( x );

  return;
}
/******************************************************************************/

void test02 ( int n )

/******************************************************************************/
/*
  Purpose:

    TEST02 runs the revised version of SGEFA in parallel.

  Modified:

    07 April 2008

  Author:

    John Burkardt
*/
{
  float *a;
  float *b;
  float err;
  int i;
  int info;
  int *ipvt;
  int job;
  int lda;
  double wtime;
  float *x;
/*
  Generate the linear system A * x = b.
*/
  lda = n;
  a = ( float * ) malloc ( lda * n * sizeof ( float ) );
  b = ( float * ) malloc ( n * sizeof ( float ) );
  x = ( float * ) malloc ( n * sizeof ( float ) );

  matgen ( lda, n, a, x, b );
/*
  Factor the linear system.
*/
  ipvt = ( int * ) malloc ( n * sizeof ( int ) );

  wtime = omp_get_wtime ( );
  info = msgefa ( a, lda, n, ipvt );
  wtime = omp_get_wtime ( ) - wtime;

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST02 - Fatal error!\n" );
    printf ( "  MSGEFA reports the matrix is singular.\n" );
    exit ( 1 );
  }
/*
  Solve the linear system.
*/
  job = 0;
  sgesl ( a, lda, n, ipvt, b, job );

  err = 0.0;
  for ( i = 0; i < n; i++ )
  {
    err = err + fabs ( x[i] - b[i] );
  }

  printf ( "  Revised     Parallel   %8d  %10.4e  %10.4e\n", n, err, wtime );

  free ( a );
  free ( b );
  free ( ipvt );
  free ( x );

  return;
}
/******************************************************************************/

void test03 ( int n )

/******************************************************************************/
/*
  Purpose:

    TEST03 runs the revised version of SGEFA in sequential mode.

  Modified:

    07 April 2008

  Author:

    John Burkardt
*/
{
  float *a;
  float *b;
  float err;
  int i;
  int info;
  int *ipvt;
  int job;
  int lda;
  double wtime;
  float *x;
/*
  Generate the linear system A * x = b.
*/
  lda = n;
  a = ( float * ) malloc ( lda * n * sizeof ( float ) );
  b = ( float * ) malloc ( n * sizeof ( float ) );
  x = ( float * ) malloc ( n * sizeof ( float ) );

  matgen ( lda, n, a, x, b );
/*
  Factor the linear system.
*/
  ipvt = ( int * ) malloc ( n * sizeof ( int ) );

  wtime = omp_get_wtime ( );
  info = msgefa2 ( a, lda, n, ipvt );
  wtime = omp_get_wtime ( ) - wtime;

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST03 - Fatal error!\n" );
    printf ( "  MSGEFA2 reports the matrix is singular.\n" );
    exit ( 1 );
  }
/*
  Solve the linear system.
*/
  job = 0;
  sgesl ( a, lda, n, ipvt, b, job );

  err = 0.0;
  for ( i = 0; i < n; i++ )
  {
    err = err + fabs ( x[i] - b[i] );
  }

  printf ( "  Revised   Sequential   %8d  %10.4e  %10.4e\n", n, err, wtime );

  free ( a );
  free ( b );
  free ( ipvt );
  free ( x );

  return;
}
/******************************************************************************/

int isamax ( int n, float x[], int incx )

/******************************************************************************/
/*
  Purpose:

    ISAMAX finds the index of the vector element of maximum absolute value.

  Discussion:

    WARNING: This index is a 1-based index, not a 0-based index!

  Modified:

    07 April 2008

  Author:

    FORTRAN77 original version by Lawson, Hanson, Kincaid, Krogh.
    C version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979,
    ISBN13: 978-0-898711-72-1,
    LC: QA214.L56.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Algorithm 539: 
    Basic Linear Algebra Subprograms for Fortran Usage,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, September 1979, pages 308-323.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, float X[*], the vector to be examined.

    Input, int INCX, the increment between successive entries of SX.

    Output, int ISAMAX, the index of the element of maximum
    absolute value.
*/
{
  float xmax;
  int i;
  int ix;
  int value;

  value = 0;

  if ( n < 1 || incx <= 0 )
  {
    return value;
  }

  value = 1;

  if ( n == 1 )
  {
    return value;
  }

  if ( incx == 1 )
  {
    xmax = fabs ( x[0] );

    for ( i = 1; i < n; i++ )
    {
      if ( xmax < fabs ( x[i] ) )
      {
        value = i + 1;
        xmax = fabs ( x[i] );
      }
    }
  }
  else
  {
    ix = 0;
    xmax = fabs ( x[0] );
    ix = ix + incx;

    for ( i = 1; i < n; i++ )
    {
      if ( xmax < fabs ( x[ix] ) )
      {
        value = i + 1;
        xmax = fabs ( x[ix] );
      }
      ix = ix + incx;
    }
  }

  return value;
}
/*******************************************************************************/

void matgen ( int lda, int n, float a[], float x[], float b[] )

/*******************************************************************************/
/* 
  Purpose:

    MATGEN generates a "random" matrix for testing.

  Modified:

    27 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, int LDA, the leading dimension of the matrix.

    Input, int N, the order of the matrix, and the length of the vector.

    Output, float A[LDA*N], the matrix.

    Output, float X[N], the solution vector.

    Output, float B[N], the right hand side vector.
*/
{
  int i;
  int j;
  int seed;
  float value;

  seed = 1325;
/*
  Set the matrix A.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      seed = ( 3125 * seed ) % 65536;
      value = ( ( float ) seed - 32768.0 ) / 16384.0;
      a[i+j*lda] = value;
    }
  }
/*
  Set x.
*/
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( float ) ( i + 1 ) / ( ( float ) n );
  }
/*
  Set b = A * x.
*/
  for ( i = 0; i < n; i++ ) 
  {
    b[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i] = b[i] + a[i+j*lda] * x[j];
    }
  }
  return;
}
/******************************************************************************/

void msaxpy ( int nr, int nc, float a[], int n, float x[], float y[] )

/******************************************************************************/
/*
  Purpose:

    MSAXPY carries out multiple "SAXPY" operations.

  Discussion:

    This routine carries out the step of Gaussian elimination where multiples
    of the pivot row are added to the rows below the pivot row.

    A single call to MSAXPY replaces multiple calls to SAXPY.

  Modified:

    07 April 2008

  Author:

    Wesley Petersen

  Parameters:

    Input, int NR, NC, the number of rows and columns in the matrix.

    Input, float A[*], ...

    Input, int N, ...

    Input, float X[*], ...

    Output, float Y[*], ...
*/
{
  int i,j;

# pragma omp parallel \
  shared ( a, nc, nr, x, y ) \
  private ( i, j )

# pragma omp for
  for ( j = 0; j < nc; j++)
  {
    for ( i = 0; i < nr; i++ )
    {
      y[i+j*n] += a[j*n] * x[i];
    }
  }
  return;
}
/******************************************************************************/

void msaxpy2 ( int nr, int nc, float a[], int n, float x[], float y[] )

/******************************************************************************/
/*
  Purpose:

    MSAXPY2 carries out multiple "SAXPY" operations.

  Discussion:

    This routine carries out the step of Gaussian elimination where multiples
    of the pivot row are added to the rows below the pivot row.

    A single call to MSAXPY replaces multiple calls to SAXPY.

  Modified:

    07 April 2008

  Author:

    Wesley Petersen

  Parameters:

    Input, int NR, NC, the number of rows and columns in the matrix.

    Input, float A[*], ...

    Input, int N, ...

    Input, float X[*], ...

    Output, float Y[*], ...
*/
{
  int i,j;

  for ( j = 0; j < nc; j++)
  {
    for ( i = 0; i < nr; i++ )
    {
      y[i+j*n] += a[j*n] * x[i];
    }
  }
  return;
}
/******************************************************************************/

int msgefa ( float a[], int lda, int n, int ipvt[] )

/******************************************************************************/
/* 
  Purpose:

    MSGEFA factors a matrix by gaussian elimination.

  Discussion:

    Matrix references which would, mathematically, be written A(I,J)
    must be written here as:
    * A[I+J*LDA], when the value is needed, or
    * A+I+J*LDA, when the address is needed.

    This variant of SGEFA uses OpenMP for improved parallel execution.
    The step in which multiples of the pivot row are added to individual
    rows has been replaced by a single call which updates the entire
    matrix sub-block.

  Modified:

    07 March 2008

  Author:

    FORTRAN77 original version by Cleve Moler.
    C version by Wesley Petersen.

  Reference:

    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979,
    ISBN13: 978-0-898711-72-1,
    LC: QA214.L56.

  Parameters:

    Input/output, float A[LDA*N].  On input, the matrix to be factored.
    On output, an upper triangular matrix and the multipliers which were 
    used to obtain it.  The factorization can be written A = L * U where
    L is a product of permutation and unit lower triangular matrices and
    U is upper triangular.

    Input, int LDA, the leading dimension of the matrix.

    Input, int N, the order of the matrix.

    Output, int IPVT[N], the pivot indices.

    Output, int MSGEFA, indicates singularity.
    If 0, this is the normal value, and the algorithm succeeded.
    If K, then on the K-th elimination step, a zero pivot was encountered.
    The matrix is numerically not invertible.
*/
{
  int info;
  int k,kp1,l,nm1;
  float t;

  info = 0;
  nm1 = n - 1;
  for ( k = 0; k < nm1; k++ )
  {
    kp1 = k + 1;
    l = isamax ( n-k, a+k+k*lda, 1 ) + k - 1;
    ipvt[k] = l + 1;

    if ( a[l+k*lda] == 0.0 )
    {
      info = k + 1;
      return info;
    }

    if ( l != k )
    {
      t          = a[l+k*lda];
      a[l+k*lda] = a[k+k*lda];
      a[k+k*lda] = t;
    }
    t = -1.0 / a[k+k*lda]; 
    sscal ( n-k-1, t, a+kp1+k*lda, 1 );
/*
  Interchange the pivot row and the K-th row.
*/
    if ( l != k )
    {
      sswap ( n-k-1, a+l+kp1*lda, lda, a+k+kp1*lda, lda );
    }
/*
  Add multiples of the K-th row to rows K+1 through N.
*/
    msaxpy ( n-k-1, n-k-1, a+k+kp1*lda, n, a+kp1+k*lda, a+kp1+kp1*lda );
  }

  ipvt[n-1] = n;

  if ( a[n-1+(n-1)*lda] == 0.0 )
  {
    info = n;
  }

  return info;
}
/******************************************************************************/

int msgefa2 ( float a[], int lda, int n, int ipvt[] )

/******************************************************************************/
/* 
  Purpose:

    MSGEFA2 factors a matrix by gaussian elimination.

  Discussion:

    Matrix references which would, mathematically, be written A(I,J)
    must be written here as:
    * A[I+J*LDA], when the value is needed, or
    * A+I+J*LDA, when the address is needed.

    This variant of SGEFA uses OpenMP for improved parallel execution.
    The step in which multiples of the pivot row are added to individual
    rows has been replaced by a single call which updates the entire
    matrix sub-block.

  Modified:

    07 March 2008

  Author:

    FORTRAN77 original version by Cleve Moler.
    C version by Wesley Petersen.

  Reference:

    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979,
    ISBN13: 978-0-898711-72-1,
    LC: QA214.L56.

  Parameters:

    Input/output, float A[LDA*N].  On input, the matrix to be factored.
    On output, an upper triangular matrix and the multipliers which were 
    used to obtain it.  The factorization can be written A = L * U where
    L is a product of permutation and unit lower triangular matrices and
    U is upper triangular.

    Input, int LDA, the leading dimension of the matrix.

    Input, int N, the order of the matrix.

    Output, int IPVT[N], the pivot indices.

    Output, int MSGEFA, indicates singularity.
    If 0, this is the normal value, and the algorithm succeeded.
    If K, then on the K-th elimination step, a zero pivot was encountered.
    The matrix is numerically not invertible.
*/
{
  int info;
  int k,kp1,l,nm1;
  float t;

  info = 0;
  nm1 = n - 1;
  for ( k = 0; k < nm1; k++ )
  {
    kp1 = k + 1;
    l = isamax ( n-k, a+k+k*lda, 1 ) + k - 1;
    ipvt[k] = l + 1;

    if ( a[l+k*lda] == 0.0 )
    {
      info = k + 1;
      return info;
    }

    if ( l != k )
    {
      t          = a[l+k*lda];
      a[l+k*lda] = a[k+k*lda];
      a[k+k*lda] = t;
    }
    t = -1.0 / a[k+k*lda]; 
    sscal ( n-k-1, t, a+kp1+k*lda, 1 );
/*
  Interchange the pivot row and the K-th row.
*/
    if ( l != k )
    {
      sswap ( n-k-1, a+l+kp1*lda, lda, a+k+kp1*lda, lda );
    }
/*
  Add multiples of the K-th row to rows K+1 through N.
*/
    msaxpy2 ( n-k-1, n-k-1, a+k+kp1*lda, n, a+kp1+k*lda, a+kp1+kp1*lda );
  }

  ipvt[n-1] = n;

  if ( a[n-1+(n-1)*lda] == 0.0 )
  {
    info = n;
  }

  return info;
}
/******************************************************************************/

void saxpy ( int n, float a, float x[], int incx, float y[], int incy )

/******************************************************************************/
/*
  Purpose:

    SAXPY computes float constant times a vector plus a vector.

  Discussion:

    This routine uses unrolled loops for increments equal to one.

  Modified:

    23 February 2006

  Author:

    FORTRAN77 original version by Dongarra, Moler, Bunch, Stewart.
    C version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979,
    ISBN13: 978-0-898711-72-1,
    LC: QA214.L56.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539, 
    ACM Transactions on Mathematical Software, 
    Volume 5, Number 3, September 1979, pages 308-323.

  Parameters:

    Input, int N, the number of elements in X and Y.

    Input, float A, the multiplier of X.

    Input, float X[*], the first vector.

    Input, int INCX, the increment between successive entries of X.

    Input/output, float Y[*], the second vector.
    On output, Y[*] has been replaced by Y[*] + A * X[*].

    Input, int INCY, the increment between successive entries of Y.
*/
{
  int i;
  int ix;
  int iy;
  int m;

  if ( n <= 0 )
  {
    return;
  }

  if ( a == 0.0 )
  {
    return;
  }
/*
  Code for unequal increments or equal increments
  not equal to 1.
*/
  if ( incx != 1 || incy != 1 )
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      y[iy] = y[iy] + a * x[ix];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
/*
  Code for both increments equal to 1.
*/
  else
  {
    m = n % 4;

    for ( i = 0; i < m; i++ )
    {
      y[i] = y[i] + a * x[i];
    }

    for ( i = m; i < n; i = i + 4 )
    {
      y[i  ] = y[i  ] + a * x[i  ];
      y[i+1] = y[i+1] + a * x[i+1];
      y[i+2] = y[i+2] + a * x[i+2];
      y[i+3] = y[i+3] + a * x[i+3];
    }
  }

  return;
}
/******************************************************************************/

float sdot ( int n, float x[], int incx, float y[], int incy )

/******************************************************************************/
/*
  Purpose:

    SDOT forms the dot product of two vectors.

  Discussion:

    This routine uses unrolled loops for increments equal to one.

  Modified:

    23 February 2006

  Author:

    FORTRAN77 original version by Dongarra, Moler, Bunch, Stewart
    C version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539, 
    ACM Transactions on Mathematical Software, 
    Volume 5, Number 3, September 1979, pages 308-323.

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, float X[*], the first vector.

    Input, int INCX, the increment between successive entries in X.

    Input, float Y[*], the second vector.

    Input, int INCY, the increment between successive entries in Y.

    Output, float SDOT, the sum of the product of the corresponding
    entries of X and Y.
*/
{
  int i;
  int ix;
  int iy;
  int m;
  float temp;

  temp = 0.0;

  if ( n <= 0 )
  {
    return temp;
  }
/*
  Code for unequal increments or equal increments
  not equal to 1.
*/
  if ( incx != 1 || incy != 1 )
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      temp = temp + x[ix] * y[iy];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
/*
  Code for both increments equal to 1.
*/
  else
  {
    m = n % 5;

    for ( i = 0; i < m; i++ )
    {
      temp = temp + x[i] * y[i];
    }

    for ( i = m; i < n; i = i + 5 )
    {
      temp = temp + x[i  ] * y[i  ] 
                  + x[i+1] * y[i+1] 
                  + x[i+2] * y[i+2] 
                  + x[i+3] * y[i+3] 
                  + x[i+4] * y[i+4];
    }
  }

  return temp;
}
/*******************************************************************************/

int sgefa ( float a[], int lda, int n, int ipvt[] )

/*******************************************************************************/
/*
  Purpose:

    SGEFA factors a matrix by gaussian elimination.

  Discussion:

    Matrix references which would, mathematically, be written A(I,J)
    must be written here as:
    * A[I+J*LDA], when the value is needed, or
    * A+I+J*LDA, when the address is needed.

  Modified:

    07 March 2008

  Author:

    FORTRAN77 original version by Cleve Moler.
    C version by John Burkardt.

  Reference:

    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979,
    ISBN13: 978-0-898711-72-1,
    LC: QA214.L56.

  Parameters:

    Input/output, float A[LDA*N].  On input, the matrix to be factored.
    On output, an upper triangular matrix and the multipliers which were 
    used to obtain it.  The factorization can be written A = L * U where
    L is a product of permutation and unit lower triangular matrices and
    U is upper triangular.

    Input, int LDA, the leading dimension of the matrix.

    Input, int N, the order of the matrix.

    Output, int IPVT[N], the pivot indices.

    Output, int SGEFA, indicates singularity.
    If 0, this is the normal value, and the algorithm succeeded.
    If K, then on the K-th elimination step, a zero pivot was encountered.
    The matrix is numerically not invertible.
*/
{
  int j;
  int info;
  int k;
  int l;
  float t;

  info = 0;

  for ( k = 1; k <= n - 1; k++ ) 
  {
/* 
  Find l = pivot index.
*/
    l = isamax ( n-k+1, &a[k-1+(k-1)*lda], 1 ) + k - 1;
    ipvt[k-1] = l;
/* 
  Zero pivot implies this column already triangularized.
*/
    if ( a[l-1+(k-1)*lda] != 0.0 ) 
    {
/* 
  Interchange if necessary.
*/
      if ( l != k ) 
      {
        t                = a[l-1+(k-1)*lda];
        a[l-1+(k-1)*lda] = a[k-1+(k-1)*lda];
        a[k-1+(k-1)*lda] = t; 
      }
/* 
  Compute multipliers.
*/
      t = - 1.0 / a[k-1+(k-1)*lda];
      sscal ( n-k, t, &a[k+(k-1)*lda], 1 );
/* 
  Row elimination with column indexing.
*/
      for ( j = k + 1; j <= n; j++ ) 
      {
        t = a[l-1+(j-1)*lda];
        if (l != k) 
        {
          a[l-1+(j-1)*lda] = a[k-1+(j-1)*lda];
          a[k-1+(j-1)*lda] = t;
        }
        saxpy ( n-k, t, &a[k+(k-1)*lda], 1, &a[k+(j-1)*lda], 1 );
      } 
    }
    else
    { 
      info = k;
    }
  } 
  ipvt[n-1] = n;

  if (a[n-1+(n-1)*lda] == 0.0 ) 
  {
    info = n - 1;
  }
  return info;
}
/******************************************************************************/

void sgesl ( float a[], int lda, int n, int ipvt[], float b[], int job )

/******************************************************************************/
/*
  Purpose:

    SGESL solves a real general linear system A * X = B.

  Discussion:

    SGESL can solve either of the systems A * X = B or A' * X = B.

    The system matrix must have been factored by SGECO or SGEFA.

    A division by zero will occur if the input factor contains a
    zero on the diagonal.  Technically this indicates singularity
    but it is often caused by improper arguments or improper
    setting of LDA.  It will not occur if the subroutines are
    called correctly and if SGECO has set 0.0 < RCOND
    or SGEFA has set INFO == 0.

  Modified:

    04 April 2006

  Author:

    FORTRAN77 original by Dongarra, Moler, Bunch and Stewart.
    C translation by John Burkardt.

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, (Society for Industrial and Applied Mathematics),
    3600 University City Science Center,
    Philadelphia, PA, 19104-2688.
    ISBN: 0-89871-172-X

  Parameters:

    Input, float A[LDA*N], the output from SGECO or SGEFA.

    Input, int LDA, the leading dimension of A.

    Input, int N, the order of the matrix A.

    Input, int IPVT[N], the pivot vector from SGECO or SGEFA.

    Input/output, float B[N].
    On input, the right hand side vector.
    On output, the solution vector.

    Input, int JOB.
    0, solve A * X = B;
    nonzero, solve A' * X = B.
*/
{
  int k;
  int l;
  float t;
/*
  Solve A * X = B.
*/
  if ( job == 0 )
  {
    for ( k = 1; k <= n-1; k++ )
    {
      l = ipvt[k-1];
      t = b[l-1];

      if ( l != k )
      {
        b[l-1] = b[k-1];
        b[k-1] = t;
      }
      saxpy ( n-k, t, a+k+(k-1)*lda, 1, b+k, 1 );
    }

    for ( k = n; 1 <= k; k-- )
    {
      b[k-1] = b[k-1] / a[k-1+(k-1)*lda];
      t = -b[k-1];
      saxpy ( k-1, t, a+0+(k-1)*lda, 1, b, 1 );
    }
  }
/*
  Solve A' * X = B.
*/
  else
  {
    for ( k = 1; k <= n; k++ )
    {
      t = sdot ( k-1, a+0+(k-1)*lda, 1, b, 1 );
      b[k-1] = ( b[k-1] - t ) / a[k-1+(k-1)*lda];
    }

    for ( k = n-1; 1 <= k; k-- )
    {
      b[k-1] = b[k-1] + sdot ( n-k, a+k+(k-1)*lda, 1, b+k, 1 );
      l = ipvt[k-1];

      if ( l != k )
      {
        t = b[l-1];
        b[l-1] = b[k-1];
        b[k-1] = t;
      }
    }
  }
  return;
}
/******************************************************************************/

void sscal ( int n, float sa, float x[], int incx )

/******************************************************************************/
/*
  Purpose:

    SSCAL scales a float vector by a constant.

  Modified:

    23 February 2006

  Author:

    Jack Dongarra
    C version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979,
    ISBN13: 978-0-898711-72-1,
    LC: QA214.L56.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, September 1979, pages 308-323.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, float SA, the multiplier.

    Input/output, float X[*], the vector to be scaled.

    Input, int INCX, the increment between successive entries of X.
*/
{
  int i;
  int ix;
  int m;

  if ( n <= 0 )
  {
  }
  else if ( incx == 1 )
  {
    m = n % 5;

    for ( i = 0; i < m; i++ )
    {
      x[i] = sa * x[i];
    }

    for ( i = m; i < n; i = i + 5 )
    {
      x[i]   = sa * x[i];
      x[i+1] = sa * x[i+1];
      x[i+2] = sa * x[i+2];
      x[i+3] = sa * x[i+3];
      x[i+4] = sa * x[i+4];
    }
  }
  else
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    for ( i = 0; i < n; i++ )
    {
      x[ix] = sa * x[ix];
      ix = ix + incx;
    }

  }

  return;
}
/******************************************************************************/

void sswap ( int n, float x[], int incx, float y[], int incy )

/******************************************************************************/
/*
  Purpose:

    SSWAP interchanges two float vectors.

  Modified:

    23 February 2006

  Author:

    C version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979,
    ISBN13: 978-0-898711-72-1,
    LC: QA214.L56.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539, 
    ACM Transactions on Mathematical Software, 
    Volume 5, Number 3, September 1979, pages 308-323.

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input/output, float X[*], one of the vectors to swap.

    Input, int INCX, the increment between successive entries of X.

    Input/output, float Y[*], one of the vectors to swap.

    Input, int INCY, the increment between successive elements of Y.
*/
{
  int i;
  int ix;
  int iy;
  int m;
  float temp;

  if ( n <= 0 )
  {
  }
  else if ( incx == 1 && incy == 1 )
  {
    m = n % 3;

    for ( i = 0; i < m; i++ )
    {
      temp = x[i];
      x[i] = y[i];
      y[i] = temp;
    }

    for ( i = m; i < n; i = i + 3 )
    {
      temp = x[i];
      x[i] = y[i];
      y[i] = temp;

      temp = x[i+1];
      x[i+1] = y[i+1];
      y[i+1] = temp;

      temp = x[i+2];
      x[i+2] = y[i+2];
      y[i+2] = temp;
    }
  }
  else
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      temp = x[ix];
      x[ix] = y[iy];
      y[iy] = temp;
      ix = ix + incx;
      iy = iy + incy;
    }
  }
  return;
}
/******************************************************************************/

void timestamp ( )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
