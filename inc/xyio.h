/*
 * XYIO stuff - modeled after miriad's xyio
 * assumed image.h has been loaded
 */
image *xyopen  ( int *, int, int *, string,int *);
void   xyclose ( int );
void   xyread  ( int, int,  real *);
void   xywrite ( int, int, real *);
void   xysetpl ( int, int, int *);
