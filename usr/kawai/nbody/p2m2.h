#define DIM (3)
#define real double

#if defined(__cplusplus)
extern "C"
{
#endif

void plgndr0(int n, real x, real *pln);
real plgndr(int l, int m, real x);
void ylm(int l, int m, real theta, real phi, real *re, real *im);
int jacobi(real (*a)[DIM], real d[DIM], real (*v)[DIM], int *nrot);
void eigenvalsort(real d[DIM], real (*v)[DIM]);

#if defined(__cplusplus)
}
#endif
