/*
 * special function and some other tools
 * for P2M2
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include  <stdiostream.h>

#define real double

#include "vector.h"
#include "nbody_particle.h"
#include "p2m2.h"
#ifdef P2M2

#include "design.c"
void nbody_system::load_design(int order, real ss, int full_dof)
{
    typedef struct design_t
    {
	int n;
	int d;
	real *pppos;
    } Design;

    static Design des[] =
    {
	{1, 1, pppos1},
	{4, 2, pppos4},
	{12, 5, pppos12},
	{24, 7, pppos24},
	{36, 8, pppos36},
	{60, 10, pppos60},
	{84, 12, pppos84},
    };

    if (full_dof == 1) {
	p2m2.full_dof = 1;
	des[1].n = 1;
	des[1].d = 0;
	des[1].pppos = NULL;
	des[2].n = 3;
	des[2].d = 0;
	des[2].pppos = NULL;
    }
    else {
	p2m2.full_dof = 0;
    }

    if (order >= sizeof(des)/sizeof(Design)) {
	cerr << "too large expansion order p: " << order;
	cerr << "  max p = " << sizeof(des)/sizeof(Design)-1 << endl;
	exit(1);
    }

    p2m2.order = order;
    p2m2.npp = des[p2m2.order].n;

    if (full_dof == 0 || order > 2) {
	p2m2.spherical_design = des[p2m2.order].d;
	p2m2.ppscale = ss;
	for (int i = 0; i < p2m2.npp; i++) {
	    for (int k = 0; k < 3; k++) {
		p2m2.pppos0[i][k] = des[p2m2.order].pppos[i*3+k]*p2m2.ppscale;
	    }
	}
    }
}

#endif // P2M2

/*
 * returns Legendre polynominals of the n-th order
 * P_0(x)+P_1(x)+..+P_{n-1}(x)
 */
void plgndr0(int n, real x, real *pln)
{
    int i;
    real p0, p1, p2;

/*
    printf("#####  %d    %f\n", n, x);
    */

    pln[0] = p0 = 1.0;
    pln[1] = p1 = x;

    if (n < 3)
    {
	return;
    }
    for (i = 2; i < n; i++)
    {
	p2 = (x*(2*i-1)*p1-(i-1)*p0)/(real)i;
	p0 = p1;
	p1 = p2;
	pln[i] = p2;
    }

#if 0
    fprintf(stderr, "P(%d, %5.3f): ", n, x);
    for (i = 0; i < n; i++)
    {
	fprintf(stderr, "%5.3f ", pln[i]);
    }
    fprintf(stderr, "\n");
#endif
}

/*
 * returns Pl,m(x)
 * 0 <= m <= l
 * -1 <= x <= 1
 */
real plgndr(int l, int m, real x)
{
    void nrerror(char error_text[]);
    real fact, pll, pmm, pmmp1, somx2;
    int i,ll;

    if (m < 0 || m > l || fabs(x) > 1.0)
    {
	fprintf(stderr, "plgndr(): bad argument\n");
	exit(1);
    }

    /* Pm,m(x) */
    pmm=1.0;
    if (m > 0)
    {
	somx2=sqrt((1.0-x)*(1.0+x));
	fact=1.0;
	for (i=1; i<=m; i++)
	{
	    pmm *= -fact*somx2;
	    fact += 2.0;
	}
    }

    if (l == m)
    {
	return pmm;
    }
    else
    {
	/* Pm+1,m(x) */
	pmmp1=x*(2*m+1)*pmm;
	if (l == (m+1))
	{
	    return pmmp1;
	}
	else
	{
	    /* Pl,m(x) */
	    for (ll=m+2; ll<=l; ll++)
	    {
		pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
		pmm=pmmp1;
		pmmp1=pll;
	    }
	    return pll;
	}
    }
}

void ylm(int l, int m, real theta, real phi,
	 real *re, real *im)
{
    int i;
    real coeff = 1.0;
    real plm, ylm;
    int sign = 1;

    if (m < 0.0)
    {
	sign *= -1;
    }
    m = abs(m);

    for (i = l-m; i > 1; i--)
    {
	coeff *= i;
    }
    for (i = l+m; i > 1; i--)
    {
	coeff /= i;
    }
    coeff *= (2.0*l+1.0)/4.0/M_PI;
    coeff = sqrt(coeff);
/*
    printf("coeff: %e\n", coeff);
    */

    plm = plgndr(l, m , cos(theta));
    ylm = coeff * plm;
    *re = ylm * cos(m*phi);
    *im = ylm * sin(m*phi);

    if (sign < 0)
    {
	*im *= -1;
	if (m % 2 == 1)
	{
	    *re *= -1;
	    *im *= -1;
	}
    }
}





#define DIM (3)
#define EPS (1e-16)
#define ITERATIONMAX (50)
#define ROTATE(a,i,j,k,l) {g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau); a[k][l]=h+s*(g-h*tau);}

#define TEST0 (0)
/*
 */

/*
 * Computes all eigenvalues and eigenvectors of a real symmetric matrix
 * a[1..n][1..n]. On output, elements of a above the diagonal are
 * destroyed. d[1..n] returns the eigenvalues of a. v[1..n][1..n] is a
 * matrix whose columns contain, on output, the normalized eigenvectors
 * of a. nrot returns the number of Jacobi rotations that were required.
 * return 0 on successful execution, otherwise -1.
 */

int jacobi(real (*a)[DIM], real d[DIM], real (*v)[DIM], int *nrot)
{
    int j,iq,ip,i;
    real tresh,theta,tau,t,sm,s,h,g,c;
    static real b[DIM];
    static real z[DIM];

    for (ip = 0; ip < DIM; ip++) {
	for (iq = 0; iq < DIM; iq++) {
	    v[ip][iq]=0.0;
	}
	v[ip][ip]=1.0;
    }

    for (ip = 0; ip < DIM; ip++) {
	b[ip] = d[ip] = a[ip][ip];
	z[ip] = 0.0;
    }

    *nrot=0;
    for (i = 0; i < ITERATIONMAX; i++) {

	sm = 0.0;
	for (ip = 0; ip < DIM-1; ip++) {
	    for (iq = ip+1; iq < DIM; iq++) {
		sm += fabs(a[ip][iq]);
	    }
	}

        /* the normal return, which relies on quadratic convergence to
	   machine underflow. */
	if (sm == 0.0) {
	    return (0);
	}
	if (i < 3) {
	    tresh=0.2*sm/(DIM*DIM);
	} else {
	    tresh=0.0;
	}
	for (ip = 0; ip < DIM-1; ip++) {
	    for (iq = ip+1; iq < DIM; iq++) {
		g=100.0*fabs(a[ip][iq]);

		if (i > 3 &&
		    (fabs(d[ip])+g) == fabs(d[ip]) &&
		    (fabs(d[iq])+g) == fabs(d[iq])) { 
		    a[ip][iq]=0.0;
		} else if (fabs(a[ip][iq]) > tresh) {
		    h=d[iq]-d[ip];
		    if ((fabs(h)+g) == fabs(h)) {
			t=(a[ip][iq])/h;
		    } else {
			theta=0.5*h/(a[ip][iq]);
			t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
			if (theta < 0.0) {
			    t = -t;
			}
		    }
		    c=1.0/sqrt(1+t*t);
		    s=t*c;
		    tau=s/(1.0+c);
		    h=t*a[ip][iq];
		    z[ip] -= h;
		    z[iq] += h;
		    d[ip] -= h;
		    d[iq] += h;
		    a[ip][iq]=0.0;
		    for (j = 0; j <= ip-1; j++) {
			ROTATE(a,j,ip,j,iq);
		    }
		    for (j = ip+1; j <= iq-1; j++) {
			ROTATE(a,ip,j,j,iq);
		    }
		    for (j = iq+1; j < DIM; j++) {
			ROTATE(a,ip,j,iq,j);
		    }
		    for (j = 0; j < DIM; j++) {
			ROTATE(v,j,ip,j,iq);
		    }
		    ++(*nrot);
		}
	    }
	}
	for (ip = 0; ip < DIM; ip++) {
	    b[ip] += z[ip];
	    d[ip] = b[ip];
	    z[ip] = 0.0;
	}
    }
    fprintf(stderr, "jacobi(): too many iterations\n");
    return (-2);
}

void eigenvalsort(real d[DIM], real (*v)[DIM])
{
    int i, j, k;
    real p;

    for (i = 0; i < DIM-1; i++) {
	k = i;
	p = d[k];
	for (j = i+1; j < DIM; j++) {
	    if (d[j] > p) {
		k = j;
		p = d[k];
	    }
	}
	if (k != i) {
	    d[k] = d[i];
	    d[i] = p;
	    for (j = 0; j < DIM; j++) {
		p = v[j][i];
		v[j][i] = v[j][k];
		v[j][k] = p;
	    }
	}
    }
}


#if TEST0

void
printmatrix(real a[3][3])
{
    int i, j;

    for (i = 0; i < 3; i++)
    {
	for (j = 0; j < 3; j++)
	{
	    printf("%10.6f ", a[i][j]);
	}
	printf("\n");
    }
}

void
main(void)
{
    int n = 3;
    int i, j, k;
    static real src[3][3] = {
	10.0, 1.0, 1.0,
	1.0, 1.0, 0.0,
	1.0, 0.0, 1.0,
    };
    static real dst[3][3], tmp[3][3];
    static real evec[3][3];
    static real eval[3];
    int nrot;

    for (i = 0; i < n; i++) {
	for (j = 0; j < n; j++) {
	    tmp[i][j] = src[i][j];
	}
    }

    printf("tmp0:\n");
    printmatrix(tmp);
    jacobi(tmp, eval, evec, &nrot);
    printf("tmp1:\n");
    printmatrix(tmp);

    for (i = 0; i < n; i++) {
	for (j = 0; j < n; j++) {
	    tmp[i][j] = 0.0;
	    for (k = 0; k < n; k++) {
		tmp[i][j] += src[i][k]*evec[k][j];
	    }
	}
    }
    printf("tmp2:\n");
    printmatrix(tmp);

#if 0
    for (i = 0; i < n; i++) {
	for (j = 0; j < n; j++) {
	    tmp[i][j] = evec[i][j]*eval[j];
	}
    }
    printf("tmp3:\n");
    printmatrix(tmp);
#endif

    for (i = 0; i < n; i++) {
	for (j = 0; j < n; j++) {
	    dst[i][j] = 0.0;
	    for (k = 0; k < n; k++) {
		dst[i][j] += evec[k][i]*tmp[k][j];
	    }
	}
    }

    printf("nrot: %d\n eval: %10.6f %10.6f %10.6f\n",
	   nrot, eval[0], eval[1], eval[2]);


    printf("src:\n");
    printmatrix(src);
    printf("\nevec:\n");
    printmatrix(evec);
    printf("\ndst:\n");
    printmatrix(dst);
    printf("\n");

}
#endif /* TEST0 */
