/*
 * special function and some mathematical tools for P2M2
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vtclocal.h"

#include "design.c" /* position sets of pseudo particles obtained from 
                       Hardin and Sloane's
                       http://www.research.att.com/~njas/sphdesigns/ */

void
load_design(int order, double ss, int full_dof, int negativemass,
	    int *npp, double (**pppos)[3])
{
    int i, k, n;

    typedef struct
    {
	int n;       /* # of pseudo particles */
	int d;       /* sperical-design order */
	double *pos; /* position coordinates */
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
	des[1].n = 1;
	des[1].d = 0;
	des[1].pos = NULL;
	des[2].n = 3;
	des[2].d = 0;
	des[2].pos = NULL;
    }
    if (order >= sizeof(des)/sizeof(Design)) {
	fprintf(stderr, "too large expansion order p: %d\n", order);
	fprintf(stderr, "max p: %d\n", sizeof(des)/sizeof(Design)-1);
	exit(1);
    }
    *npp = n = des[order].n;
    if (negativemass) {
	(*npp) *= 2;
    }
    *pppos = (double (*)[3])malloc(sizeof(double)*3*(*npp));
    if (full_dof == 0 || order > 2) {
	for (i = 0; i < *npp; i++) {
	    for (k = 0; k < 3; k++) {
		(*pppos)[i][k] = des[order].pos[(i%n)*3+k]*ss;
	    }
	}
    }
}

/*
 * returns Legendre polynominals of the n-th order
 * P_0(x)+P_1(x)+..+P_{n-1}(x)
 */
void plgndr0(int n, double x, double *pln)
{
    int i;
    double p0, p1, p2;

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
	p2 = (x*(2*i-1)*p1-(i-1)*p0)/(double)i;
	p0 = p1;
	p1 = p2;
	pln[i] = p2;
    }
}

/*
 * returns Pl,m(x)
 * 0 <= m <= l
 * -1 <= x <= 1
 */
double plgndr(int l, int m, double x)
{
    void nrerror(char error_text[]);
    double fact, pll, pmm, pmmp1, somx2;
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

void ylm(int l, int m, double theta, double phi,
	 double *re, double *im)
{
    int i;
    double coeff = 1.0;
    double plm, ylm;
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

#define ITERATIONMAX (50)
#define ROTATE(a,i,j,k,l) {g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau); a[k][l]=h+s*(g-h*tau);}

#define TEST0 (0)

/*
 * Computes all eigenvalues and eigenvectors of a double symmetric matrix
 * a[1..n][1..n]. On output, elements of a above the diagonal are
 * destroyed. d[1..n] returns the eigenvalues of a. v[1..n][1..n] is a
 * matrix whose columns contain, on output, the normalized eigenvectors
 * of a. nrot returns the number of Jacobi rotations that were required.
 * return 0 on successful execution, otherwise -1.
 */

int jacobi(double (*a)[3], double d[3], double (*v)[3], int *nrot)
{
    int j,iq,ip,i;
    double tresh,theta,tau,t,sm,s,h,g,c;
    static double b[3];
    static double z[3];

    for (ip = 0; ip < 3; ip++) {
	for (iq = 0; iq < 3; iq++) {
	    v[ip][iq]=0.0;
	}
	v[ip][ip]=1.0;
    }

    for (ip = 0; ip < 3; ip++) {
	b[ip] = d[ip] = a[ip][ip];
	z[ip] = 0.0;
    }

    *nrot=0;
    for (i = 0; i < ITERATIONMAX; i++) {

	sm = 0.0;
	for (ip = 0; ip < 3-1; ip++) {
	    for (iq = ip+1; iq < 3; iq++) {
		sm += fabs(a[ip][iq]);
	    }
	}

        /* the normal return, which relies on quadratic convergence to
	   machine underflow. */
	if (sm == 0.0) {
	    return (0);
	}
	if (i < 3) {
	    tresh=0.2*sm/(3*3);
	} else {
	    tresh=0.0;
	}
	for (ip = 0; ip < 3-1; ip++) {
	    for (iq = ip+1; iq < 3; iq++) {
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
		    for (j = iq+1; j < 3; j++) {
			ROTATE(a,ip,j,iq,j);
		    }
		    for (j = 0; j < 3; j++) {
			ROTATE(v,j,ip,j,iq);
		    }
		    ++(*nrot);
		}
	    }
	}
	for (ip = 0; ip < 3; ip++) {
	    b[ip] += z[ip];
	    d[ip] = b[ip];
	    z[ip] = 0.0;
	}
    }
    fprintf(stderr, "jacobi(): too many iterations\n");
    return (-1);
}

void eigenvalsort(double d[3], double (*v)[3])
{
    int i, j, k;
    double p;

    for (i = 0; i < 3-1; i++) {
	k = i;
	p = d[k];
	for (j = i+1; j < 3; j++) {
	    if (d[j] > p) {
		k = j;
		p = d[k];
	    }
	}
	if (k != i) {
	    d[k] = d[i];
	    d[i] = p;
	    for (j = 0; j < 3; j++) {
		p = v[j][i];
		v[j][i] = v[j][k];
		v[j][k] = p;
	    }
	}
    }
}

void eigenvalsort_reverse(double d[3], double (*v)[3])
{
    int i, j, k;
    double p;

    for (i = 0; i < 3-1; i++) {
	k = i;
	p = d[k];
	for (j = i+1; j < 3; j++) {
	    if (d[j] < p) {
		k = j;
		p = d[k];
	    }
	}
	if (k != i) {
	    d[k] = d[i];
	    d[i] = p;
	    for (j = 0; j < 3; j++) {
		p = v[j][i];
		v[j][i] = v[j][k];
		v[j][k] = p;
	    }
	}
    }
}


#if TEST0

void
printmatrix(double a[3][3])
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
    static double src[3][3] = {
	10.0, 1.0, 1.0,
	1.0, 1.0, 0.0,
	1.0, 0.0, 1.0,
    };
    static double dst[3][3], tmp[3][3];
    static double evec[3][3];
    static double eval[3];
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
