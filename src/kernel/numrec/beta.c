/* NEMO version of some NR gamma functions */
/* float's are now "real", so depending on the compilation float or double */

#include <stdinc.h>

extern real gammln(real xx);

#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30



real betacf(real a, real b, real x)
{
	int m,m2;
	real aa,c,d,del,h,qab,qam,qap;

	qab=a+b;
	qap=a+1.0;
	qam=a-1.0;
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=MAXIT;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (m > MAXIT) error("betacf: a or b too big, or MAXIT too small");
	return h;
}


real betai(real a, real b, real x)
{
	real bt;

	if (x < 0.0 || x > 1.0) error("Bad x=%g in routine betai",x);
	if (x == 0.0 || x == 1.0) bt=0.0;
	else
		bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0))
		return bt*betacf(a,b,x)/a;
	else
		return 1.0-bt*betacf(b,a,1.0-x)/b;
}

#undef MAXIT
#undef EPS
#undef FPMIN

/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */

