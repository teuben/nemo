/* NEMO version of some NR gamma functions */
/* float's are now "real", so depending on the compilation float or double */

#include <stdinc.h>


#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30


real gammln(real xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

void gcf(real *gammcf, real a, real x, real *gln)
{
	int i;
	real an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) error("i=%d too large, ITMAX=%d too small in gcf,i,ITMAX");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}

void gser(real *gamser, real a, real x, real *gln)
{
	int n;
	real sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
	        if (x < 0.0) error("x=%g less than 0 in routine gser",x);
	        *gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		error("a too large, ITMAX=%d too small in routine gser",ITMAX);
		return;
	}
}

real gammq(real a, real x)
{
	real gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) error("a=%g, x=%g: Invalid arguments in routine gammq",a,x);
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}



#undef ITMAX
#undef EPS
#undef FPMIN


