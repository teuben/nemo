
/*
 *   pow.c:  some specialized & optimized versions of the pow() function
 *
 *   powi(double,int)      especially for low order
 *   powd(double,double)   to handle NaN's
 *     
 */

#include <stdinc.h>

/* 
 * powi(x,p):   like pow(x,p) except p is a lowish order integer
 *    this has the advantage that:
 *     1) x < 0  and it does the right thing. For odd p's  sign(powi) = sign(x)
 *     2) p < 1  
 *  Note: returns 0.0 if x=0 and p<0
 *        returns 1.0 if p=0
 *     
 */

double powi(double x, int p)
{
  int i;
  real prod;

  if (p==0) return 1;
  if (x==0.0) return 0.0;
  if (p<0) {
    p = -p;
    x = 1.0/x;
  } 
  if (p==1) return x;
  prod = x*x;
  for (i=2; i<p; i++)
    prod *= x;
  return prod;
}


/* 
 * powd(x,p):  we need  our own, to prevent NaN's for x < 0
 */

double powd(double x, double p)
{
  if (x<=0) return 0.0;
  return pow(x,p);
}


#ifdef TESTBED

#include <getparam.h>

string defv[] = {
  "x=2.0\n       x from powi(x,p)",
  "p=1\n         p from powi(x,p)",
  "powi=t\n      Use pow() or our powi() ? ",
  "n1=1\n        Inner iteration count for benchmark",
  "n2=1\n        Outer Iteration count",
  "VERSION=1\n   25-jan-05 pjt",
  NULL,
};

string usage="test pow";

void nemo_main(void)
{
  double sum, xp, x = getdparam("x");
  int p = getiparam("p");
  int i1,n1=getiparam("n1");
  int i2,n2=getiparam("n2");
  bool Qpowi=getbparam("powi");

  xp = (double) p;
  if (n1==1 && n2==1) 
    sum = powi(x,p);
  else {
    sum = 0.0;
    if (Qpowi) {
      for (i2=0; i2<n2; i2++)
	for (i1=0; i1<n1; i1++)
	  sum += powi(x,p);
    } else {
      for (i2=0; i2<n2; i2++)
	for (i1=0; i1<n1; i1++)
	  sum += pow(x,xp);
    }
  }
  printf("%lg\n",sum);
    
}

#endif
