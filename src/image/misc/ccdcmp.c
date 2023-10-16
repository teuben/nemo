/*
 *  Student T and F test
 *	8-oct-97	0.2
 *	28-jan-98	0.2a   isolated beta functions in numrec
 *      12-jan-99     0.3 - cloned off tabcmp for images
 *      29-sep-01     0.3a - don't fatal if two images identical
 *       5-dec-2020   0.4 - plotting
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#include <image.h>


string defv[] = {
    "in1=???\n		First image dataset",
    "in2=???\n		Second image dataset",
    "xmin=\n            In case minimum used (need both minmax)",
    "xmax=\n            In case maximum used (need both minmax)",
    "ymin=\n            Min",
    "ymax=\n            Max",
    "mode=0\n           Output mode  0: table  1=histo 2=d1 vs. d2 3=....",
    "bins=\n            Select histogram of in1-in2 with these #bins",
    "tab=\n             Output table",
    "test=f\n           Report on F/T/test",
    "VERSION=0.4\n	5-dec-2020 PJT",
    NULL,
};

string usage="Student T and F test for mean and variance comparisons";

#ifndef MAXHIST
#define MAXHIST 1024
#endif


/* local NR routines */

void ftest(real data1[], unsigned long n1, real data2[], unsigned long n2,
	real *f, real *prob);
void tptest(real data1[], real data2[], unsigned long n, real *t,
	real *prob);
void ttest(real data1[], unsigned long n1, real data2[], unsigned long n2,
	real *t, real *prob);
void tutest(real data1[], unsigned long n1, real data2[], unsigned long n2,
	real *t, real *prob);
void avevar(real data[], unsigned long n, real *ave, real *var);

/* NEMO's NR routines in src/kernel/numrec */

extern real betai(real a, real b, real x);

void nemo_main(void)
{
    int i, j, k, l, npt, npt1, npt2, nsteps, mode = getiparam("mode");
    real xrange[2];
    real *x1 = NULL, *x2 = NULL;
    real t, f, tu, tp, tprob, fprob, tuprob, tpprob;
    real xmin, xmax;
    string input1 = getparam("in1");
    string input2 = getparam("in2");
    stream instr1, instr2;
    imageptr iptr1=NULL, iptr2=NULL;
    real *d1, *d2;
    bool Qauto, Qminmax, Qtest = getbparam("test");

    Qminmax = hasvalue("xmin") && hasvalue("xmax");
    if (Qminmax) {
      xmin = getdparam("xmin");
      xmax = getdparam("xmax");
    }

    instr1 = stropen(input1,"r");
    instr2 = stropen(input2,"r");

    read_image(instr2, &iptr2);
    read_image(instr1, &iptr1);

    x1 = Frame(iptr1);      
    x2 = Frame(iptr2);
    npt1 = Nx(iptr1) * Ny(iptr1) * Nz(iptr1);
    npt2 = Nx(iptr2) * Ny(iptr2) * Nz(iptr2);
    if (npt1 != npt2)
        error("fornow--Images differ in size: %d != %d", npt1, npt2);

    if (Qminmax) {
      d1 = (real *) allocate(npt1*sizeof(real));
      d2 = (real *) allocate(npt2*sizeof(real));
      for (i=0, l=0; i<npt1; i++) {
	if (x1[i]==0.0 || x2[i]==0.0) continue;
	if (xmin<x1[i]<xmax || xmin<x2[i]<xmax) {
	  d1[l] = x1[i];
	  d2[l] = x2[i];
	  l++;
	}
      }
      npt = l;
      dprintf(0,"Found %d/%d between %g %g\n",npt,npt1,xmin,xmax);      
    } else {
      /* find min max from both cubes */
      xmin = xmax = CubeValue(iptr1,0,0,0);
      for (k=0; k<Nz(iptr1); k++)
        for (j=0; j<Ny(iptr1); j++)
	  for (i=0; i<Nx(iptr1); i++) {
	    xmin = MIN(xmin, CubeValue(iptr1,i,j,k));
	    xmax = MAX(xmax, CubeValue(iptr1,i,j,k));
	  }
      for (k=0; k<Nz(iptr2); k++)
        for (j=0; j<Ny(iptr2); j++)
	  for (i=0; i<Nx(iptr2); i++) {
	    xmin = MIN(xmin, CubeValue(iptr2,i,j,k));
	    xmax = MAX(xmax, CubeValue(iptr2,i,j,k));
	  }
      dprintf(0,"[Datamin/max found: %g %g]\n",xmin,xmax);
      d1 = x1;
      d2 = x2;
      npt = npt1;
    }

    /* raw comparisons on the two maps */
    if (Qminmax && Qtest) {
      ftest(d1, npt, d2, npt, &f, &fprob);
      printf("F-test: f=%g  prob=%g\n",f,fprob);

      ttest(d1, npt, d2, npt , &t, &tprob);
      printf("T-test: t=%g  prob=%g\n",t,tprob);

      tutest(d1, npt, d2, npt, &tu, &tuprob);
      printf("TU-test: tu=%g  prob=%g\n",tu,tuprob);

      tptest(d1, d2, npt, &tp, &tpprob);
      printf("TP-test: tp=%g  prob=%g\n",tp,tpprob);
    } else if (Qtest) {
      ftest(x1, npt1, x2, npt2, &f, &fprob);
      printf("F-test: f=%g  prob=%g\n",f,fprob);

      ttest(x1, npt1, x2, npt2, &t, &tprob);
      printf("T-test: t=%g  prob=%g\n",t,tprob);

      tutest(x1, npt1, x2, npt2, &tu, &tuprob);
      printf("TU-test: tu=%g  prob=%g\n",tu,tuprob);

      if (npt1 == npt2) {
        tptest(x1, x2, npt1, &tp, &tpprob);
        printf("TP-test: tp=%g  prob=%g\n",tp,tpprob);
      }
    }

    if (hasvalue("bins")) {
      nsteps=getiparam("bins");
      if (nsteps > MAXHIST)
        error("bins=%d too large; MAXHIST=%d",nsteps,MAXHIST);
      dprintf(0,"Histogram with bins=%d\n",nsteps);
    } else {
      xrange[0] = xmin;
      xrange[1] = xmax;
    }

    if (hasvalue("tab")) {
      dprintf(0,"Output table mode=%d\n",mode);
      stream tstr = stropen(getparam("tab"),"w");

      if (mode == 0) {                       //   0: d1-d2 (histogram)
	for (i=0; i<npt; i++)
	  fprintf(tstr,"%g\n",d1[i]-d2[i]);
      } else if (mode == 1) {                //   1: d1, d2
	for (i=0; i<npt; i++)
	  fprintf(tstr,"%g %g\n",d1[i],d2[i]);
      } else if (mode == 2) {                //   2: d1, d1-d2
	for (i=0; i<npt; i++)
	  fprintf(tstr,"%g %g\n",d1[i],d1[i]-d2[i]);
      } else if (mode == 3) {                //   3: d1, 2*(d1-d2)/(d1+d2)
	for (i=0; i<npt; i++)
	  fprintf(tstr,"%g %g\n",d1[i],2*(d1[i]-d2[i])/(d1[i]+d2[i]));
      }
      strclose(tstr);
    }
      

}


/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */
/* adapted for NEMO					  */


void ftest(real data1[], unsigned long n1, real data2[], unsigned long n2,
	real *f, real *prob)
{
	real var1,var2,ave1,ave2,df1,df2;

	avevar(data1,n1,&ave1,&var1);
	avevar(data2,n2,&ave2,&var2);
	if (var1 > var2) {
		*f=var1/var2;
		df1=n1-1;
		df2=n2-1;
	} else {
		*f=var2/var1;
		df1=n2-1;
		df2=n1-1;
	}
	*prob = 2.0*betai(0.5*df2,0.5*df1,df2/(df2+df1*(*f)));
	if (*prob > 1.0) *prob=2.0-*prob;
}

void tptest(real data1[], real data2[], unsigned long n, real *t,
	real *prob)
{
	unsigned long j;
	real var1,var2,ave1,ave2,sd,df,cov=0.0;

	avevar(data1,n,&ave1,&var1);
	avevar(data2,n,&ave2,&var2);
	for (j=1;j<=n;j++)
		cov += (data1[j]-ave1)*(data2[j]-ave2);
	cov /= df=n-1;
	sd=sqrt((var1+var2-2.0*cov)/n);	
	dprintf(1,"sd = %g\n",sd);
	if (sd > 0.0) {
	  *t=(ave1-ave2)/sd;
	  *prob=betai(0.5*df,0.5,df/(df+(*t)*(*t)));
	} else {
	  *t = 0.0;
	  *prob = 1.0;
	}
}

void ttest(real data1[], unsigned long n1, real data2[], unsigned long n2,
	real *t, real *prob)
{
	real var1,var2,svar,df,ave1,ave2;

	avevar(data1,n1,&ave1,&var1);
	avevar(data2,n2,&ave2,&var2);
	df=n1+n2-2;
	svar=((n1-1)*var1+(n2-1)*var2)/df;
	*t=(ave1-ave2)/sqrt(svar*(1.0/n1+1.0/n2));
	*prob=betai(0.5*df,0.5,df/(df+(*t)*(*t)));
}

void tutest(real data1[], unsigned long n1, real data2[], unsigned long n2,
	real *t, real *prob)
{
	real var1,var2,df,ave1,ave2;

	avevar(data1,n1,&ave1,&var1);
	avevar(data2,n2,&ave2,&var2);
	*t=(ave1-ave2)/sqrt(var1/n1+var2/n2);
	df=sqr(var1/n1+var2/n2)/(sqr(var1/n1)/(n1-1)+sqr(var2/n2)/(n2-1));
	*prob=betai(0.5*df,0.5,df/(df+sqr(*t)));
}


void avevar(real data[], unsigned long n, real *ave, real *var)
{
	unsigned long j;
	real s,ep;

	for (*ave=0.0,j=1;j<=n;j++) *ave += data[j];
	*ave /= n;
	*var=ep=0.0;
	for (j=1;j<=n;j++) {
		s=data[j]-(*ave);
		ep += s;
		*var += s*s;
	}
	*var=(*var-ep*ep/n)/(n-1);
}

