/*
 *  Student T and F test; and KS2 test
 *	8-oct-97	0.2
 *	28-jan-98	0.2a   isolated beta functions in numrec 
 *       5-oct-00       0.3    added KStwo test
 *      23-mar-01       0.3a   more sanity check of input
 *      29-sep-01       0.3c   handle identical tables
 *      23-nov-02       0.4    bug calling NR, was my brain on steroids?
 *
 * TODO?
 *	- xmin/xmax for range selection
 */

#include <stdinc.h>
#include <getparam.h>

string defv[] = {
    "in1=???\n		First dataset",
    "in2=???\n		Second dataset",
    "col1=1\n           Column from 1st dataset",
    "col2=1\n           Column from 2nd dataset",
    "nmax=10000\n       Max lines in data, if pipe",
    "VERSION=0.4\n	23-nov-02 PJT",
    NULL,
};

string usage="KS2, Student T and F test for mean and variance comparisons";

/* local NR routines ; now properly 0 based */

void ftest(real data1[], int  n1, real data2[], int  n2,
	real *f, real *prob);
void tptest(real data1[], real data2[], int  n, real *t,
	real *prob);
void ttest(real data1[], int  n1, real data2[], int  n2,
	real *t, real *prob);
void tutest(real data1[], int  n1, real data2[], int  n2,
	real *t, real *prob);
void avevar(real data[], int  n, real *ave, real *var);

/* NEMO's NR routines in src/kernel/numrec */

extern real betai(real a, real b, real x);

extern void kstwo(real *data1, int n1, real *data2, int n2, real *d, real *p);

nemo_main()
{
    int i, npt, npt1, npt2, nmax = getiparam("nmax");
    int col1 = getiparam("col1");
    int col2 = getiparam("col2");
    real *coldat1[2], *coldat2[2];
    int colnr1[2], colnr2[2];
    real *x1 = NULL, *x2 = NULL;
    real t, f, tu, tp, tprob, fprob, tuprob, tpprob, ks2, ks2prob;
    string input1 = getparam("in1");
    string input2 = getparam("in2");
    stream instr1, instr2;


    instr1 = stropen(input1,"r");
    if (!streq(input1,input2)) {
        dprintf(0,"Two different files\n");
        instr2 = stropen(input2,"r");
    } else {
        dprintf(0,"One file, different columns\n");
        instr2 = NULL;
    }

    npt1 = nemo_file_lines(input1,nmax);
    if (instr2)
      npt2 = nemo_file_lines(input2,nmax);
    else
      npt2 = npt1;
    x1 = (real *) allocate(npt1 * sizeof(real));
    x2 = (real *) allocate(npt2 * sizeof(real));


    coldat1[0] = x1;    colnr1[0] = col1;
    if (instr2) {
        coldat2[0] = x2;    colnr2[0] = col2;
        npt1 = get_atable(instr1,1,colnr1,coldat1,npt1);
        npt2 = get_atable(instr2,1,colnr2,coldat2,npt2);
    } else {
        coldat1[1] = x2;    colnr1[1] = col2;
        npt1 = get_atable(instr1,2,colnr1,coldat1,npt1);
        npt2 = npt1;
    }
    
    dprintf(0,"Npt1=%d Npt2=%d\n",npt1,npt2);
    if (npt1 < 0) {
      warning("Could not read all data from %s",input1);
      npt1 = -npt1;
    }
    if (npt2 < 0) {
      warning("Could not read all data from %s",input2);
      npt2 = -npt2;
    }

    strclose(instr1);
    if (instr2)
        strclose(instr2);

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

    /* arrays are now going to be sorted and changed by kstwo !! */

    kstwo(x1, npt1, x2, npt2, &ks2, &ks2prob);
    printf("KS2-test: ks2=%g  prob=%g\n",ks2,ks2prob);
}


/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */
/* adapted for NEMO					  */


void ftest(real data1[], int  n1, real data2[], int n2,
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

void tptest(real data1[], real data2[], int  n, real *t,
	real *prob)
{
	int  j;
	real var1,var2,ave1,ave2,sd,df,cov=0.0;

	avevar(data1,n,&ave1,&var1);
	avevar(data2,n,&ave2,&var2);
	for (j=0;j<n;j++)
		cov += (data1[j]-ave1)*(data2[j]-ave2);
	cov /= df=n-1;
	sd=sqrt((var1+var2-2.0*cov)/n);
	if (sd > 0.0) {
	  *t=(ave1-ave2)/sd;
	  *prob=betai(0.5*df,0.5,df/(df+(*t)*(*t)));
	} else {
	  *t = 0.0;
	  *prob = 1.0;
	}
}

void ttest(real data1[], int  n1, real data2[], int  n2,
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

void tutest(real data1[], int n1, real data2[], int n2,
	real *t, real *prob)
{
	real var1,var2,df,ave1,ave2;

	avevar(data1,n1,&ave1,&var1);
	avevar(data2,n2,&ave2,&var2);
	*t=(ave1-ave2)/sqrt(var1/n1+var2/n2);
	df=sqr(var1/n1+var2/n2)/(sqr(var1/n1)/(n1-1)+sqr(var2/n2)/(n2-1));
	*prob=betai(0.5*df,0.5,df/(df+sqr(*t)));
}


void avevar(real data[], int  n, real *ave, real *var)
{
	int j;
	real s,ep;

	for (*ave=0.0,j=0;j<n;j++) *ave += data[j];
	*ave /= n;
	*var=ep=0.0;
	for (j=0;j<n;j++) {
		s=data[j]-(*ave);
		ep += s;
		*var += s*s;
	}
	*var=(*var-ep*ep/n)/(n-1);
}








