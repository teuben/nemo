/*
 *  MODE: determines the mode of a continuous distribution of data
 *	  For a proper description see NumRec pp462
 *
 *	Peter Teuben    29-nov-86
 */

mode (x, n, j, xmax, sigma)
double *x[];			/* pointer (!) array must contain N data  */
int    n,j;			/* J is the window size                   */
double *xmax;
double *sigma;
{
	int i;
	double f, p, pmax;
	double sqrt();
	void   sort();
			
	if (j<=0)
		j=3;		/* this is really the minimum value       */

	sort (x, n);     	/* If the array is not sorted, do it now  */

	pmax = 0;
	f = (double) j / (double) n;
	for (i=1; i<n-j; i++) {
		p = f / (*x[i+j] - *x[i]);
		if (p>pmax) {
			pmax = p;
			*xmax = 0.5*( *x[i+j] + *x[i]);
		}
	}
	*sigma = pmax / sqrt ( (double)j );
}

#ifdef TESTBED

double  xxx[20]= {0,5,7,8,9,9.5,10,11,13,16,20,25,27,28,28.5,29,30,32,35,40};
main()
{
	double *x[20], xmax, sigma;
	int i,j;
	
	for (i=0; i<20; i++) 
		x[i] = &xxx[i];
	
	/** no sorting otherwise:    sort(x,20);	**/
	
	for (j=3; j<10; j++) {
		mode (x, 20, j, &xmax, &sigma);
		printf ("window %d   xmax=%f  sigma=%f\n",j,xmax,sigma);
	}
}
#endif

