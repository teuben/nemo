/* NEMO version of nrutil.c */

#include <stdinc.h>
#include <malloc.h>


void nrerror(error_text)
char error_text[];
/* Numerical Recipes standard error handler */
{
	warning("Numerical Recipes run-time error...");
	warning("%s",error_text);
	error("...now exiting to system...");
}

float *fvector(nl,nh)
int nl,nh;
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float))-nl;
	if (!v) nrerror("allocation failure in fvector()");
	return v;
}

int *ivector(nl,nh)
int nl,nh;
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int))-nl;
	if (!v) nrerror("allocation failure in ivector()");
	return v;
}

double *dvector(nl,nh)
int nl,nh;
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double))-nl;
	if (!v) nrerror("allocation failure in dvector()");
	return v;
}

float **fmatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*))-nrl;
	if (!m) nrerror("allocation failure 1 in fmatrix()");

	/* allocate rows and set pointers to them */
	for(i=nrl;i<=nrh;i++) {
		m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float))-ncl;
		if (!m[i]) nrerror("allocation failure 2 in fmatrix()");
	}
	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*))-nrl;
	if (!m) nrerror("allocation failure 1 in dmatrix()");

	/* allocate rows and set pointers to them */
	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double))-ncl;
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
	}
	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
/* allocate an int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i,**m;

	/* allocate pointers to rows */
	m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*))-nrl;
	if (!m) nrerror("allocation failure 1 in imatrix()");

	/* allocate rows and set pointers to them */
	for(i=nrl;i<=nrh;i++) {
		m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int))-ncl;
		if (!m[i]) nrerror("allocation failure 2 in imatrix()");
	}
	/* return pointer to array of pointers to rows */
	return m;
}

float **submatrix(a,oldrl,oldrh,oldcl,oldch,newrl,newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
float **a;
int oldrl,oldrh,oldcl,oldch,newrl,newcl;
{
	int i,j;
	float **m;

	/* allocate array of pointers to rows */
	m=(float **) malloc((unsigned) (oldrh-oldrl+1)*sizeof(float*))-newrl;
	if (!m) nrerror("allocation failure in submatrix()");

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+oldcl-newcl;

	/* return pointer to array of pointers to rows */
	return m;
}

float **convert_matrix(a,nrl,nrh,ncl,nch)
float *a;
int nrl,nrh,ncl,nch;
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	int i,j,nrow,ncol;
	float **m;

	nrow=nrh-nrl+1;
	ncol=nch-ncl+1;
	/* allocate pointers to rows */
	if ((m=(float **) malloc((unsigned) (nrow)*sizeof(float*))) == NULL)
		nrerror("allocation failure in convert_matrix()");
	m -= nrl;

	/* set pointers to rows */
	for(i=0,j=nrl;i<=nrow-1;i++,j++) m[j]=a+ncol*i-ncl;
	/* return pointer to array of pointers to rows */
	return m;
}

void free_fvector(v,nl,nh)
float *v;
int nl,nh;
/* free a float vector allocated with fvector() */
{
	free((char*) (v+nl));
}

void free_ivector(v,nl,nh)
int *v,nl,nh;
/* free an int vector allocated with ivector() */
{
	free((char*) (v+nl));
}

void free_dvector(v,nl,nh)
double *v;
int nl,nh;
/* free a double vector allocated with dvector() */
{
	free((char*) (v+nl));
}

void free_fmatrix(m,nrl,nrh,ncl,nch)
float **m;
int nrl,nrh,ncl,nch;
/* free a float matrix allocated by fmatrix() */
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_dmatrix(m,nrl,nrh,ncl,nch)
double **m;
int nrl,nrh,ncl,nch;
/* free a double matrix allocated by dmatrix() */
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
int nrl,nrh,ncl,nch;
/* free an int matrix allocated by imatrix() */
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_submatrix(b,nrl,nrh,ncl,nch)
/* free a submatrix allocated by submatrix() */
float **b;
int nrl,nrh,ncl,nch;
{
	free((char*) (b+nrl));
}

void free_convert_matrix(b,nrl,nrh,ncl,nch)
/* free a matrix allocated by convert_matrix() */
float **b;
int nrl,nrh,ncl,nch;
{
	free((char*) (b+nrl));
}
