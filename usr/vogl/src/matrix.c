
#include <stdio.h>
#include "vogl.h"

static	Mstack	*msfree = (Mstack *)NULL;

/*
 * copyvector
 *
 * Copy the 4 vector b to a.
 *
 */
void
copyvector(a, b)
	register	Vector	a, b;
{
	a[0] = b[0];
	a[1] = b[1];
	a[2] = b[2];
	a[3] = b[3];
}

/*
 * copymatrix
 *
 * Copy the  4 x 4 matrix b to the 4 x 4 matrix a
 *
 */
void
copymatrix(a, b)
	register	Matrix	a, b;
{
	register int	i;
	register float	*pa, *pb;

	pa = (float *)a;
	pb = (float *)b;
	for(i = 0; i < 16; i++)
		*(pa++) = *(pb++);
}

/*
 * copytranspose
 *
 *	copy the transpose of the 4 x 4 matrix b to the 4 x 4 matrix a.
 */
void
copytranspose(a, b)
	register Matrix	a, b;
{
	register int	i, j;

	for(i = 0; i < 4; i++)
		for(j = 0; j < 4; j++)
			a[i][j] = b[j][i];
}

/*
 * Retreive the top matrix on the stack and place it in m
 */
void
getmatrix(m)
	Matrix m;
{
	copymatrix(m, vdevice.transmat->m);
}

/*
 * pushmatrix
 *
 * Push the top matrix of the stack down, placing a copy
 * of it on the top of the stack.
 *
 */
void
pushmatrix()
{
	Mstack	*tmpmat;
	Token	*p;

	if (vdevice.inobject) {
		p = newtokens(1);

		p->i = PUSHMATRIX;

		return;
	}

	if (msfree != (Mstack *)NULL) {
		tmpmat = vdevice.transmat;
		vdevice.transmat = msfree;
		msfree = msfree->back;
		vdevice.transmat->back = tmpmat;
		copymatrix(vdevice.transmat->m, tmpmat->m);
	} else {
		tmpmat = (Mstack *)vallocate(sizeof(Mstack));
		tmpmat->back = vdevice.transmat;
		copymatrix(tmpmat->m, vdevice.transmat->m);
		vdevice.transmat = tmpmat;
	}
}

/*
 * popmatrix
 *
 * Pop the top matrix from the stack.
 *
 */
void
popmatrix()
{
	Token	*p;
	Mstack	*oldtop;

	if (vdevice.inobject) {
		p = newtokens(1);

		p->i = POPMATRIX;

		return;
	}

	if (vdevice.transmat->back == (Mstack *)NULL)
		verror("popmatrix: matrix stack empty");
	else {
		oldtop = vdevice.transmat;
		vdevice.transmat = vdevice.transmat->back;
		oldtop->back = msfree;
		msfree = oldtop;
	}

	vdevice.cpVvalid = 0;	/* may have changed mapping from world to device coords */
}

/*
 * loadmatrix
 *
 * Replace the top matrix on the stack
 *
 */
void
loadmatrix(mat)
	Matrix	mat;
{
	register int	i;
	register float	*cm, *mp;
	Token		*p;

	if (!vdevice.initialised)
		verror("loadmatrix: vogl not initialised");

	if (vdevice.inobject) {
		p = newtokens(17);

		p[0].i = LOADMATRIX;
		cm = (float *)mat;
		for (i = 0; i < 16; i++)
			(++p)->f = *cm++;

		return;
	}

	cm = (float *)vdevice.transmat->m;
	mp = (float *)mat;
	for (i = 0; i < 16; i++)
		*cm++ = *mp++;

	vdevice.cpVvalid = 0;		/* may have changed mapping from world to device coords */
}

/*
 * mult4x4
 *
 *	multiply 4 x 4 matrices b and c assigning them into a. Readers are
 * reminded that speed can be important here.
 *
 */
void
mult4x4(a, b, c)
	register Matrix	a, b, c;
{
	a[0][0] = b[0][0] * c[0][0] + b[0][1] * c[1][0] + b[0][2] * c[2][0] + b[0][3] * c[3][0];
	a[0][1] = b[0][0] * c[0][1] + b[0][1] * c[1][1] + b[0][2] * c[2][1] + b[0][3] * c[3][1];
	a[0][2] = b[0][0] * c[0][2] + b[0][1] * c[1][2] + b[0][2] * c[2][2] + b[0][3] * c[3][2];
	a[0][3] = b[0][0] * c[0][3] + b[0][1] * c[1][3] + b[0][2] * c[2][3] + b[0][3] * c[3][3];

	a[1][0] = b[1][0] * c[0][0] + b[1][1] * c[1][0] + b[1][2] * c[2][0] + b[1][3] * c[3][0];
	a[1][1] = b[1][0] * c[0][1] + b[1][1] * c[1][1] + b[1][2] * c[2][1] + b[1][3] * c[3][1];
	a[1][2] = b[1][0] * c[0][2] + b[1][1] * c[1][2] + b[1][2] * c[2][2] + b[1][3] * c[3][2];
	a[1][3] = b[1][0] * c[0][3] + b[1][1] * c[1][3] + b[1][2] * c[2][3] + b[1][3] * c[3][3];

	a[2][0] = b[2][0] * c[0][0] + b[2][1] * c[1][0] + b[2][2] * c[2][0] + b[2][3] * c[3][0];
	a[2][1] = b[2][0] * c[0][1] + b[2][1] * c[1][1] + b[2][2] * c[2][1] + b[2][3] * c[3][1];
	a[2][2] = b[2][0] * c[0][2] + b[2][1] * c[1][2] + b[2][2] * c[2][2] + b[2][3] * c[3][2];
	a[2][3] = b[2][0] * c[0][3] + b[2][1] * c[1][3] + b[2][2] * c[2][3] + b[2][3] * c[3][3];

	a[3][0] = b[3][0] * c[0][0] + b[3][1] * c[1][0] + b[3][2] * c[2][0] + b[3][3] * c[3][0];
	a[3][1] = b[3][0] * c[0][1] + b[3][1] * c[1][1] + b[3][2] * c[2][1] + b[3][3] * c[3][1];
	a[3][2] = b[3][0] * c[0][2] + b[3][1] * c[1][2] + b[3][2] * c[2][2] + b[3][3] * c[3][2];
	a[3][3] = b[3][0] * c[0][3] + b[3][1] * c[1][3] + b[3][2] * c[2][3] + b[3][3] * c[3][3];
}

/*
 * multmatrix
 *
 * Premultipy the top matrix on the stack by "mat"
 *
 */
void
multmatrix(mat)
	Matrix	mat;
{
	Matrix	prod;
	float	*m;
	Token	*p;
	int	i;

	if (vdevice.inobject) {
		p = newtokens(17);

		p[0].i = MULTMATRIX;
		m = (float *)mat;
		for (i = 0; i < 16; i++)
			(++p)->f = *m++;

		return;
	}

	mult4x4(prod, mat, vdevice.transmat->m);
	loadmatrix(prod);
}

/*
 * identmatrix
 *
 * Return a 4 x 4 identity matrix
 *
 */
void
identmatrix(a)
	Matrix 	 a;
{
	register float	*p;

	for (p = (float *)a; p != (float *)a + 16; p++)
		*p = 0;

	a[0][0] = a[1][1] = a[2][2] = a[3][3] = 1;
}

/*
 * multvector
 *
 * Multiply the vector a and the matrix b to form v. Need it to be snappy again.
 * 
 */
void
multvector(v, a, b)
	register	Vector	v, a;
	register	Matrix	b;
{
	v[0] = a[0] * b[0][0] + a[1] * b[1][0] + a[2] * b[2][0] + a[3] * b[3][0];
	v[1] = a[0] * b[0][1] + a[1] * b[1][1] + a[2] * b[2][1] + a[3] * b[3][1];
	v[2] = a[0] * b[0][2] + a[1] * b[1][2] + a[2] * b[2][2] + a[3] * b[3][2];
	v[3] = a[0] * b[0][3] + a[1] * b[1][3] + a[2] * b[2][3] + a[3] * b[3][3];
}

/*
 * premultvector
 *
 * PreMultiply the vector a and the matrix b to form v. 
 * Need it to be snappy again.
 * 
 */
void
premultvector(v, a, b)
	Vector	v, a;
	Matrix	b;
{
	v[0] = a[0] * b[0][0] + a[1] * b[0][1] + a[2] * b[0][2] + a[3] * b[0][3];
	v[1] = a[0] * b[1][0] + a[1] * b[1][1] + a[2] * b[1][2] + a[3] * b[1][3];
	v[2] = a[0] * b[2][0] + a[1] * b[2][1] + a[2] * b[2][2] + a[3] * b[2][3];
	v[3] = a[0] * b[3][0] + a[1] * b[3][1] + a[2] * b[3][2] + a[3] * b[3][3];
}

#ifdef DEBUG 

/*
 * printmat
 *
 *	print s and then dump matrix m. Useful for debugging, you get
 * sick of typing in the print loop otherwise.
 */
printmat(s, m)
	char	*s;
	Matrix	m;
{
	int	i, j;

	printf("%s\n", s);
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++)
			printf("%f ",m[i][j]);
		printf("\n");
	}
}
printvect(s, v)
	char	*s;
	Vector v;
{
	printf("%s %f %f %f %f\n", s, v[0], v[1], v[2], v[3]);
}

#endif
