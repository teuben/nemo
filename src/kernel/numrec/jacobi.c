#include <math.h>

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

extern float *fvector();
extern void nrerror(),free_fvector();


void jacobi(a,n,d,v,nrot)
float **a,d[],**v;
int n,*nrot;
/*209101*/
{
	int j,iq,ip,i;
	float tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

	b=fvector(1,n);
	z=fvector(1,n);
	for (ip=1;ip<=n;ip++) {
/*509102 2 5 0 0*/
		for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=1;ip<=n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
/*509103 2 0 0 0*/
		z[ip]=0.0;
/*509104 2 3 0 0*/
	}
	*nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=1;ip<=n-1;ip++) {
/*509105 2 3 4 0*/
			for (iq=ip+1;iq<=n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
/*509106 2 3 2 0*/
			free_fvector(z,1,n);
			free_fvector(b,1,n);
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
/*509107 2 0 0 0*/
		else
			tresh=0.0;
/*509108 2 3 4 4*/
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++) {
				g=100.0*fabs(a[ip][iq]);
/*309109*/
				if (i > 4 && fabs(d[ip])+g == fabs(d[ip])
					&& fabs(d[iq])+g == fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if (fabs(h)+g == fabs(h))
						t=(a[ip][iq])/h;
/*509110 4 1 6 0*/
					else {
						theta=0.5*h/(a[ip][iq]);
/*509111 6 5 7 6*/
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
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
					for (j=1;j<=ip-1;j++) {
/*509112 6 3 3 6*/
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
/*509113 6 3 3 6*/
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n;j++) {
/*509114 6 3 3 6*/
						ROTATE(a,ip,j,iq,j)
					}
					for (j=1;j<=n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
/*509115 2 0 0 0*/
			z[ip]=0.0;
/*509116 2 0 8 0*/
		}
	}
	nrerror("Too many iterations in routine JACOBI");
}

#undef ROTATE
