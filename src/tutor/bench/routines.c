

#define for_all(i)	register int i; for (i = 0; i < n; i++)

int itor(int n)	{for_all(i) v1[i]  = iv2[i];}
int rtoi(int n)	{for_all(i) iv1[i] = v2[i];}
int iadd(int n)	{for_all(i) iv1[i] = iv2[i] + iv3[i];}

int vmove(int n)	{for_all(i) v1[i] = v2[i];}

int ssum1(int n)	{for_all(i) s1 += v1[i];}
int ssum2(int n)	{for_all(i) s1 += v1[i] + v2[i];}
int ssum3(int n)	{for_all(i) s1 += v1[i] * v2[i];}

int vsadd1(int n)	{for_all(i) v1[i] += s1;}
int vsmul1(int n)	{for_all(i) v1[i] *= s1;}
int vsdiv1(int n)	{for_all(i) v1[i] /= s2;}

int vsmul1a(int n)
{
    register int i;
    for (i = 0; i < n; i += 4) {
        v1[i] *= s1;
        v1[i+1] *= s1;
        v1[i+2] *= s1;
        v1[i+3] *= s1;
    }
}

int vsadd2(int n)	{for_all(i) v1[i] = v2[i] + s1;}
int vsmul2(int n)	{for_all(i) v1[i] = v2[i] * s1;}
int vsdiv2(int n)	{for_all(i) v1[i] = v2[i] / s2;}

int vsum1(int n)	{for_all(i) v1[i] += v2[i];}
int vsum2(int n)	{for_all(i) v1[i] = v2[i] + v3[i];}
int vmul1(int n)	{for_all(i) v1[i] *= v2[i];}
int vmul2(int n)	{for_all(i) v1[i] = v2[i] * v3[i];}
int vdiv1(int n)	{for_all(i) v1[i] /= v2[i];}
int vdiv2(int n)	{for_all(i) v1[i] = v3[i] / v2[i];}

int saxpy1(int n)	{for_all(i) v1[i] = s1*v2[i] + s2;}
int saxpy2(int n)	{for_all(i) v1[i] = s1*v2[i] + v3[i];}
int saxpy3(int n)	{for_all(i) v1[i] = v2[i]*v3[i] + v4[i];}

int vsqrt(int n)	{for_all(i) v1[i] = sqrt(v2[i]);}
int vabs(int n)	{for_all(i) v1[i] = fabs(v3[i]);}
int vsin(int n)	{for_all(i) v1[i] = sin(v3[i]);}
int vexp(int n)	{for_all(i) v1[i] = exp(v2[i]);}
int vpow(int n)	{for_all(i) v1[i] = pow(v2[i], s2);}

int scatter1(int n)	{for_all(i) v1[ind1[i]] = v2[i];}
int scatter2(int n)	{for_all(i) v1[ind2[i]] = v2[i];}
int gather1(int n)	{for_all(i) v1[i] = v2[ind1[i]];}
int gather2(int n)	{for_all(i) v1[i] = v2[ind2[i]];}

int vif(int n)		{for_all(i) if (v4[i] > 0.0) v1[i] = v4[i];}

int vforce(int n)
{
    real d0, d1, d2, rij2, fac;

    for_all(i) {
	d0 = v2[i] - x;
	d1 = v3[i] - y;
	d2 = v4[i] - z;
	rij2 = d0*d0 + d1*d1 + d2*d2 + 0.1;
	fac = v1[i] / (rij2*sqrt(rij2));
	v1[0] += d0*fac;
	v1[1] += d1*fac;
	v1[2] += d2*fac;
    }
}

