/*
 *  Standard error functions and some related to that : (accuracy ~10^-7)
 *  Note that libm.a has erf() and erfc() these days.
 *
 *  Public routines:
 *	erf(x)		error function
 *	erfc(x)		complementary error function
 *	fdawson(x)	dawson's function
 *
 *	dark-ages   Created             PJT
 *      sep-90      f -> fdawson        PJT
 *	feb-95      ansi & prototypes   pjt
 */

#include <math.h>

static 
double  sqrtpi = 1.772453851,
	p  =  0.3275911,
	a1 =  0.254829592,
	a2 = -0.284496736,
	a3 =  1.421413741,
	a4 = -1.453152027,
	a5 =  1.061405429;

/*

    Calculate error function,

                           2     /x        2
               erf (x) = -----  (   exp (-t ) dt
                           1/2   )
                         pi     /0

    from Abramowitz and Stegun, Handbook of Math. Functions, N.B.S. (1964) 299.
                                ---------------------------
    Approximated by form. 7.1.26 from Abramowitz & Stegun (acc. 1.5 10^(-7)

 */

 
double erf(double x)
{
	double t;

	t = 1.0 / (1.0 + p*x);
	return 1.0 - ( t*(a1+t*(a2+t*(a3+t*(a4+t*a5)))) ) * exp(-x*x);
}

double erfc(double x)
{
	return 1.0-erf(x);
}

double fdawson(double x)
{
	double e1, e2;

	e1 = exp(-x*x);
	e2 = exp(-1.0/(x*x));
	return  e1*(x*e2 - sqrtpi*erfc(x));
}

#ifdef TESTBED

#include <stdinc.h>

double Erf(double x )     /* From JH - see cora code */
{
    int i;
    double z, zn, den, den16;
    static double a[6] = {.0705230784, .0422820123, .0092705272, 
                          .0001520143, .0002765672, .0000430638 };

    z = fabs( x );
    if( z > 1.0E+5 )
        return 1.0;
    else {
        zn = 1.0;
        den = 1.0;
        for( i=0; i<6; i++ ) {
            zn = zn*z;
            den = den + a[i]*zn;
        }
        if( den > 100.0 )
            return 1.0;
        else {
            den16 = den*den*den*den;
            den16 = den16*den16*den16*den16;
            return  1.0 - 1.0/den16;
        }
    }
}

/*******************************************************************************    inverf

    Inverse error function.  (Hastings Approx. for Digital Computers)
    Good to about 4 parts in 10**4.
*******************************************************************************/
double inverf( double x )
{
    double eta, num, den, q, sgn;

    if( x > 0.0 )
        sgn = 1.0;
    else
        sgn = -1.0;
    q = 1.0 - fabs( x );
    if( q == 0.0 ) return( 100.0 );
    eta = sqrt( log( 4.0 / (q*q) ) );
    num = (0.010328*eta + 0.802853)*eta + 2.515517;
    den = ((0.001308*eta + 0.189269)*eta + 1.432788)*eta + 1.0;
    return( sgn * 0.70710678 * (eta - num/den) );
}




main()
{
	double x, y, y1, z;
	
	printf ("TESTBED for erf(x) and Dawson's F(x)\n");
	do {
		printf ("Enter x: "); scanf ("%lf",&x);
		if (x<0) break;
		y = erf(x);
                y1 = Erf(x);
		z = fdawson(x);
		printf ("x = %f  erf(x) = %f   Erf(x) = %f Dawson F(x)= %f\n",
                        x,y,y1,z);
	} while (1);
}
#endif


	
	 

