/*
 *  TABGSL:   some simple interactions with GSL
 */

#include <nemo.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_sf.h>

string defv[] = {
    "x=1.0\n             Evaluate functions",
    "astro=t\n           show ASTRO constants",
    "VERSION=0.1\n       12-May-2022 XYZ",
    NULL,
};

string usage="show GSL constants";

extern double bessi0(double);   // in misc/besselfunc.c

void nemo_main()
{
    bool Qastro = getbparam("astro");
    double x = getdparam("x");

    warning("Just some GSL testing");

    if (Qastro) {
      printf("c    = %f km / s\n", GSL_CONST_MKSA_SPEED_OF_LIGHT);
      printf("AU   = %f km\n", GSL_CONST_MKSA_ASTRONOMICAL_UNIT);
      printf("G    = %f m^3 / kg s^2 \n", GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT);
    }

    double y = gsl_sf_bessel_I0(x);
    double z = bessi0(x);
    printf("bessel_I0 %f %f %f\n",x,y,z);

}

